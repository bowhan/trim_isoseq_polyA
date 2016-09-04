// Copyright (c) 2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Bo Han

#include <stdio.h>
#include <assert.h>
#include <boost/thread.hpp>
#include "fasta.hpp"
#include "thread.hpp"
#include "polyA_hmm_model.hpp"
#include "kernel_color.h"
#include "common.hpp"
#include "flnc.hpp" // for adjustHeader

using namespace std;

/* send every 100 fasta entries to a thread each time */
const int default_bulk_size = 100;

/* buffer to hold std output in each thread, assuming each fasta (name + seq) is shorter than 50000 */
constexpr int stdout_buffer_size = default_bulk_size * 50000;

/* buffer to hold log in each thread, assuming each fasta (name + polyA length) is shorter than 1000 */
constexpr int stderr_buffer_size = default_bulk_size * 1000;

/* default # of threads */
#define DEFAULT_NUM_THREADS 8

/* mutex for io */
boost::mutex k_io_mx;

using fasta_t = Fasta<caseInsensitiveString>;

void setDefaultHMM(PolyAHmmMode&);

/* thread worker*/
template <class MTQ, bool showColor, bool isoSeqFormat>
class Worker {
    using multi_thread_safe_queue_type = MTQ;
    using container_type = typename multi_thread_safe_queue_type::container_type;
public:
    Worker(const PolyAHmmMode& hmm, multi_thread_safe_queue_type& producer)
        : hmm_(hmm), producer_(producer) {}

    Worker(const Worker& other) = delete;

    Worker(Worker&& other)
        : hmm_(other.hmm_), producer_(other.producer_) {}

    Worker& operator=(const Worker&) = delete;

    void operator()() {
        auto data = producer_.get();
        char stdout_buf[stdout_buffer_size];
        char stderr_buf[stderr_buffer_size];
        size_t stdout_buff_off{0}, stderr_buff_off{0};
        while (!data.empty()) {
            size_t polyalen;
            for (auto& fa : data) {
                const Matrix<int>& path = hmm_.calculateVirtabi(fa.seq_.rbegin(), fa.seq_.size());
                for (polyalen = 0; polyalen < path.size(); ++polyalen) {
                    /* cannot use binary search because polyA might appear in the middle */
                    if (path[polyalen] == PolyAHmmMode::States::NONPOLYA) {
                        break;
                    }
                }
                if (isoSeqFormat) { // static decision
                    if (polyalen)
                        adjustHeader(fa.name_, polyalen);
                }
                stderr_buff_off += sprintf(stderr_buf + stderr_buff_off, "%s\t%zu\n", fa.name_.c_str(), polyalen);
                if (stderr_buff_off * 5 > stderr_buffer_size * 4) {
                    /* manually flush stderr */
                    boost::lock_guard<boost::mutex> lock(k_io_mx);
                    fwrite(stderr_buf, 1, stderr_buff_off, stderr);
                    stderr_buff_off = 0;
                }
                if (showColor) { // static decision; always print
                    stdout_buff_off +=
                        sprintf(stdout_buf + stdout_buff_off
                                , ">%s\n%s" KERNAL_RED "%s \n" KERNAL_RESET
                                , fa.name_.c_str()
                                , fa.seq_.substr(0, fa.seq_.size() - polyalen).c_str()
                                , fa.seq_.substr(fa.seq_.size() - polyalen).c_str()
                        );
                }
                if (!showColor) { // static decision
                    if (polyalen < fa.size()) { // print only when there are at least some non-polyA region
                        stdout_buff_off += sprintf(stdout_buf + stdout_buff_off, ">%s\n%s\n"
                                                   , fa.name_.c_str()
                                                   , fa.seq_.substr(0, fa.seq_.size() - polyalen).c_str()
                        );
                    }
                }
                if (stdout_buff_off * 5 > stdout_buffer_size * 4) {
                    /* manually flush stdout */
                    boost::lock_guard<boost::mutex> lock(k_io_mx);
                    fwrite(stdout_buf, 1, stdout_buff_off, stdout);
                    stdout_buff_off = 0;
                }
            } /* end of for loop to process each fasta in data */
            {
                /* flush buffer */
                boost::lock_guard<boost::mutex> lock(k_io_mx);
                fwrite(stdout_buf, 1, stdout_buff_off, stdout);
                fwrite(stderr_buf, 1, stderr_buff_off, stderr);
            }
            stdout_buff_off = stderr_buff_off = 0;
            data = producer_.get(); /* get new chulk of data */
        }
    }

private:
    PolyAHmmMode hmm_;
    /* keep a COPY of the HMM model since it does mutable calculation inside the class */
    multi_thread_safe_queue_type& producer_;
};

void
ArgumentParse(int argc
              , char** argv
              , string& input
              , string& model_file
              , int& num_of_threads
              , bool& show_color
              , bool& generic_format
             ) {
    string usage =
        KERNAL_GREEN
            "This program trim polyA from Iso-Seq flnc fasta\n"
            "\nusage: [options] flnc.fa\n\n"
            "options:\n"
            KERNAL_CYAN
            "\n[optional]\n"
            "\t-t      number of threads to use, default: 8 \n"
            "\t-m      HMM model file to use, default: None (use default HMM model)\n"
            "\t-c      Instead of trimming polyA, show them in red, default: Off\n"
            "\t-G      Generic fasta file, will NOT adjust the header, default: Off\n"
            KERNAL_RESET;
    int c;
    while ((c = getopt(argc, argv, "Gt:m:ch")) != -1) {
        switch (c) {
            case 'G':generic_format = true;
                break;
            case 't':num_of_threads = atoi(optarg);
                break;
            case 'm':model_file = optarg;
                break;
            case 'c':show_color = true;
                break;
            case 'h':
            default:cerr << usage;
                exit(EXIT_FAILURE);
        }
    }
    if (optind == argc) {
        Error("Please provide input fasta files");
    }
    input = argv[optind];
}


int main(int argc, char** argv) {
    /** options **/
    std::string input_fa_file;
    std::string model_file = "";
    int num_thread = DEFAULT_NUM_THREADS;
    bool show_color = false;
    bool generic_format = false;
    ArgumentParse(argc, argv, input_fa_file, model_file, num_thread, show_color, generic_format);
    if (input_fa_file.empty()) {
        Error("Please provide input fasta file");
        exit(EXIT_FAILURE);
    }
    PolyAHmmMode hmm;
    // initializing HMM model
    if (!model_file.empty()) {
        // read saved HMM model from file
        if (!hmm.read(model_file)) {
            return EXIT_FAILURE;
        }
    } else {
        // using default HMM model
        setDefaultHMM(hmm);
    }
    // trim
    FastaReader<> reader(input_fa_file);
    MultiThreadSafeQueue<fasta_t, std::vector> producer(reader.begin(), reader.end(), default_bulk_size);
    std::vector<boost::thread> threads;
    if (show_color) {
        if (generic_format) /* generic fasta */
            for (int i = 0; i < num_thread; ++i)
                threads.emplace_back(Worker<decltype(producer), true, false>(hmm, producer));
        else /* Iso-Seq FLNC specific fasta, need to adjust some coordinates in the header */
            for (int i = 0; i < num_thread; ++i)
                threads.emplace_back(Worker<decltype(producer), true, true>(hmm, producer));
    } else { // don't show color
        if (generic_format) /* generic fasta */
            for (int i = 0; i < num_thread; ++i)
                threads.emplace_back(Worker<decltype(producer), false, false>(hmm, producer));
        else /* Iso-Seq FLNC specific fasta, need to adjust some coordinates in the header */
            for (int i = 0; i < num_thread; ++i)
                threads.emplace_back(Worker<decltype(producer), false, true>(hmm, producer));
    }
    for (auto& t : threads)
        if (t.joinable())
            t.join();
    return EXIT_SUCCESS;
}

void setDefaultHMM(PolyAHmmMode& hmm) {
    hmm.initialProb(PolyAHmmMode::States::POLYA) = 0.99283668;
    hmm.initialProb(PolyAHmmMode::States::NONPOLYA) = 0.00716332;
    hmm.transProb(PolyAHmmMode::States::POLYA, PolyAHmmMode::States::NONPOLYA) = 3.16493e-07;
    hmm.transProb(PolyAHmmMode::States::POLYA, PolyAHmmMode::States::POLYA) =
        1.0 - hmm.transProb(PolyAHmmMode::States::POLYA, PolyAHmmMode::States::NONPOLYA);
    hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::POLYA) = 2.74842e-09;
    hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::NONPOLYA) =
        1.0 - hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::POLYA);
    hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['A']) = 0.928165;
    hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['C']) = 0.025917;
    hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['G']) = 0.024170;
    hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['T']) = 0.021748;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['A']) = 0.271806;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['C']) = 0.249539;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['G']) = 0.281787;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['T']) = 0.196867;
}
