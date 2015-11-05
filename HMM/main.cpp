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
#include "fasta.hpp"
#include "polyA_hmm_model.hpp"

#define MYDEBUG

void usage(const char*);
void setDefaultHMM(PolyAHmmMode&);

int main(int argc, const char * argv[])
{
	if(argc < 2)
	{
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	PolyAHmmMode hmm;
	setDefaultHMM(hmm);
	
	FastaReader<> fa_reader{argv[1]};
	for(auto& fa : fa_reader)
	{
		const Matrix<int>& path = hmm.calculateVirtabi(fa.seq_.rbegin(), fa.seq_.size());
		size_t polyalen = 0;
		for(polyalen = 0; polyalen < path.size(); ++polyalen)
		{
			if(path[polyalen] == PolyAHmmMode::States::NONPOLYA)
			{
				break;
			}
		}
#ifdef MYDEBUG
		fprintf(stderr, "%s\t%zu\n", fa.name_.c_str(), polyalen);
#endif
		if(polyalen < fa.size())
		{
			fprintf(stdout, ">%s\n%s\n",
					fa.name_.c_str(),
					fa.seq_.substr(0, fa.seq_.size()-polyalen).c_str());
		}
	}
    return EXIT_SUCCESS;
}

void usage(const char* p)
{
	
	fprintf(stderr,
			"this program trims polyA from Iso-Seq classify output"
			"usage:  %s  to_be_trimmed.fa \n", p);
}

void setDefaultHMM(PolyAHmmMode& hmm)
{
	hmm.initialProb(PolyAHmmMode::States::POLYA)    = 0.5;
	hmm.initialProb(PolyAHmmMode::States::NONPOLYA) = 0.5;
	
	hmm.transProb(PolyAHmmMode::States::POLYA,     PolyAHmmMode::States::POLYA)   = 0.7;
	hmm.transProb(PolyAHmmMode::States::POLYA,    PolyAHmmMode::States::NONPOLYA) = 0.3;
	hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::POLYA)    = 0.0;
	hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::NONPOLYA) = 1.0;
	
    hmm.emitProb(PolyAHmmMode::States::POLYA,    to_idx['A']) = 0.96;
    hmm.emitProb(PolyAHmmMode::States::POLYA,    to_idx['C']) = 0.01;
    hmm.emitProb(PolyAHmmMode::States::POLYA,    to_idx['G']) = 0.01;
    hmm.emitProb(PolyAHmmMode::States::POLYA,    to_idx['T']) = 0.01;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['A']) = 0.3;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['C']) = 0.2;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['G']) = 0.2;
    hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['T']) = 0.3;
}