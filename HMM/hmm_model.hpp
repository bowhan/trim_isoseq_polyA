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
#ifndef hmm_model_hpp
#define hmm_model_hpp

#include "matrix.hpp"
#include <fstream>

class HmmModeBase {
    // type
protected:
    using value_type    = double;
    using matrix_type   = Matrix<value_type>;
    using pointer       = value_type*;
    using const_pointer = const pointer;
    using reference     = value_type&;

    // methods
public:
    explicit HmmModeBase(int sta = 0, int sym = 0)
        : no_states_(sta)
        , no_symbol_(sym)
        , init_(sta, 1)
        , tran_(sta, sta)
        , emit_(sta, sym)
    {
    }

    virtual ~HmmModeBase() {}

    HmmModeBase(const HmmModeBase&) = delete;

    HmmModeBase(HmmModeBase&& other)
        : no_states_(other.no_states_)
        , no_symbol_(other.no_symbol_)
        , init_(std::move(other.init_))
        , tran_(std::move(other.tran_))
        , emit_(std::move(other.emit_))
    {
    }

    HmmModeBase& operator=(const HmmModeBase&) = delete;

    HmmModeBase& operator=(HmmModeBase&& other)
    {
        if (this != &other) {
            no_states_ = other.no_states_;
            no_symbol_ = other.no_symbol_;
            init_ = std::move(other.init_);
            tran_ = std::move(other.tran_);
            emit_ = std::move(other.emit_);
        }
        return *this;
    }

    size_t states() const { return no_states_; }
    size_t symbols() const { return no_symbol_; }

    reference initialProb(size_t i)
    {
        assert(i < no_states_);
        return init_(i, 0);
    }

    void initialProb(size_t i, value_type v)
    {
        assert(i < no_states_);
        init_(i, 0) = v;
    }

    reference transProb(size_t i, size_t j)
    {
        assert(i < no_states_);
        assert(j < no_states_);
        return tran_(i, j);
    }

    void transProb(size_t i, size_t j, value_type v)
    {
        assert(i < no_states_);
        assert(j < no_states_);
        tran_(i, j) = v;
    }

    reference emitProb(size_t i, size_t j)
    {
        assert(i < no_states_);
        assert(j < no_symbol_);
        return emit_(i, j);
    }

    void emitProb(size_t i, size_t j, value_type v)
    {
        assert(i < no_states_);
        assert(j < no_symbol_);
        emit_(i, j) = v;
    }

    virtual bool read(const std::string& filename)
    { /* serialization is overkill */
        std::ifstream ifs;
        try {
            ifs.open(filename);
            if (!ifs)
                return false;
            ifs >> no_states_ >> no_symbol_;
            // reading init_
            init_.reSize(no_states_, 1);
            for (size_t i = 0; i < no_states_; ++i) {
                ifs >> init_(i, 0);
            }
            // reading trans_
            tran_.reSize(no_states_, no_states_);
            for (size_t i = 0; i < no_states_; ++i) {
                for (size_t j = 0; j < no_states_; ++j) {
                    ifs >> tran_(i, j);
                }
            }
            // reading emit_
            emit_.reSize(no_states_, no_symbol_);
            for (size_t i = 0; i < no_states_; ++i) {
                for (size_t j = 0; j < no_symbol_; ++j) {
                    ifs >> emit_(i, j);
                }
            }
        }
        catch (std::ifstream::failure) {
            return false;
        }
        return true;
    }

    virtual bool write(const std::string& filename)
    {
        std::ofstream ofs;
        try {
            ofs.open(filename);
            if (!ofs)
                return false;
            ofs << no_states_ << ' ' << no_symbol_ << ' ';
            // reading init_
            for (size_t i = 0; i < no_states_; ++i) {
                ofs << init_(i, 0) << ' ';
            }
            // reading trans_
            for (size_t i = 0; i < no_states_; ++i) {
                for (size_t j = 0; j < no_states_; ++j) {
                    ofs << tran_(i, j) << ' ';
                }
            }
            // reading emit_
            for (size_t i = 0; i < no_states_; ++i) {
                for (size_t j = 0; j < no_symbol_; ++j) {
                    ofs << emit_(i, j) << ' ';
                }
            }
        }
        catch (std::ofstream::failure) {
            return false;
        }
        return true;
    }

    // data
protected:
    size_t no_states_;
    size_t no_symbol_;
    matrix_type init_;
    matrix_type tran_;
    matrix_type emit_;
};

#endif
