#ifndef polyA_hmm_model_hpp
#define polyA_hmm_model_hpp

#include <string>
#include <cmath>
#include "hmm_model.hpp"
#include "hmm_utilities.h"

template<class TSequence = std::string>
class PolyAHmmMode : public HmmModeBase
{
    // types
private:
    using _base = HmmModeBase;
    using _self = PolyAHmmMode;
protected:
    using value_type    = typename HmmModeBase::value_type;
    using pointer       = typename HmmModeBase::pointer;
    using matrix_type   = typename HmmModeBase::matrix_type;
    using path_type     = Matrix<int>;
    using const_pointer = typename std::add_const<pointer>::type;
    using seq_t         = TSequence;
public:
    enum States
    {
        POLYA = 0, NONPOLYA = 1
    };
    
    // methods
public:
    PolyAHmmMode() : _base(2, 4) { }
    virtual ~PolyAHmmMode() { }
    PolyAHmmMode(const PolyAHmmMode&) = delete;
    PolyAHmmMode(PolyAHmmMode&& other):
        _base(std::move(other))
        ,forw_(std::move(other.forw_))
        ,back_(std::move(other.back_))
        ,post_(std::move(other.post_))
        ,path_(std::move(other.path_))
    { }
    
    PolyAHmmMode& operator=(const PolyAHmmMode&) = delete;

    
    PolyAHmmMode& operator=(PolyAHmmMode&& other)
    {
        if(this != &other)
        {
            _base::operator=(std::move(other));
            forw_ = std::move(other.forw_);
            back_ = std::move(other.back_);
            post_ = std::move(other.post_);
            path_ = std::move(other.path_);
        }
        return *this;
    }
    
    const matrix_type& calculateForward(const seq_t&);
    const matrix_type& calculateBackward(const seq_t&);
    const path_type&   calculateVirtabi(const seq_t&);
//    const matrix_type& calculatePosterior(const seq_t&); // TODO: implement it
    
protected:
    matrix_type forw_;
    matrix_type back_;
    matrix_type post_;
    path_type   path_;
};

template<class TSequence>
auto PolyAHmmMode<TSequence>::calculateVirtabi(const seq_t& seq) -> const path_type&
{
    size_t N = seq.size();
    matrix_type prob(no_states_, N);
    prob = 0.0;
    path_.reSize(no_states_, N);
    path_ = -1;
    // fill first sequence
    auto striter = seq.cbegin();
    double curmax = 0.0, tmp;
    for(int i = 0; i < no_states_; ++i)
    {
        prob(i, 0) = std::log2(init_(i, 0) * emit_(i, to_idx[*striter]));
    }
    // dynamically fill d
    ++striter; // at seq[1]
    for(int j = 1; j < N; ++j, ++striter)
    { // j is the observe sequence index
        for(int i = 0; i < no_states_; ++i)
        { // current state index i
            curmax = -INFINITY;
            for(int k = 0; k < no_states_; ++k)
            { // previous state index k
                tmp = prob(k, j-1) + std::log2(tran_(k, i)) + std::log2(emit_(i, to_idx[*striter]));
                 if(tmp > curmax)
                 {
                     curmax = tmp;
                     path_(i, j) = k;
                 }
            }
            prob(i, j) = curmax;
        }
    }
    for(int i = 0; i < no_states_; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            fprintf(stderr, "%f\t", prob(i,j));
        }
        fprintf(stderr, "\n");
    }
    for(int i = 0; i < no_states_; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            fprintf(stderr, "%d\t", path_(i,j));
        }
        fprintf(stderr, "\n");
    }
    return path_;
}

template<class TSequence>
auto PolyAHmmMode<TSequence>::calculateForward(const seq_t& seq) -> const matrix_type&
{
    size_t N = seq.size();
    forw_.reSize(no_states_, N);
    forw_ = 0.0;
    // fill first sequence
    auto striter = seq.cbegin();
    for(int i = 0; i < no_states_; ++i)
    {
        forw_(i, 0) = std::log2(init_(i, 0) * emit_(i, to_idx[*striter]));
    }
    
    // dynamically fill d
    ++striter; // at seq[1]
    for(int j = 1; j < N; ++j, ++striter)
    { // j is the observe sequence index
        for(int i = 0; i < no_states_; ++i)
        { // current state index i
            for(int k = 0; k < no_states_; ++k)
            { // previous state index k
                forw_(i,j) += std::exp2(forw_(k, j-1)) * tran_(k, i);
            }
            forw_(i,j) = std::log2(forw_(i,j) * emit_(i, to_idx[*striter]));
        }
    }
    return forw_;
}

template<class TSequence>
auto PolyAHmmMode<TSequence>::calculateBackward(const seq_t& seq) -> const matrix_type&
{
    size_t N = seq.size();
    back_.reSize(no_states_, N);
    back_ = 0.0;
    // fill first sequence
    auto strriter = seq.crbegin();
    for(int i = 0; i < no_states_; ++i)
    {
        back_(i, N-1) = 0.0 /* std::log2(1.) */;
    }
    for(int j = int(N-2); j >= 0; --j, ++strriter)
    { // state at j
        for(int i = 0; i < no_states_; ++i)
        { // i in on j
            for(int k = 0; k < no_states_; ++k)
            { // k is on the j + 1
                back_(i,j) += std::exp2(back_(k, j+1)) * tran_(i, k) * emit_(k, to_idx[*strriter]);
            }
            back_(i,j) = std::log2(back_(i,j));
        }
    }
    return back_;
}


#endif