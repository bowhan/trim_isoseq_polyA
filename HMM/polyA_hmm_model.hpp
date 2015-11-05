#ifndef polyA_hmm_model_hpp
#define polyA_hmm_model_hpp

#include <string>
#include <iostream>
#include <cmath>
#include "sequence.hpp" // policy strsize<>::size()
#include "hmm_model.hpp"
#include "hmm_utilities.h"

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
public:
    enum States { POLYA = 0, NONPOLYA = 1, UNKNOWN = 2 };
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
	
	/* evaluating algorithms */
public:
    template<class TSequence> const matrix_type& calculateForward (const TSequence&);
	template<class TIter>     const matrix_type& calculateForward (TIter, size_t);
    template<class TSequence> const matrix_type& calculateBackward(const TSequence&);
	template<class TIter>     const matrix_type& calculateBackward(TIter, size_t);
	
//	template<class TSequence> const matrix_type& calculatePosterior(const TSequence&, size_t); // TODO: implement
//	template<class TIter>     const matrix_type& calculatePosterior(TIter, size_t, size_t); // TODO: implement
public:
	/* decoding algorithms */
    template<class TSequence> const path_type&   calculateVirtabi (const TSequence&);
	template<class TIter>     const path_type&   calculateVirtabi (TIter, size_t);
	

	/* training algorithms */
	
	//void maximumLikelihoodEstimation();
	//void ViterbiTraining();
	//void BaumWelch();
	
protected:
    matrix_type forw_;
    matrix_type back_;
    matrix_type post_;
    path_type   path_;
};

template<class TIterator>
auto PolyAHmmMode::calculateVirtabi(TIterator striter, size_t N) -> const path_type&
{
	matrix_type prob(no_states_, N);
	prob = 0.0;
	double curmax, tmp;
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
				tmp = prob(k, j-1) + std::log2(tran_(k, i));
				if(tmp > curmax)
				{
					curmax = tmp;
				}
			}
			prob(i, j) = curmax + std::log2(emit_(i, to_idx[*striter]));
		}
	}
	path_.reSize(1, N);
	path_ = States::UNKNOWN;
	// determine the max ending
	curmax = -INFINITY;
	for(int i = 0; i < no_states_; ++i)
	{
		if(prob(i, N-1) > curmax)
		{
			curmax = prob(i, N-1);
			path_[N-1] = i;
		}
	}
	for(int j = int(N-2); j >= 0; --j)
	{
		curmax = -INFINITY;
		for(int i = 0; i < no_states_; ++i)
		{ // from pos j (state: i) to pos j + 1 (state: path_(0, j+1) )
			tmp = prob(i, j) + std::log2(tran_(i, path_[j+1]));
			if(tmp > curmax)
			{
				curmax = tmp;
				path_[j] = i;
			}
		}
	}
//	for(int i = 0; i < no_states_; ++i)
//	{
//		for(int j = 0; j < N; ++j)
//		{
//			fprintf(stderr, "%f\t", prob(i,j));
//		}
//		fprintf(stderr, "\n");
//	}
//	
//	for(int j = 0; j < N; ++j)
//	{
//		fprintf(stderr, "%d\t", path_(0,j));
//	}
//	fprintf(stderr, "\n");
	return path_;
}

template<class TSequence>
auto PolyAHmmMode::calculateVirtabi(const TSequence& seq) -> const path_type&
{
	size_t N = strsize<TSequence>::size(seq);
	return PolyAHmmMode::calculateVirtabi(std::begin(seq), N);
}

template<class TIterator>
auto PolyAHmmMode::calculateForward(TIterator striter, size_t N) -> const matrix_type&
{
	forw_.reSize(no_states_, N);
	forw_ = 0.0;
	double prevmin{1.0}, curmin{1.0};
	// fill first sequence
	for(int i = 0; i < no_states_; ++i)
	{
		forw_(i, 0) = std::log2(init_(i, 0) * emit_(i, to_idx[*striter]));
		if(prevmin > forw_(i, 0)) prevmin = forw_(i, 0); // record the min from previous step
	}
	// dynamically fill d
	++striter; // at seq[1]
	for(int j = 1; j < N; ++j, ++striter)
	{ // j is the observe sequence index
		for(int i = 0; i < no_states_; ++i)
		{ // current state index i
			for(int k = 0; k < no_states_; ++k)
			{ // previous state index
				forw_(i,j) += std::exp2(forw_(k, j-1) - prevmin) * tran_(k, i);
			}
			forw_(i,j) = std::log2(forw_(i,j)) + std::log2(emit_(i, to_idx[*striter])) + prevmin;
			if(curmin > forw_(i,j)) curmin = forw_(i,j);
		}
		prevmin = curmin;
	}
	return forw_;
}

template<class TSequence>
auto PolyAHmmMode::calculateForward(const TSequence& seq) -> const matrix_type&
{
	size_t N = strsize<TSequence>::size(seq);
	return PolyAHmmMode::calculateForward(std::begin(seq), N);
}

template<class TIterator>
auto PolyAHmmMode::calculateBackward(TIterator strriter, size_t N) -> const matrix_type&
{
	back_.reSize(no_states_, N);
	back_ = 0.0;
	double prevmin{1.0}, curmin{1.0};
	// fill first sequence
	for(int i = 0; i < no_states_; ++i)
	{
		back_(i, N-1) = 0.0 /* std::log2(1.) */;
	}
	for(int j = int(N-2); j >= 0; --j, --strriter)
	{ // state at j
		for(int i = 0; i < no_states_; ++i)
		{ // i is on position j
			for(int k = 0; k < no_states_; ++k)
			{ // k is on position j + 1
				back_(i,j) += std::exp2(back_(k, j+1) - prevmin) * tran_(i, k) * emit_(k, to_idx[*strriter]);
			}
			back_(i,j) = std::log2(back_(i,j)) + prevmin;
			if(curmin > back_(i,j)) curmin = back_(i,j);
		}
		prevmin = curmin;
	}
	return back_;
}

template<class TSequence>
auto PolyAHmmMode::calculateBackward(const TSequence& seq) -> const matrix_type&
{
	size_t N = strsize<TSequence>::size(seq);
	// no corespoding std::rbegin()
	// for container, should use rbegin()
	// but for char*, has to overwrite operator++... TODO
	auto iter = std::begin(seq);
	std::advance(iter, N-1);
	return PolyAHmmMode::calculateBackward(iter, N);
}


#endif