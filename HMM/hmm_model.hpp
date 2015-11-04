#ifndef hmm_model_hpp
#define hmm_model_hpp

#include "matrix.hpp"

class HmmModeBase
{
// type
protected:
    using value_type      = double;
    using matrix_type     = Matrix<value_type>;
    using pointer         = value_type*;
    using const_pointer   = const pointer;
    using reference       = value_type&;
    
// methods
public:
    explicit HmmModeBase(int sta = 0, int sym = 0) :
        no_states_(sta),
        no_symbol_(sym),
        init_(sta, 1),
        tran_(sta, sta),
        emit_(sta, sym)
    { }
    
    virtual ~HmmModeBase() { }
    
    HmmModeBase(const HmmModeBase&) = delete;
    
    HmmModeBase(HmmModeBase&& other):
        no_states_(other.no_states_),
        no_symbol_(other.no_symbol_),
        init_(std::move(other.init_)),
        tran_(std::move(other.tran_)),
        emit_(std::move(other.emit_))
    { }
    
    HmmModeBase& operator=(const HmmModeBase&) = delete;
    
    HmmModeBase& operator=(HmmModeBase&& other)
    {
        if(this != &other)
        {
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
    
    // data
protected:
    size_t no_states_;
    size_t no_symbol_;
    matrix_type init_;
    matrix_type tran_;
    matrix_type emit_;
};


#endif
