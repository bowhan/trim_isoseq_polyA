#ifndef sequence_h
#define sequence_h

#include "char_traits.hpp"

using caseInsensitiveString = std::basic_string<char, CaseInsensitiveCharTrait<char>, std::allocator<char> > ;

template <class T = caseInsensitiveString>
struct Sequence
{
// type
    using seq_type = T;
// methods
    Sequence() : seq_{""} { }
    Sequence(const seq_type& other) : seq_(other) { }
    //    template<class Iter> Sequence(Iter a, Iter b) : seq_(a, b) { }
    Sequence(const Sequence<T>& other) : seq_(other.seq_) { }
    Sequence(Sequence<T>&& other) : seq_(std::move(other.seq_)) { }
    Sequence& operator=(const Sequence& other)
    {
        assert(this != &other);
        seq_ = other.seq_;
        return *this;
    }
    Sequence& operator=(Sequence&& other)
    {
        if(this != &other)
        {
            seq_ = std::move(other.seq_);
        }
        return *this;
    }
    
    Sequence& reverse();
    Sequence& complement();
    Sequence& reverse_complement();

    Sequence reverse_copy() const;
    Sequence complement_copy() const;
    Sequence reverse_complement_copy() const;
    size_t size() const { return seq_.size(); }
// data
    seq_type seq_;
};


template <class T>
Sequence<T>& Sequence<T>::reverse()
{
    auto begin    = seq_.begin();
    auto end      = seq_.end();
    auto distance = std::distance(begin, end);
    std::advance(end, -1);
    while (distance > 0) {
        std::swap(*begin, *end);
        std::advance(begin, 1);
        std::advance(end, -1);
        distance -= 2;
    }
    return *this;
}

template <class T>
Sequence<T>& Sequence<T>::complement()
{
    for (auto& x : seq_)
    {
        switch (x)
        {
            case 'A': case 'a': x = 'T'; break;
            case 'C': case 'c': x = 'G'; break;
            case 'G': case 'g': x = 'C'; break;
            case 'T': case 't':
            case 'U': case 'u': x = 'A'; break;
            case 'N': case 'n': x = 'N'; break;
            default: break;
        }
    }
    return *this;
}

template <class T>
Sequence<T>& Sequence<T>::reverse_complement()
{
    this->reverse().complement();
    return *this;
}

template <class T>
Sequence<T> Sequence<T>::reverse_copy() const
{
    return Sequence<T>{ seq_type{ seq_.crbegin(), seq_.crend() } };
}

template <class T>
Sequence<T> Sequence<T>::complement_copy() const
{
    Sequence<T> ret{ *this };
    ret.complement();
    return ret;
}

template <class T>
Sequence<T> Sequence<T>::reverse_complement_copy() const
{
    Sequence<T> ret{ seq_type{ seq_.crbegin(), seq_.crend() } };
    ret.complement();
    return ret;
}

template<class T, class U>
bool operator==(const Sequence<T>& lhs, const Sequence<U>& rhs)
{
    return lhs.seq_ == rhs.seq_;
}

#endif /* sequence_h */
