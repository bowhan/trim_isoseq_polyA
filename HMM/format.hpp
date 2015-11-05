#ifndef format_h
#define format_h

#include <unistd.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "type_policy.h"

template <class T> class FormatReader;

/* define iterator */
template<class T>
class FormatReaderIter : public boost::iterator_facade<FormatReaderIter<T>, T, std::input_iterator_tag>
{
public:
    // default constructr, meant to be called by FormatReader.end();
    FormatReaderIter() : fp_(nullptr), data_{ new T{} } { }
    
    explicit FormatReaderIter(boost::iostreams::filtering_istream* fp) : fp_(fp)
    {
        this->increment();
    }
    
private:
    friend class boost::iterator_core_access;
    
    void increment()
    {
        data_.reset(new T{ read_policy<T>::read(fp_) });
        // the policy should return default T{} to indicate EOF or ill-formated file
    }
    
    bool equal(const FormatReaderIter& other) const
    {
        return *data_ == *(other.data_); // only called by != ...end()
    }
    
    T& dereference() const
    {
        return *data_;
    }
    
    boost::iostreams::filtering_istream* fp_ = nullptr;
    std::shared_ptr<T> data_ = nullptr;
}; /* end of iterator definition */

template <class T>
class FormatReader
{
public:
    using iterator = FormatReaderIter<T>;
    using const_iterator = const iterator;
    
    explicit FormatReader(const std::string& file_name)
    {
        std::istream* p_ist_in{ &std::cin };
        if (file_name != "stdin" && file_name != "-")
        {
            if(access(file_name.c_str(), R_OK) != 0)
            {
                fprintf(stderr, "error, cannot read file %s. Please double check.\n", file_name.c_str());
                exit(EXIT_FAILURE);
            }
            p_ist_in = new std::ifstream{ file_name };
        }
        else
        {
            std::ios::sync_with_stdio(false);
            std::cin.tie(nullptr);
            std::cerr.tie(nullptr);
        }
        
        char magic_number[4] = "\0\0\0";
        p_ist_in->get(magic_number, 3);
        
        if (magic_number[0] == '\037' && magic_number[1] == (char)'\213')
        {
            ins_.push(boost::iostreams::gzip_decompressor());
        }
        else if (magic_number[0] == 'B' && magic_number[1] == 'Z')
        {
            ins_.push(boost::iostreams::bzip2_decompressor());
        }
        else
        {
            
        }
        p_ist_in->seekg(0, p_ist_in->beg);
        ins_.push(*p_ist_in);
    }
    
    ~FormatReader() {}
    
    iterator begin()
    {
        return iterator{ &ins_ };
    }
    
    iterator end()
    {
        return iterator{};
    }
    
    const_iterator cbegin()
    {
        return iterator{ &ins_ };
    }
    
    const_iterator cend()
    {
        return iterator{};
    }

protected:
    boost::iostreams::filtering_istream ins_;
};


#endif /* format_h */
