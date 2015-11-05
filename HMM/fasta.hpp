#ifndef fasta_h
#define fasta_h

#include <string>
#include <iostream>
#include <cassert>
#include <fstream>
#include <memory>
#include "format.hpp"
#include "sequence.hpp"
#include "type_policy.h"

template<class T = caseInsensitiveString>
struct Fasta : public Sequence<T>
{
    using seq_type = Sequence<T>;
    std::string name_; // sequnece can take different representation such as char*, std::string or twoBits; name is just string
    
};

/* reading policy */
template <class T>
struct read_policy<Fasta<T> >
{
    using string_type = T;
    using fasta_type  = Fasta<T>;
    // reading policy for fasta in the ifstream
    static fasta_type read(std::istream* ins)
    {
        const size_t bufsize = 102400;
        fasta_type fa{};
        if (ins->peek() != '>' || ins->eof())
        {
            return fa; // returning an empty Fa meant EOF or ill-formated file
        }
        ins->get(); // consume '>'
        std::getline(*ins, fa.name_); // read name, which is always std::string
        char buffer[bufsize];
        while (ins->peek() != '>' && ins->good())
        {
            
            ins->getline(buffer, bufsize);
            fa.seq_ += buffer;
        }
        return fa;
    }
};


template<class T> struct FastaSupported
{ const static bool value = false;};

template<> struct FastaSupported<std::string>
{ const static bool value = true; }; // currently string is supported

template<> struct FastaSupported<caseInsensitiveString>
{ const static bool value = true; }; // currently caseInsensitiveString is supported


template <class T = caseInsensitiveString, class = typename std::enable_if<FastaSupported<T>::value>::type>
using FastaReader = FormatReader<Fasta<T> >;

#endif