
#ifndef char_traits_h
#define char_traits_h

#include <string>

template<class TChar = char >
struct CaseInsensitiveCharTrait : public std::char_traits<TChar>
{
    static bool eq(char c1, char c2)
    {
        return toupper(c1) == toupper(c2);
    }
    
    static bool ne(char c1, char c2)
    {
        return toupper(c1) != toupper(c2);
    }
    
    static bool lt(char c1, char c2)
    {
        return toupper(c1) < toupper(c2);
    }
    
    static int compare(const char* s1, const char* s2, size_t n)
    {
        for(size_t i = 0; i < n; ++i)
        {
            if(toupper(*s1) < toupper(*s2)) return -1;
            else if(toupper(*s1) > toupper(*s2)) return 1;
        }
        return 0;
    }
    
    static const char* find(const char* s, int n, char a)
    {
        while(n-- > 0 && toupper(*s) != toupper(a))
        {
            ++s;
        }
        return s; // returning \0 means not found
    }
    
};


#endif /* char_traits_h */
