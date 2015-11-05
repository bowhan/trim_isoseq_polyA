#ifndef type_policy_h
#define type_policy_h

#include <string.h>

template <class T>
struct read_policy { };

template <class T>
struct write_policy { };

/* Sequence size calculation */
template<class T>
struct strsize{};

template<>
struct strsize<std::string>
{
	static size_t size(const std::string& s) { return s.size(); }
};

template<size_t N>
struct strsize<char [N]>
{
	static size_t size(const char* s) { return strlen(s) /* necessarilly N-1? */; }
};

template<>
struct strsize<const char *>
{
	static size_t size(const char* s) { return strlen(s); }
};

namespace std
{
	
}


#endif /* type_policy_h */
