#ifndef UTIL_HH
#define UTIL_HH

#include<iterator>

template<typename It>
using DerefType = typename std::iterator_traits<It>::value_type;

template<typename T> struct DivT{
   T quot, rem;
};

#endif // UTIL_HH
