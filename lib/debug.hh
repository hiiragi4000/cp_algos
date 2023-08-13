#ifndef DEBUG_HH
#define DEBUG_HH

#include<deque>
#include<forward_list>
#include<list>
#include<map>
#include<ostream>
#include<set>
#include<tuple>
#include<unordered_map>
#include<unordered_set>
#include<utility>
#include<vector>

template<typename T, typename U>
std::ostream &operator<<(std::ostream &os, std::pair<T, U> const &p){
   return os << '(' << p.first << ", " << p.second << ')';
}

#define PRINT_ELEMENTS(bracket, delimiter){\
   os << bracket[0];\
   bool first = true;\
   for(auto const &ai: a){\
      if(first){\
         first = false;\
      }else os << delimiter;\
      os << ai;\
   }\
   return os << bracket[1];\
}
#define DEF_PRINT_CONTAINER(container, bracket, delimiter) template<typename T>\
std::ostream &operator<<(std::ostream &os, container<T> const &a) PRINT_ELEMENTS(bracket, delimiter)
DEF_PRINT_CONTAINER(std::vector, "[]", ", ")
DEF_PRINT_CONTAINER(std::deque, "[]", ", ")
DEF_PRINT_CONTAINER(std::list, "[]", " <-> ")
DEF_PRINT_CONTAINER(std::forward_list, "[]", " -> ")
DEF_PRINT_CONTAINER(std::unordered_set, "{}", ", ")
DEF_PRINT_CONTAINER(std::unordered_multiset, "{}", ", ")
DEF_PRINT_CONTAINER(std::set, "{}", ", ")
DEF_PRINT_CONTAINER(std::multiset, "{}", ", ")
#undef DEF_PRINT_CONTAINER
template<typename T, size_t N>
std::ostream &operator<<(std::ostream &os, std::array<T, N> const &a) PRINT_ELEMENTS("[]", ", ")
#undef PRINT_ELEMENTS

#define DEF_PRINT_DICT(dict) template<typename K, typename V>\
std::ostream &operator<<(std::ostream &os, dict<K, V> const &a){\
   os << '{';\
   bool first = true;\
   for(auto const &p: a){\
      if(first){\
         first = false;\
      }else os << ", ";\
      os << p.first << ": " << p.second;\
   }\
   return os << '}';\
}
DEF_PRINT_DICT(std::unordered_map)
DEF_PRINT_DICT(std::unordered_multimap)
DEF_PRINT_DICT(std::map)
DEF_PRINT_DICT(std::multimap)
#undef DEF_PRINT_DICT

namespace impl{
template<size_t I, typename ...Ts> struct TupleFormatter{
   std::ostream &operator()(std::ostream &os, std::tuple<Ts...> const &t){
      if constexpr(I >= sizeof...(Ts)){
         return os;
      }else{
         if constexpr(I > 0){
            os << ", ";
         }
         os << std::get<I>(t);
         return TupleFormatter<I+1, Ts...>()(os, t);
      }
   }
};
} // namespace impl

template<typename ...Ts>
std::ostream &operator<<(std::ostream &os, std::tuple<Ts...> const &t){
   return impl::TupleFormatter<0, Ts...>()(os<<'(', t) << ')';
}

#endif
