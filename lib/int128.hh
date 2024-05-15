#ifndef INT128_HH
#define INT128_HH

#include<algorithm>
#include<stdexcept>
#include<string>
#include<cctype>
#ifdef _MSC_BUILD
#include<type_traits>
#include<climits>
#include<intrin.h>
#endif

#define U8 unsigned char
#define I64 long long
#define U64 unsigned long long

#ifdef _MSC_BUILD
#define DEF_TYPE(T)\
   U64 lo, hi;\
   T() = default;\
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>\
   constexpr T(INT a) noexcept: lo(a), hi(-(a<0)){}\
   explicit constexpr operator bool() const noexcept{\
      return lo || hi;\
   }\
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>\
   explicit constexpr operator INT() const noexcept{\
      return static_cast<INT>(lo);\
   }\
   constexpr T operator+() const noexcept{\
      return *this;\
   }\
   constexpr T &into_bitwise_not() noexcept{\
      lo = ~lo; hi = ~hi;\
      return *this;\
   }\
   constexpr T operator~() const noexcept{\
      T res = *this;\
      return res.into_bitwise_not();\
   }\
   constexpr T &operator++() noexcept{\
      hi += (++lo == 0);\
      return *this;\
   }\
   constexpr T operator++(int) noexcept{\
      T res = *this;\
      ++*this;\
      return res;\
   }\
   constexpr T &operator--() noexcept{\
      hi -= (lo-- == 0);\
      return *this;\
   }\
   constexpr T operator--(int) noexcept{\
      T res = *this;\
      --*this;\
      return res;\
   }\
   constexpr T &into_opp() noexcept{\
      return ++into_bitwise_not();\
   }\
   constexpr T operator-() const noexcept{\
      T res = *this;\
      return res.into_opp();\
   }\
   constexpr T &operator+=(T const &rhs) noexcept{\
      lo += rhs.lo;\
      hi += rhs.hi + (lo<rhs.lo);\
      return *this;\
   }\
   constexpr T &operator-=(T const &rhs) noexcept{\
      U64 lo_bak = lo;\
      lo -= rhs.lo;\
      hi -= rhs.hi + (lo>lo_bak);\
      return *this;\
   }\
   constexpr T &operator&=(T const &rhs) noexcept{\
      lo &= rhs.lo;\
      hi &= rhs.hi;\
      return *this;\
   }\
   constexpr T &operator^=(T const &rhs) noexcept{\
      lo ^= rhs.lo;\
      hi ^= rhs.hi;\
      return *this;\
   }\
   constexpr T &operator|=(T const &rhs) noexcept{\
      lo |= rhs.lo;\
      hi |= rhs.hi;\
      return *this;\
   }\
   constexpr T &operator<<=(U8 shift) noexcept{\
      switch(shift>>6){\
      case 0:\
         hi = __shiftleft128(lo, hi, shift);\
         lo <<= shift;\
         break;\
      case 1:\
         hi = lo << (shift & 0x3F);\
         lo = 0;\
         break;\
      default:\
         lo = hi = 0;\
      }\
      return *this;\
   }
struct i128;
struct u128{
   DEF_TYPE(u128)
   constexpr u128 &operator>>=(U8 shift) noexcept{
      switch(shift>>6){
      case 0:
         lo = __shiftright128(lo, hi, shift);
         hi >>= shift;
         break;
      case 1:
         lo = hi >> (shift & 0x3F);
         hi = 0;
         break;
      default:
         lo = hi = 0;
      }
      return *this;
   }
   explicit constexpr operator i128() const noexcept;
};
struct i128{
   DEF_TYPE(i128)
   constexpr i128 &operator>>=(U8 shift) noexcept{
      switch(shift>>6){
      case 0:
         lo = __shiftright128(lo, hi, shift);
         hi >>= shift;
         break;
      case 1:
         lo = hi >> (shift & 0x3F);
         hi = -(hi>LLONG_MAX);
      default:
         lo = hi = -(hi>LLONG_MAX);
      }
      return *this;
   }
   explicit constexpr operator u128() const noexcept{
      u128 res;
      res.lo = lo; res.hi = hi;
      return res;
   }
};
#undef DEF_TYPE
constexpr u128::operator i128() const noexcept{
   i128 res;
   res.lo = lo; res.hi = hi;
   return res;
}
#define DEF_BIOP(T, OP)\
constexpr T operator OP(T lhs, T const &rhs) noexcept{\
   return lhs OP##= rhs;\
}
DEF_BIOP(u128, +)
DEF_BIOP(u128, -)
DEF_BIOP(u128, &)
DEF_BIOP(u128, ^)
DEF_BIOP(u128, |)
DEF_BIOP(i128, +)
DEF_BIOP(i128, -)
DEF_BIOP(i128, &)
DEF_BIOP(i128, ^)
DEF_BIOP(i128, |)
#undef DEF_BIOP
constexpr u128 operator<<(u128 self, U8 shift) noexcept{
   return self <<= shift;
}
constexpr u128 operator>>(u128 self, U8 shift) noexcept{
   return self >>= shift;
}
constexpr i128 operator<<(i128 self, U8 shift) noexcept{
   return self <<= shift;
}
constexpr i128 operator>>(i128 self, U8 shift) noexcept{
   return self >>= shift;
}
constexpr bool operator==(u128 const &lhs, u128 const &rhs) noexcept{
   return lhs.lo==rhs.lo && lhs.hi==rhs.hi;
}
constexpr bool operator!=(u128 const &lhs, u128 const &rhs) noexcept{
   return lhs.lo!=rhs.lo || lhs.hi!=rhs.hi;
}
constexpr bool operator==(i128 const &lhs, i128 const &rhs) noexcept{
   return lhs.lo==rhs.lo && lhs.hi==rhs.hi;
}
constexpr bool operator!=(i128 const &lhs, i128 const &rhs) noexcept{
   return lhs.lo!=rhs.lo || lhs.hi!=rhs.hi;
}
#define DEF_U128_CMP(OP)\
constexpr bool operator OP(u128 const &lhs, u128 const &rhs) noexcept{\
   return lhs.hi==rhs.hi? lhs.lo OP rhs.lo: lhs.hi OP rhs.hi;\
}
DEF_U128_CMP(<)
DEF_U128_CMP(<=)
DEF_U128_CMP(>)
DEF_U128_CMP(>=)
#undef DEF_U128_CMP
#define DEF_I128_CMP(OP)\
constexpr bool operator OP(i128 const &lhs, i128 const &rhs) noexcept{\
   constexpr U64 sign_bit = 0x8000'0000'0000'0000;\
   if(lhs.hi & sign_bit){\
      if((rhs.hi & sign_bit) == 0){\
         return -1 OP 0;\
      }\
      return lhs.hi==rhs.hi? lhs.lo OP rhs.lo: lhs.hi OP rhs.hi;\
   }\
   if(rhs.hi & sign_bit){\
      return -1 OP 0;\
   }\
   return lhs.hi==rhs.hi? lhs.lo OP rhs.lo: lhs.hi OP rhs.hi;\
}
DEF_I128_CMP(<)
DEF_I128_CMP(<=)
DEF_I128_CMP(>)
DEF_I128_CMP(>=)
#undef DEF_I128_CMP
#else // !defined(_MSC_BUILD)
using u128 = __uint128_t;
using i128 = __int128_t;
#endif

inline u128 &u128_mul_eq_u64(u128 &self, U64 m){
#ifdef _MSC_BUILD
   U64 carry;
   self.lo = _umul128(self.lo, m, &carry);
   self.hi = self.hi*m + carry;
#else
   self *= m;
#endif
   return self;
}

inline u128 &u128_div_eq_u64(u128 &self, U64 m, U64 *r=nullptr){
#ifdef _MSC_BUILD
   if(r != nullptr){
      *r = self.hi % m;
      self.hi /= m;
      self.lo = _udiv128(*r, self.lo, m, r);
   }else{
      U64 rem = self.hi % m;
      self.hi /= m;
      self.lo = _udiv128(rem, self.lo, m, &rem);
   }
#else
   if(r != nullptr){
      *r = self % m;
   }
   self /= m;
#endif
   return self;
}

struct u64div_t{
   U64 quot, rem;
};

inline u64div_t u64_mul_div(U64 a, U64 b, U64 m){
   u64div_t res;
#ifdef _MSC_BUILD
   U64 h, l = _umul128(a, b, &h);
   res.quot = _udiv128(h, l, m, &res.rem);
#else
   u128 p = static_cast<u128>(a)*b;
   res.quot = p/m, res.rem = p%m;
#endif
   return res;
}

struct i64div_t{
   I64 quot, rem;
};

inline i64div_t i64_mul_div(I64 a, I64 b, I64 m){
   i64div_t res;
#ifdef _MSC_BUILD
   I64 h, l = _mul128(a, b, &h);
   res.quot = _div128(h, l, m, &res.rem);
#else
   i128 p = static_cast<i128>(a)*b;
   res.quot = p/m, res.rem = p%m;
#endif
   return res;
}

#define DEF_STO128(S) inline S##128 sto##S##128(char const *s, size_t *pos=nullptr, int base=10){\
   if(base<0 || base==1 || base>36){\
      throw std::invalid_argument("expect base = 0, 2, 3, ..., 36");\
   }\
   S##128 res = 0;\
   char const *p = s;\
   int base_bak = base;\
   bool neg;\
   while(std::isspace(*p)){\
      ++p;\
   }\
   neg = false;\
   for(; ; ++p){\
      if(*p == '+'){\
         continue;\
      }\
      if(*p == '-'){\
         neg = !neg;\
         continue;\
      }\
      break;\
   }\
   if(base == 0){\
      if(!std::isdigit(*p)){\
         goto err_no_conversion;\
      }\
      if(*p == '0'){\
         if((p[1]=='X' || p[1]=='x') && std::isxdigit(p[2])){\
            p += 2;\
            base = 16;\
         }else{\
            base = 8;\
         }\
      }else{\
         base = 10;\
      }\
   }{\
      auto digit = [base](char c){\
         if(base > 10){\
            if(std::isdigit(c)){\
               return c-'0';\
            }\
            if('A'<=c && c<'A'+(base-10)){\
               return c-'A'+10;\
            }\
            if('a'<=c && c<'a'+(base-10)){\
               return c-'a'+10;\
            }\
            return -1;\
         }\
         if('0'<=c && c<'0'+base){\
            return c-'0';\
         }\
         return -1;\
      };\
      if(digit(*p) == -1){\
         goto err_no_conversion;\
      }\
      if(base_bak==16 && p[0]=='0' && (p[1]=='X' || p[1]=='x') && std::isxdigit(p[2])){\
         p += 2;\
      }\
      for(; ; ++p){\
         int d = digit(*p);\
         if(d == -1){\
            break;\
         }\
         /*As i128 and u128 share the same structure, the reinterpret_cast is safe.*/\
         u128_mul_eq_u64(*reinterpret_cast<u128*>(&res), base);\
         res += (neg? -d: d);\
      }\
   }\
   if(pos != nullptr){\
      *pos = p-s;\
   }\
   return res;\
err_no_conversion:\
   throw std::invalid_argument("no conversion");\
}
DEF_STO128(u)
DEF_STO128(i)
#undef DEF_STO128

namespace std{
inline string to_string(u128 value){
   string res;
   while(value){
      U64 d;
      u128_div_eq_u64(value, 10, &d);
      res += static_cast<char>('0'+d);
   }
   if(res.empty()){
      res = "0";
   }else{
      reverse(res.begin(), res.end());
   }
   return res;
}
inline string to_string(i128 value){
   string res;
   if(value >= 0){
      while(value){
         U64 d;
         u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
         res += static_cast<char>('0'+d);
      }
      if(res.empty()){
         res = "0";
      }
   }else{
      // As -i128::MIN = i128::MIN = 2^127, the negation is safe.
      value = static_cast<i128>(-*reinterpret_cast<u128*>(&value));
      while(value){
         U64 d;
         u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
         res += static_cast<char>('0'+d);
      }
      res += '-';
   }
   reverse(res.begin(), res.end());
   return res;
}
}

#undef U64
#undef I64
#undef U8

#endif // INT128_HH
