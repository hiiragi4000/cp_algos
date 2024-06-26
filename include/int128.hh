#ifndef INT128_HH
#define INT128_HH

#include"assert_platform.hh"
#include"util.hh"
#include<algorithm>
#include<stdexcept>
#include<string>
#include<type_traits>
#include<cctype>
#include<climits>
#include<cwctype>

#ifdef _MSC_BUILD
#include<intrin.h>
#endif

#define U8 unsigned char
#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

#ifdef _MSC_BUILD
#define DEF_TYPE(T)\
   U64 lo, hi;\
   T() = default;\
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>\
   constexpr T(INT a) noexcept: lo(a), hi(-(a<0)){}\
   constexpr T(U64 a, U64 b) noexcept: lo(a), hi(b){}\
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
   constexpr T &operator*=(T const &rhs) noexcept;\
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
      if(shift == 0){\
         ; /*null statement*/\
      }else if(shift>>6 == 0){\
         hi = hi<<shift | lo>>(64-shift);\
         lo <<= shift;\
      }else if(shift>>6 == 1){\
         hi = lo << (shift&63);\
         lo = 0;\
      }else{\
         hi = lo = 0;\
      }\
      return *this;\
   }
struct i128;
struct u128{
   DEF_TYPE(u128)
   constexpr u128 &operator>>=(U8 shift) noexcept{
      if(shift == 0){
         ; // null statement
      }else if(shift>>6 == 0){
         lo = ((hi&((1ull<<shift)-1))<<(64-shift)) | (lo>>shift);
         hi >>= shift;
      }else if(shift>>6 == 1){
         lo = hi >> (shift&63);
         hi = 0;
      }else{
         lo = hi = 0;
      }
      return *this;
   }
   explicit constexpr operator i128() const noexcept;
};
struct i128{
   DEF_TYPE(i128)
   constexpr i128 &operator>>=(U8 shift) noexcept{
      constexpr U64 sign_bit = 0x8000'0000'0000'0000;
      if(shift == 0){
         ; // null statement
      }else if(shift>>6 == 0){
         lo = ((hi&((1ull<<shift)-1))<<(64-shift)) | (lo>>shift);
         hi = static_cast<U64>(static_cast<I64>(hi)>>shift);
      }else if(shift>>6 == 1){
         lo = static_cast<U64>(static_cast<I64>(hi)>>(shift&63));
         hi = (hi&sign_bit)? -1: 0;
      }else{
         lo = hi = (hi&sign_bit)? -1: 0;
      }
      return *this;
   }
   explicit constexpr operator u128() const noexcept{
      return u128(lo, hi);
   }
};
constexpr u128::operator i128() const noexcept{
   return i128(lo, hi);
}
#undef DEF_TYPE

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
      return 0 OP -1;\
   }\
   return lhs.hi==rhs.hi? lhs.lo OP rhs.lo: lhs.hi OP rhs.hi;\
}
DEF_I128_CMP(<)
DEF_I128_CMP(<=)
DEF_I128_CMP(>)
DEF_I128_CMP(>=)
#undef DEF_I128_CMP
static_assert(std::is_standard_layout_v<u128>, "u128 is not of standard layout QQ");
static_assert(std::is_standard_layout_v<i128>, "i128 is not of standard layout QQ");

#else // !defined(_MSC_BUILD)
using u128 = __uint128_t;
using i128 = __int128_t;
#endif

constexpr u128 c_mul128(U64 a, U64 b) noexcept{
#ifdef _MSC_BUILD
   U32 al = a&UINT_MAX, ah = a>>32;
   U32 bl = b&UINT_MAX, bh = b>>32;
   u128 res{static_cast<U64>(al)*bh, 0};
   res += static_cast<U64>(ah)*bl;
   res <<= 32;
   res += static_cast<U64>(al)*bl;
   res.hi += static_cast<U64>(ah)*bh;
   return res;
#else
   return static_cast<u128>(a)*b;
#endif
}
#ifdef _MSC_BUILD
constexpr u128 &u128::operator*=(u128 const &rhs) noexcept{
   u128 res = c_mul128(lo, rhs.lo);
   res.hi += lo*rhs.hi + hi*rhs.lo;
   return *this = res;
}
constexpr i128 &i128::operator*=(i128 const &rhs) noexcept{
   u128 res = c_mul128(lo, rhs.lo);
   res.hi += lo*rhs.hi + hi*rhs.lo;
   return *this = static_cast<i128>(res);
}
constexpr u128 operator*(u128 lhs, u128 const &rhs) noexcept{
   return lhs *= rhs;
}
constexpr i128 operator*(i128 lhs, i128 const &rhs) noexcept{
   return lhs *= rhs;
}
#endif

namespace impl{
constexpr int u64_mul_10_carry(U64 a) noexcept{
   constexpr U64 u = ULLONG_MAX/5;
   int res = a>>63? 5: 0;
   a &= LLONG_MAX;
   if(a <= u){
      res += a > u/2;
   }else if(a <= u+u/2){
      res += 2;
   }else if(a <= 2*u){
      res += 3;
   }else{
      res += 4;
   }
   return res;
}
constexpr u128 &u128_mul_eq_base(u128 &self, int base) noexcept{
#ifdef _MSC_BUILD
   if(base == 2){
      self <<= 1;
   }else if(base == 8){
      self <<= 3;
   }else if(base == 16){
      self <<= 4;
   }else{
      int carry = u64_mul_10_carry(self.lo);
      self.lo *= 10;
      self.hi = 10*self.hi + carry;
   }
#else
   self *= base;
#endif
   return self;
}
constexpr i128 &i128_mul_eq_base(i128 &self, int base) noexcept{
#ifdef _MSC_BUILD
   u128 uself = static_cast<u128>(self);
   self = static_cast<i128>(u128_mul_eq_base(uself, base));
#else
   self *= base;
#endif
   return self;
}
} // namespace impl

#define DEF_LITERAL(T)\
T operator""_##T(long double) = delete;\
constexpr T operator""_##T(char const *s) noexcept{\
   int base = 10;\
   if(*s == '0'){\
      if('0'<=s[1] && s[1]<='7'){\
         base = 8;\
         ++s;\
      }else if(s[1]=='B' || s[1]=='b'){\
         base = 2;\
         s += 2;\
      }else{\
         base = 16;\
         s += 2;\
      }\
   }\
   T res = 0;\
   auto char_to_digit = [](char c){\
      if(c <= '9'){\
         return c-'0';\
      }\
      if(c <= 'F'){\
         return c-'A'+10;\
      }\
      return c-'a'+10;\
   };\
   for(; *s; ++s){\
      if(*s == '\''){\
         continue;\
      }\
      impl::T##_mul_eq_base(res, base) += char_to_digit(*s);\
   }\
   return res;\
}
DEF_LITERAL(u128)
DEF_LITERAL(i128)
#undef DEF_LITERAL

inline DivT<U64> u64_mul_div(U64 a, U64 b, U64 m){
   DivT<U64> res;
#ifdef _MSC_BUILD
   U64 h, l = _umul128(a, b, &h);
   res.quot = _udiv128(h, l, m, &res.rem);
#else
   u128 p = static_cast<u128>(a)*b;
   res.quot = p/m, res.rem = p%m;
#endif
   return res;
}

inline DivT<I64> i64_mul_div(I64 a, I64 b, I64 m){
   DivT<I64> res;
#ifdef _MSC_BUILD
   I64 h, l = _mul128(a, b, &h);
   res.quot = _div128(h, l, m, &res.rem);
#else
   i128 p = static_cast<i128>(a)*b;
   res.quot = p/m, res.rem = p%m;
#endif
   return res;
}

namespace impl{
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
} // namespace impl

#define DEF_STOI(T) inline T sto##T(char const *s, size_t *pos=nullptr, int base=10){\
   if(base<0 || base==1 || base>36){\
      throw std::invalid_argument("expect base = 0, 2, 3, ..., 36");\
   }\
   T res = 0;\
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
         impl::u128_mul_eq_u64(*reinterpret_cast<u128*>(&res), base);\
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
DEF_STOI(u128)
DEF_STOI(i128)
#undef DEF_STOI

#define DEF_WSTOI(T) inline T sto##T(wchar_t const *s, size_t *pos=nullptr, int base=10){\
   if(base<0 || base==1 || base>36){\
      throw std::invalid_argument("expect base = 0, 2, 3, ..., 36");\
   }\
   T res = 0;\
   wchar_t const *p = s;\
   int base_bak = base;\
   bool neg;\
   while(std::iswspace(*p)){\
      ++p;\
   }\
   neg = false;\
   for(; ; ++p){\
      if(*p == L'+'){\
         continue;\
      }\
      if(*p == L'-'){\
         neg = !neg;\
         continue;\
      }\
      break;\
   }\
   if(base == 0){\
      if(!std::iswdigit(*p)){\
         goto err_no_conversion;\
      }\
      if(*p == L'0'){\
         if((p[1]==L'X' || p[1]==L'x') && std::iswxdigit(p[2])){\
            p += 2;\
            base = 16;\
         }else{\
            base = 8;\
         }\
      }else{\
         base = 10;\
      }\
   }{\
      auto digit = [base](wchar_t c){\
         if(base > 10){\
            if(std::iswdigit(c)){\
               return c-L'0';\
            }\
            if(L'A'<=c && c<L'A'+(base-10)){\
               return c-L'A'+10;\
            }\
            if(L'a'<=c && c<L'a'+(base-10)){\
               return c-L'a'+10;\
            }\
            return -1;\
         }\
         if(L'0'<=c && c<L'0'+base){\
            return c-L'0';\
         }\
         return -1;\
      };\
      if(digit(*p) == -1){\
         goto err_no_conversion;\
      }\
      if(base_bak==16 && p[0]==L'0' && (p[1]==L'X' || p[1]==L'x') && std::iswxdigit(p[2])){\
         p += 2;\
      }\
      for(; ; ++p){\
         int d = digit(*p);\
         if(d == -1){\
            break;\
         }\
         impl::u128_mul_eq_u64(*reinterpret_cast<u128*>(&res), base);\
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
DEF_WSTOI(u128)
DEF_WSTOI(i128)
#undef DEF_WSTOI

namespace std{
inline string to_string(u128 value){
   string res;
   while(value){
      U64 d;
      impl::u128_div_eq_u64(value, 10, &d);
      res += static_cast<char>('0'+d);
   }
   if(res.empty()){
      res = "0";
   }else{
      reverse(res.begin(), res.end());
   }
   return res;
}
inline wstring to_wstring(u128 value){
   wstring res;
   while(value){
      U64 d;
      impl::u128_div_eq_u64(value, 10, &d);
      res += static_cast<wchar_t>(L'0'+d);
   }
   if(res.empty()){
      res = L"0";
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
         impl::u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
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
         impl::u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
         res += static_cast<char>('0'+d);
      }
      res += '-';
   }
   reverse(res.begin(), res.end());
   return res;
}
inline wstring to_wstring(i128 value){
   wstring res;
   if(value >= 0){
      while(value){
         U64 d;
         impl::u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
         res += static_cast<wchar_t>(L'0'+d);
      }
      if(res.empty()){
         res = L"0";
      }
   }else{
      value = static_cast<i128>(-*reinterpret_cast<u128*>(&value));
      while(value){
         U64 d;
         impl::u128_div_eq_u64(*reinterpret_cast<u128*>(&value), 10, &d);
         res += static_cast<wchar_t>(L'0'+d);
      }
      res += L'-';
   }
   reverse(res.begin(), res.end());
   return res;
}
}

#undef U64
#undef I64
#undef U32
#undef U8

#endif // INT128_HH
