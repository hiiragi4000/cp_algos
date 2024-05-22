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
   constexpr DivT<T> div(T const &rhs) const noexcept;\
   constexpr T &operator/=(T const &rhs) noexcept;\
   constexpr T &operator%=(T const &rhs) noexcept;\
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
   constexpr bool neg() const noexcept{
      return hi & 0x8000'0000'0000'0000;
   }
   constexpr i128 &operator>>=(U8 shift) noexcept{
      if(shift == 0){
         ; // null statement
      }else if(shift>>6 == 0){
         lo = ((hi&((1ull<<shift)-1))<<(64-shift)) | (lo>>shift);
         hi = static_cast<U64>(static_cast<I64>(hi)>>shift);
      }else if(shift>>6 == 1){
         lo = static_cast<U64>(static_cast<I64>(hi)>>(shift&63));
         hi = neg()? -1: 0;
      }else{
         lo = hi = neg()? -1: 0;
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
   if(lhs.neg()){\
      if(!rhs.neg()){\
         return -1 OP 0;\
      }\
      return lhs.hi==rhs.hi? lhs.lo OP rhs.lo: lhs.hi OP rhs.hi;\
   }\
   if(rhs.neg()){\
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

constexpr u128 mul_u64_u64(U64 a, U64 b) noexcept{
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
   u128 res = mul_u64_u64(lo, rhs.lo);
   res.hi += lo*rhs.hi + hi*rhs.lo;
   return *this = res;
}
constexpr i128 &i128::operator*=(i128 const &rhs) noexcept{
   u128 res = mul_u64_u64(lo, rhs.lo);
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

#ifdef _MSC_BUILD
namespace impl{
#define SET_LOW_32BITS(SRC64, DEST32) (SRC64 = (SRC64 & static_cast<U64>(UINT_MAX)<<32) | static_cast<U32>(DEST32))
#define SET_HIGH_32BITS(SRC64, DEST32) (SRC64 = (SRC64 & UINT_MAX) | (static_cast<U64>(DEST32) << 32))
constexpr U32 diveq_u128_u32(u128 &self, U32 rhs) noexcept{
   U64 rem = self.hi>>32;
   SET_HIGH_32BITS(self.hi, rem/rhs);
   rem = (rem%rhs)<<32 | (self.hi&UINT_MAX);
   SET_LOW_32BITS(self.hi, rem/rhs);
   rem = (rem%rhs)<<32 | (self.lo>>32);
   SET_HIGH_32BITS(self.lo, rem/rhs);
   rem = (rem%rhs)<<32 | (self.lo&UINT_MAX);
   SET_LOW_32BITS(self.lo, rem/rhs);
   rem %= rhs;
   return static_cast<U32>(rem);
}
constexpr U32 muleq_u128_u32(u128 &self, U32 rhs) noexcept{
   U64 temp = (self.lo&UINT_MAX)*rhs;
   SET_LOW_32BITS(self.lo, temp);
   temp = (self.lo>>32)*rhs+(temp>>32);
   SET_HIGH_32BITS(self.lo, temp);
   temp = (self.hi&UINT_MAX)*rhs+(temp>>32);
   SET_LOW_32BITS(self.hi, temp);
   temp = (self.hi>>32)*rhs+(temp>>32);
   SET_HIGH_32BITS(self.hi, temp);
   return static_cast<U32>(temp>>32);
}
constexpr u128 mul_u64_u32(U64 a, U32 b) noexcept{
   u128 res{0, 0};
   U64 temp = (a&UINT_MAX)*b;
   SET_LOW_32BITS(res.lo, temp);
   temp = (a>>32)*b+(temp>>32);
   SET_HIGH_32BITS(res.lo, temp);
   res.hi = temp>>32;
   return res;
}
constexpr U64 diveq_u128_u64(u128 &self, U64 rhs) noexcept{
   if(rhs>>32 == 0){
      return diveq_u128_u32(self, static_cast<U32>(rhs));
   }
   U32 norm = static_cast<U32>((1ull<<32)/((rhs>>32)+1));
   rhs *= norm;
   U32 carry = muleq_u128_u32(self, norm);
   u128 rem{self.hi, carry};
   auto dwim = [&rem, rhs](){
      U32 q = static_cast<U32>((rem.hi<<32|(rem.lo>>32))/((rhs>>32)+1));
      rem -= mul_u64_u32(rhs, q);
      if(rem >= rhs){
         rem -= rhs;
         if(rem >= rhs){
            rem -= rhs;
            q += 2;
         }else{
            ++q;
         }
      }
      return q;
   };
   self.hi = dwim();
   rem = (rem<<32)|(self.lo>>32);
   SET_HIGH_32BITS(self.lo, dwim());
   rem = (rem<<32)|(self.lo&UINT_MAX);
   SET_LOW_32BITS(self.lo, dwim());
   return rem.lo / norm;
}
constexpr u128 mul_u128_u32(u128 lhs, U32 rhs) noexcept{
   muleq_u128_u32(lhs, rhs);
   return lhs;
}
constexpr u128 diveq_u128_u128(u128 &self, u128 const &rhs) noexcept{
   if(rhs.hi == 0){
      return diveq_u128_u64(self, rhs.lo);
   }
   if(rhs.hi >> 32){
      if(self.hi <= rhs.hi){
         u128 res = self;
         if(res >= rhs){
            res -= rhs;
            self = 1;
         }else{
            self = 0;
         }
         return res;
      }
      U32 q = static_cast<U32>(self.hi/(rhs.hi+1));
      u128 res = self-mul_u128_u32(rhs, q);
      if(res >= rhs){
         res -= rhs;
         ++q;
      }
      self = q;
      return res;
   }
   U32 norm = static_cast<U32>((1ull<<32)/((rhs.hi&UINT_MAX)+1));
   u128 nr = mul_u128_u32(rhs, norm);
   U32 carry = muleq_u128_u32(self, norm);
   u128 rem{(self.hi&UINT_MAX)<<32|(self.lo>>32), static_cast<U64>(carry)<<32|(self.hi>>32)};
   auto dwim = [&rem, &nr](){
      U32 q = static_cast<U32>(rem.hi/(nr.hi+1));
      rem -= mul_u128_u32(nr, q);
      if(rem >= nr){
         rem -= nr;
         if(rem >= nr){
            rem -= nr;
            q += 2;
         }else{
            ++q;
         }
      }
      return q;
   };
   self.hi = 0;
   SET_HIGH_32BITS(self.lo, dwim());
   rem = (rem<<32)|(self.lo&UINT_MAX);
   SET_LOW_32BITS(self.lo, dwim());
   diveq_u128_u32(rem, norm);
   return rem;
}
#undef SET_HIGH_32BITS
#undef SET_LOW_32BITS
} // namespace impl

constexpr DivT<u128> u128::div(u128 const &rhs) const noexcept{
   u128 q = *this, r = impl::diveq_u128_u128(q, rhs);
   return {q, r};
}
constexpr u128 &u128::operator/=(u128 const &rhs) noexcept{
   impl::diveq_u128_u128(*this, rhs);
   return *this;
}
constexpr u128 &u128::operator%=(u128 const &rhs) noexcept{
   return *this = impl::diveq_u128_u128(*this, rhs);
}
constexpr u128 operator/(u128 lhs, u128 const &rhs) noexcept{
   return lhs /= rhs;
}
constexpr u128 operator%(u128 lhs, u128 const &rhs) noexcept{
   return lhs %= rhs;
}

constexpr DivT<i128> i128::div(i128 const &rhs) const noexcept{
   u128 a = static_cast<u128>(*this), b = static_cast<u128>(rhs);
   if(neg()){
      a.into_opp();
   }
   if(rhs.neg()){
      b.into_opp();
   }
   auto [q, r] = a.div(b);
   if(neg()){
      if(!rhs.neg()){
         q.into_opp();
      }
      r.into_opp();
   }else if(rhs.neg()){
      q.into_opp();
   }
   return {static_cast<i128>(q), static_cast<i128>(r)};
}
constexpr i128 &i128::operator/=(i128 const &rhs) noexcept{
   return *this = div(rhs).quot;
}
constexpr i128 &i128::operator%=(i128 const &rhs) noexcept{
   return *this = div(rhs).rem;
}
constexpr i128 operator/(i128 lhs, i128 const &rhs) noexcept{
   return lhs /= rhs;
}
constexpr i128 operator%(i128 lhs, i128 const &rhs) noexcept{
   return lhs %= rhs;
}
#endif // _MSC_BUILD

constexpr DivT<U64> c_u64_mul_div(U64 a, U64 b, U64 m) noexcept{
   u128 c = mul_u64_u64(a, b);
#ifdef _MSC_BUILD
   U64 r = impl::diveq_u128_u64(c, m);
#else
   U64 r = c%m; c /= m;
#endif
   return {static_cast<U64>(c), r};
}

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
inline u128 &u128_div_eq_u64(u128 &self, U64 m, U64 *r){
#ifdef _MSC_BUILD
   *r = self.hi % m;
   self.hi /= m;
   self.lo = _udiv128(*r, self.lo, m, r);
#else
   *r = self % m;
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
