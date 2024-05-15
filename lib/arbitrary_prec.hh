#ifndef ARBITRARY_PREC_HH
#define ARBITRARY_PREC_HH

#include"convolution.hh"
#include<algorithm>
#include<array>
#include<stdexcept>
#include<string>
#include<type_traits>
#include<utility>
#include<vector>
#include<cctype>
#include<csignal>
#include<cstdlib>
#include<cwctype>

#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

#define LOWBIT(N) ((N)&(~(N)+1))

namespace impl{
// Convert a[0: n-1] from base |b1| to base |b2| in-place.
// Require
// 1. n = 2^k for some k >= 1.
// 2. 0 <= a[i] < |b1|.
// 3. 1 < |b1| < |b2|.
inline void bigint_change_base(int *a, size_t n, U32 b1, U32 b2){
   Fft fft;
   std::vector<int> pb(n/2);
   std::vector<std::complex<double>> fa(n), fpb(n);
   std::vector<U32> na(n), npb(n);
   pb[0] = b1;
   for(size_t l=1; ; l*=2){
      for(size_t i=0; i<l; ++i){
         fpb[i] = npb[i] = pb[i];
      }
      fft.transform_in_place(false, fpb.begin(), 2*l);
      NttU32::transform_in_place(false, npb.begin(), 2*l);
      for(size_t i=0; i<n; i+=2*l){
         for(size_t j=0; j<l; ++j){
            fa[j] = na[j] = a[i+l+j];
         }
         std::fill(fa.begin()+l, fa.begin()+2*l, 0);
         std::fill(na.begin()+l, na.begin()+2*l, 0);
         fft.transform_in_place(false, fa.begin(), 2*l);
         NttU32::transform_in_place(false, na.begin(), 2*l);
         for(size_t j=0; j<2*l; ++j){
            fa[j] *= fpb[j];
            na[j] = static_cast<U64>(na[j])*npb[j]%NttU32::prime();
         }
         fft.transform_in_place(true, fa.begin(), 2*l);
         NttU32::transform_in_place(true, na.begin(), 2*l);
         std::fill(a+i+l, a+i+2*l, 0);
         U64 carry = 0;
         for(size_t j=0; j<2*l; ++j){
            U64 q = static_cast<U64>(std::round((fa[j].real()-na[j])/NttU32::prime()));
            U32 q1 = static_cast<U32>(q/b2), r1 = q%b2;
            U32 q2 = NttU32::prime()/b2, r2 = NttU32::prime()%b2;
            U32 q3 = na[j]/b2, r3 = na[j]%b2;
            U32 q4 = static_cast<U32>(carry/b2), r4 = carry%b2;
            carry = static_cast<U64>(q1)*q2*b2 + static_cast<U64>(q1)*r2 + q2*r1 + q3 + q4;
            U64 temp = static_cast<U64>(r1)*r2 + r3 + r4 + a[i+j];
            carry += temp/b2;
            a[i+j] = temp%b2;
         }
      }
      if(l == n/2){
         break;
      }
      for(size_t i=0; i<2*l; ++i){
         fpb[i] *= fpb[i];
         npb[i] = static_cast<U64>(npb[i])*npb[i]%NttU32::prime();
      }
      fft.transform_in_place(true, fpb.begin(), 2*l);
      NttU32::transform_in_place(true, npb.begin(), 2*l);
      U64 carry = 0;
      for(size_t i=0; i<2*l; ++i){
         U64 q = static_cast<U64>(std::round((fpb[i].real()-npb[i])/NttU32::prime()));
         U32 q1 = static_cast<U32>(q/b2), r1 = q%b2;
         U32 q2 = NttU32::prime()/b2, r2 = NttU32::prime()%b2;
         U32 q3 = npb[i]/b2, r3 = npb[i]%b2;
         U32 q4 = static_cast<U32>(carry/b2), r4 = carry%b2;
         carry = static_cast<U64>(q1)*q2*b2 + static_cast<U64>(q1)*r2 + q2*r1 + q3 + q4;
         U64 temp = static_cast<U64>(r1)*r2 + r3 + r4;
         carry += temp/b2;
         pb[i] = temp%b2;
      }
   }
}
} // namespace impl

template<typename T> struct DivT{
   T quot, rem;
};

struct BigInt: private std::vector<int>{
   constexpr static int BASE = 100'000'000;
   constexpr static int LG10_BASE = 8;
   constexpr static int SHRT_BASE = 1'000'000;
   constexpr static int LG10_SHRT_BASE = 6;
   constexpr static double LGBASE_2 = .03762874945799764940;
   constexpr static std::array<int, LG10_BASE+1> EXP10 = [](){
      std::array<int, LG10_BASE+1> res{1};
      for(int i=1; i<=LG10_BASE; ++i){
         res[i] = 10*res[i-1];
      }
      return res;
   }();
   constexpr static std::array<int, 37> LGD_BASE = [](){
      std::array<int, 37> res{};
      for(int i=2; i<=36; ++i){
         int b = BASE, ri = 0;
         while(b >= i){
            b /= i; ++ri;
         }
         res[i] = ri;
      }
      return res;
   }();
   constexpr static std::array<int, 37> BASE_CHG = [](){
      std::array<int, 37> res{};
      for(int i=2; i<=36; ++i){
         res[i] = 1;
         for(int t=LGD_BASE[i]; t-->0; ){
            res[i] *= i;
         }
      }
      return res;
   }();
   using Base = std::vector<int>;
   using Base::capacity;
   using Base::get_allocator;
   using Base::reserve;
   using Base::shrink_to_fit;
   using Base::size;
   BigInt() = default;
   BigInt(U64 a){
      while(a > 0){
         push_back(a%BASE);
         a /= BASE;
      }
   }
   BigInt(I64 a){
      int neg = a<0;
      while(a != 0){
         push_back(a%BASE);
         a /= BASE;
      }
      if(neg){
         for(size_t i=0; i<size()-1; ++i){
            data()[i] = -data()[i];
         }
      }
   }
   BigInt(unsigned long a): BigInt(static_cast<U64>(a)){}
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>
   BigInt(INT a): BigInt(static_cast<I64>(a)){}
   std::vector<int> into_vec() noexcept{
      std::vector<int> res = std::move(*this);
      return res;
   }
   explicit operator bool() const noexcept{
      return !empty();
   }
   explicit operator std::string() const{
      if(empty()){
         return "0";
      }
      std::string res;
      int t;
      for(size_t i=0; i<size()-1; ++i){
         t = data()[i];
         for(int j=LG10_BASE; j-->0; ){
            res += static_cast<char>('0'+t%10);
            t /= 10;
         }
      }
      t = std::abs(back());
      while(t > 0){
         res += static_cast<char>('0'+t%10);
         t /= 10;
      }
      if(back() < 0){
         res += '-';
      }
      std::reverse(res.begin(), res.end());
      return res;
   }
   std::string to_string() const{
      return static_cast<std::string>(*this);
   }
   std::string to_string(int base) const{
      if(base == 10){
         return static_cast<std::string>(*this);
      }
      if(base<2 || base>36){
         throw std::out_of_range("base out of range [2, 36]");
      }
      std::vector<int> v = to_vec_in_base(base);
      if(v.empty()){
         return "0";
      }
      std::string res;
      auto digit = [](int d){
         return static_cast<char>(d<10? '0'+d: 'A'+(d-10));
      };
      for(size_t i=0; i<v.size()-1; ++i){
         for(int j=LGD_BASE[base]; j-->0; ){
            res += digit(v[i]%base);
            v[i] /= base;
         }
      }
      bool neg = false;
      if(v.back() < 0){
         neg = true;
         v.back() = -v.back();
      }
      while(v.back() > 0){
         res += digit(v.back()%base);
         v.back() /= base;
      }
      if(neg){
         res += '-';
      }
      std::reverse(res.begin(), res.end());
      return res;
   }
   explicit operator std::wstring() const{
      if(empty()){
         return L"0";
      }
      std::wstring res;
      int t;
      for(size_t i=0; i<size()-1; ++i){
         t = data()[i];
         for(int j=LG10_BASE; j-->0; ){
            res += static_cast<wchar_t>(L'0'+t%10);
            t /= 10;
         }
      }
      t = std::abs(back());
      while(t > 0){
         res += static_cast<wchar_t>(L'0'+t%10);
         t /= 10;
      }
      if(back() < 0){
         res += L'-';
      }
      std::reverse(res.begin(), res.end());
      return res;
   }
   std::wstring to_wstring() const{
      return static_cast<std::wstring>(*this);
   }
   std::wstring to_wstring(int base) const{
      if(base == 10){
         return static_cast<std::wstring>(*this);
      }
      if(base<2 || base>36){
         throw std::out_of_range("base out of range [2, 36]");
      }
      std::vector<int> v = to_vec_in_base(base);
      if(v.empty()){
         return L"0";
      }
      std::wstring res;
      auto digit = [](int d){
         return static_cast<wchar_t>(d<10? L'0'+d: L'A'+(d-10));
      };
      for(size_t i=0; i<v.size()-1; ++i){
         for(int j=LGD_BASE[base]; j-->0; ){
            res += digit(v[i]%base);
            v[i] /= base;
         }
      }
      bool neg = false;
      if(v.back() < 0){
         neg = true;
         v.back() = -v.back();
      }
      while(v.back() > 0){
         res += digit(v.back()%base);
         v.back() /= base;
      }
      if(neg){
         res += L'-';
      }
      std::reverse(res.begin(), res.end());
      return res;
   }
   int sgn() const noexcept{
      return empty()? 0: back()>0? 1: -1;
   }
   BigInt operator+() const{
      return *this;
   }
   BigInt &into_opp() noexcept{
      if(!empty()){
         back() = -back();
      }
      return *this;
   }
   BigInt operator-() const{
      BigInt res = *this;
      return res.into_opp();
   }
   BigInt &into_abs() noexcept{
      if(!empty()){
         back() = std::abs(back());
      }
      return *this;
   }
   BigInt &operator++(){
      if(sgn() == -1){
         bool borrow = true;
         for(size_t i=0; i<size()-1; ++i){
            --data()[i];
            if(data()[i] == -1){
               data()[i] = BASE-1;
            }else{
               borrow = false;
               break;
            }
         }
         if(borrow){
            ++back();
            if(back() == 0){
               pop_back();
               into_opp();
            }
         }
      }else{
         bool carry = true;
         for(size_t i=0; i<size(); ++i){
            ++data()[i];
            if(data()[i] == BASE){
               data()[i] = 0;
            }else{
               carry = false;
               break;
            }
         }
         if(carry){
            push_back(1);
         }
      }
      return *this;
   }
   BigInt operator++(int){
      BigInt res = *this;
      ++*this;
      return res;
   }
   BigInt &operator--(){
      if(sgn() == -1){
         bool carry = true;
         for(size_t i=0; i<size()-1; ++i){
            ++data()[i];
            if(data()[i] == BASE){
               data()[i] = 0;
            }else{
               carry = false;
               break;
            }
         }
         if(carry){
            --back();
            if(back() == -BASE){
               back() = 0;
               push_back(-1);
            }
         }
      }else if(sgn() == 0){
         push_back(-1);
      }else{
         for(size_t i=0; i<size(); ++i){
            --data()[i];
            if(data()[i] == -1){
               data()[i] = BASE-1;
            }else{
               break;
            }
         }
         if(back() == 0){
            pop_back();
         }
      }
      return *this;
   }
   BigInt operator--(int){
      BigInt res = *this;
      --*this;
      return res;
   }
   BigInt &operator+=(BigInt const &rhs){
      if(sgn() == 0){
         return *this = rhs;
      }
      if(rhs.sgn() == 0){
         return *this;
      }
      if(sgn() == 1){
         if(rhs.sgn() == 1){
            return add_std(rhs, false);
         }
         if(cmp_abs(rhs) >= 0){
            return sub_opp(rhs, false);
         }
         return sub_from_opp(rhs, false);
      }
      if(rhs.sgn() == -1){
         return into_opp().add_std(rhs, true).into_opp();
      }
      if(cmp_abs(rhs) >= 0){
         return into_opp().sub_opp(rhs, true).into_opp();
      }
      return into_opp().sub_from_opp(rhs, true).into_opp();
   }
   BigInt &operator-=(BigInt const &rhs){
      return (into_opp()+=rhs).into_opp();
   }
#define THRES_MUL 128
   BigInt &operator*=(BigInt const &rhs){
      if(sgn() == 0){
         return *this;
      }
      if(rhs.sgn() == 0){
         clear();
         return *this;
      }
      int sign = sgn() * rhs.sgn();
      if(sgn() == -1){
         into_opp();
      }
      if(size()<=THRES_MUL || rhs.size()<=THRES_MUL){
         naive_mul(rhs);
      }else{
         fft_mul(rhs);
      }
      if(sign == -1){
         into_opp();
      }
      return *this;
   }
#undef THRES_MUL
#define THRES_DIV 64
   DivT<BigInt> div(BigInt const &rhs) const{
      if(rhs.sgn() == 0){
         std::raise(SIGFPE);
      }
      if(size() < rhs.size()){
         return {0, *this};
      }
      BigInt b = rhs;
      b.into_abs();
      DivT<BigInt> qr;
      if(b.size()<=THRES_DIV || size()-b.size()<=THRES_DIV){
         qr = naive_div(b);
      }else{
         qr = newton_div(b);
      }
      if(rhs.sgn() == -1){
         qr.quot.into_opp();
      }
      return qr;
   }
#undef THRES_DIV
   BigInt &operator/=(BigInt const &rhs){
      return *this = div(rhs).quot;
   }
   BigInt &operator%=(BigInt const &rhs){
      return *this = div(rhs).rem;
   }
   BigInt &into_bitwise_not(){
      return --into_opp();
   }
   BigInt operator~() const{
      BigInt res = *this;
      return res.into_bitwise_not();
   }
   BigInt pow(size_t rhs) const{
      if(rhs == 0){
         return 1;
      }
      if(rhs==1 || empty() || (size()==1 && back()==1)){
         return *this;
      }
      if(size()==1 && back()==-1){
         return rhs%2==0? 1: -1;
      }
      std::vector<size_t> p{rhs};
      while(p.back() > 1){
         p.push_back(p.back()/2);
      }
      std::reverse(p.begin(), p.end());
      BigInt res = *this;
      for(size_t i=1; i<p.size(); ++i){
         res *= res;
         if(p[i]%2 == 1){
            res *= *this;
         }
      }
      return res;
   }
   BigInt &operator&=(BigInt const &rhs){
      if(empty()){
         return *this;
      }
      if(rhs.empty()){
         clear();
         return *this;
      }
      if(size()==1 && back()==-1){
         return *this = rhs;
      }
      if(rhs.size()==1 && rhs.back()==-1){
         return *this;
      }
      if(*this == rhs){
         return *this;
      }
      if(size() < rhs.size()){
         BigInt res = rhs;
         return *this = res.nontrivial_and_eq(*this);
      }
      return nontrivial_and_eq(rhs);
   }
   BigInt &operator^=(BigInt const &rhs){
      if(empty()){
         return *this = rhs;
      }
      if(rhs.empty()){
         return *this;
      }
      if(size()==1 && back()==-1){
         return *this = ~rhs;
      }
      if(rhs.size()==1 && rhs.back()==-1){
         return into_bitwise_not();
      }
      if(*this == rhs){
         clear();
         return *this;
      }
      if(size() < rhs.size()){
         BigInt res = rhs;
         return *this = res.nontrivial_xor_eq(*this);
      }
      return nontrivial_xor_eq(rhs);
   }
   BigInt &operator|=(BigInt const &rhs){
      return (into_bitwise_not() &= ~rhs).into_bitwise_not();
   }
   BigInt &operator<<=(size_t rhs){
      if(empty()){
         return *this;
      }
      return *this *= BigInt(2).pow(rhs);
   }
   BigInt &operator>>=(size_t rhs){
      if(empty()){
         return *this;
      }
      if(size() < std::floor(LGBASE_2*rhs)+1){
         if(back() > 0){
            clear();
            return *this;
         }
         return *this = -1;
      }
      auto [q, r] = div(BigInt(2).pow(rhs));
      if(r.sgn() == -1){
         --q;
      }
      return *this = std::move(q);
   }
private:
   BigInt &trunc() noexcept{
      while(!empty() && back()==0){
         pop_back();
      }
      return *this;
   }
   int cmp_abs(BigInt const &rhs) const noexcept{
      if(size() != rhs.size()){
         return (size()>rhs.size())-(size()<rhs.size());
      }
      if(empty()){
         return 0;
      }
      int l = std::abs(back()), r = std::abs(rhs.back());
      if(l != r){
         return (l>r)-(l<r);
      }
      for(size_t i=size()-1; i-->0; ){
         if(data()[i] != rhs[i]){
            return (data()[i]>rhs[i])-(data()[i]<rhs[i]);
         }
      }
      return 0;
   }
   BigInt &add_std(BigInt const &rhs, bool this_opp){
      if(size() < rhs.size()){
         resize(rhs.size());
      }
      bool carry = false;
      for(size_t i=0; i<rhs.size(); ++i){
         data()[i] += (this_opp && i==rhs.size()-1? -rhs[i]: rhs[i]) + carry;
         if(data()[i] >= BASE){
            data()[i] -= BASE;
            carry = true;
         }else{
            carry = false;
         }
      }
      if(carry){
         for(size_t i=rhs.size(); i<size(); ++i){
            ++data()[i];
            if(data()[i] < BASE){
               carry = false;
               break;
            }
            data()[i] = 0;
         }
         if(carry){
            push_back(1);
         }
      }
      return *this;
   }
   BigInt &sub_opp(BigInt const &rhs, bool this_opp){
      for(size_t i=0; i<rhs.size()-1; ++i){
         data()[i] -= rhs[i];
         if(data()[i] < 0){
            data()[i] += BASE;
            --data()[i+1];
         }
      }
      data()[rhs.size()-1] += this_opp? -rhs.back(): rhs.back();
      for(size_t i=rhs.size()-1; data()[i]<0; ++i){
         data()[i] += BASE;
         --data()[i+1];
      }
      return trunc();
   }
   BigInt &sub_from_opp(BigInt const &rhs, bool this_opp){
      resize(rhs.size());
      bool borrow = false;
      for(size_t i=0; i<size()-1; ++i){
         data()[i] = rhs[i]-data()[i]-borrow;
         if(data()[i] < 0){
            data()[i] += BASE;
            borrow = true;
         }else{
            borrow = false;
         }
      }
      back() = (this_opp? rhs.back(): -rhs.back())-back()-borrow;
      return trunc().into_opp();
   }
   BigInt &naive_mul(BigInt const &rhs){
      BigInt res;
      res.resize(size()+rhs.size());
      U64 carry, new_digit;
      for(size_t i=0; i<size(); ++i){
         carry = 0;
         for(size_t j=0; j<rhs.size()-1; ++j){
            new_digit = res[i+j]+static_cast<U64>(data()[i])*rhs[j]+carry;
            res[i+j] = new_digit % BASE;
            carry = new_digit / BASE;
         }
         new_digit = res[i+rhs.size()-1]+static_cast<U64>(data()[i])*std::abs(rhs.back())+carry;
         res[i+rhs.size()-1] = new_digit % BASE;
         carry = new_digit / BASE;
         if(carry){
            res[i+rhs.size()] = static_cast<int>(carry);
         }
      }
      if(res.back() == 0){
         res.pop_back();
      }
      return *this = res;
   }
   BigInt &fft_mul(BigInt const &rhs){
      BigInt res;
      res.resize(size()+rhs.size());
      size_t n = size()+rhs.size()-1;
      while(n != LOWBIT(n)){
         n += LOWBIT(n);
      }
      std::vector<std::complex<double>> z(n);
      for(size_t i=0; i<size(); ++i){
         z[i] = data()[i];
      }
      for(size_t i=0; i<rhs.size()-1; ++i){
         z[i] += std::complex<double>(0, rhs[i]);
      }
      z[rhs.size()-1] += std::complex<double>(0, std::abs(rhs.back()));
      Fft fft;
      fft.transform_in_place(false, z.begin(), n);
      for(size_t i=0; i<n; ++i){
         z[i] *= z[i];
      }
      fft.transform_in_place(true, z.begin(), n);
      std::vector<U32> x(n), y(n);
      for(size_t i=0; i<size(); ++i){
         x[i] = data()[i];
      }
      for(size_t i=0; i<rhs.size()-1; ++i){
         y[i] = rhs[i];
      }
      y[rhs.size()-1] = std::abs(rhs.back());
      NttU32::transform_in_place(false, x.begin(), n);
      NttU32::transform_in_place(false, y.begin(), n);
      for(size_t i=0; i<n; ++i){
         x[i] = static_cast<U64>(x[i])*y[i]%NttU32::prime();
      }
      NttU32::transform_in_place(true, x.begin(), n);
      U64 carry = 0;
      for(size_t i=0; i<size()+rhs.size()-1; ++i){
         U64 q = static_cast<U64>(std::round((z[i].imag()/2-x[i])/NttU32::prime()));
         U32 q1 = static_cast<U32>(q/BASE), r1 = q%BASE;
         constexpr U32 q2 = NttU32::prime()/BASE, r2 = NttU32::prime()%BASE;
         U32 q3 = x[i]/BASE, r3 = x[i]%BASE;
         U32 q4 = static_cast<U32>(carry/BASE), r4 = carry%BASE;
         carry = static_cast<U64>(q1)*q2*BASE + static_cast<U64>(q1)*r2 + q2*r1 + q3 + q4;
         U64 temp = static_cast<U64>(r1)*r2 + r3 + r4;
         carry += temp/BASE;
         res[i] = temp%BASE;
      }
      if(carry){
         res.back() = static_cast<int>(carry);
      }else{
         res.pop_back();
      }
      return *this = res;
   }
   int base_div_eq(int b){
      if(empty()){
         return 0;
      }
      U64 r = 0;
      for(size_t i=size(); i-->0; ){
         r = BASE*r + data()[i];
         U64 qi = r/b;
         data()[i] = static_cast<int>(qi);
         r -= qi*b;
      }
      if(back() == 0){
         pop_back();
      }
      return static_cast<int>(r);
   }
   DivT<BigInt> naive_div(BigInt &b) const{
      BigInt a = *this; a.into_abs();
      if(b.size() == 1){
         int r = a.base_div_eq(b[0]);
         if(sgn() == -1){
            a.into_opp(); r = -r;
         }
         return {std::move(a), r};
      }
      int c = BASE/(b.back()+1);
      a *= c; b *= c;
      BigInt q, qib;
      q.resize(a.size()-b.size()+1);
      U64 h = 0;
      auto a_sub_at = [a=a.data()](BigInt const &rhs, size_t i){
         for(size_t j=0; j<rhs.size(); ++j){
            a[i+j] -= rhs[j];
            if(a[i+j] < 0){
               a[i+j] += BASE;
               --a[i+j+1];
            }
         }
      };
      auto a_cmp_at = [&a](BigInt const &rhs, size_t i){
         if(a.size()>i+rhs.size() && a[i+rhs.size()]>0){
            return 1;
         }
         for(size_t j=rhs.size(); j-->0; ){
            if(a[i+j] != rhs[j]){
               return (a[i+j]>rhs[j]) - (a[i+j]<rhs[j]);
            }
         }
         return 0;
      };
      for(size_t i=q.size(); i-->0; ){
         h = BASE*h + a[i+b.size()-1];
         int qi = static_cast<int>(h/(b.back()+1));
         qib = b; qib *= qi;
         a_sub_at(qib, i);
         if(a_cmp_at(b, i) >= 0){
            a_sub_at(b, i);
            if(a_cmp_at(b, i) >= 0){
               a_sub_at(b, i);
               qi += 2;
            }else{
               ++qi;
            }
         }
         q[i] = qi;
         h = a[i+b.size()-1];
      }
      if(q.back() == 0){
         q.pop_back();
      }
      a.trunc().base_div_eq(c);
      if(sgn() == -1){
         q.into_opp(); a.into_opp();
      }
      return {std::move(q), std::move(a)};
   }
   DivT<BigInt> newton_div(BigInt &b) const{
      BigInt a = *this; a.into_abs();
      int c = BASE/(b.back()+1);
      a *= c; b *= c;
      BigInt x = static_cast<U64>(BASE)*BASE/(b.back()+1);
      std::vector<size_t> d{a.size()-b.size()};
      while(d.back() > 1){
         d.push_back((d.back()+1)/2);
      }
      std::reverse(d.begin(), d.end());
      size_t prev = 0;
      for(size_t di: d){
         BigInt y;
         y.resize(di+1);
         std::copy_n(b.crbegin(), std::min(b.size(), y.size()), y.rbegin());
         ++y;
         BigInt z = x;
         z *= y;
         z.erase(z.cbegin(), z.cbegin()+prev+1);
         ++z;
         BigInt t;
         t.resize(di+2); t.back() = 2;
         t -= z;
         x *= t;
         x.erase(x.cbegin(), x.cbegin()+prev+1);
         prev = di;
      }
      BigInt q = a;
      q *= x;
      q.erase(q.cbegin(), q.cbegin()+a.size()+1);
      BigInt r = -q;
      r *= b; r += a;
      if(r.cmp_abs(b) >= 0){
         ++q; r -= b;
      }
      r.base_div_eq(c);
      if(sgn() == -1){
         q.into_opp(); r.into_opp();
      }
      return {std::move(q), std::move(r)};
   }
   // Convert |*this| to an std::vector<int> in base |BASE_CHG[base]|.
   // Require
   // 1. |base| = 2, 3, ..., 9, 11, 12, ..., 36.
   // 2. |*this| is well-formed.
   std::vector<int> to_vec_in_base(int base) const{
      if(empty() || (size()==1 && std::abs(data()[0])<BASE_CHG[base])){
         return *this;
      }
      bool neg = false;
      if(sgn() == -1){
         neg = true;
      }
      std::vector<int> a;
      I64 d = 0;
      int rem = 0;
      for(size_t i=0; i<size(); ++i){
         d += static_cast<I64>(EXP10[rem])*std::abs(data()[i]);
         rem += LG10_BASE;
         while(rem >= LG10_SHRT_BASE){
            a.push_back(d%SHRT_BASE);
            d /= SHRT_BASE;
            rem -= LG10_SHRT_BASE;
         }
      }
      if(rem > 0){
         a.push_back(static_cast<int>(d));
      }
      size_t n = a.size();
      while(n != LOWBIT(n)){
         n += LOWBIT(n);
      }
      a.resize(n);
      impl::bigint_change_base(a.data(), n, SHRT_BASE, BASE_CHG[base]);
      while(a.back() == 0){
         a.pop_back();
      }
      if(neg){
         a.back() = -a.back();
      }
      return a;
   }
   // Move the sequence a[] in base |BASE_CHG[base]| to |*this|.
   // Require
   // 1. base = 2, 3, ..., 9, 11, 12, ..., 36.
   // 2. 0 <= a[i] < BASE_CHG[base].
   // Note that a.back() is allowed to be 0.
   BigInt &move_from_vec_in_base(std::vector<int> &&a, int base){
      while(!a.empty() && a.back()==0){
         a.pop_back();
      }
      if(a.size() >= 2){
         size_t n = a.size();
         while(n != LOWBIT(n)){
            n += LOWBIT(n);
         }
         a.resize(n);
         impl::bigint_change_base(a.data(), n, BASE_CHG[base], BASE);
         while(a.back() == 0){
            a.pop_back();
         }
      }
      *static_cast<std::vector<int>*>(this) = std::move(a);
      return *this;
   }
   BigInt &nontrivial_and_eq(BigInt const &rhs){
      if(sgn() == 1){
         if(rhs.sgn() == 1){
            BigInt p2 = 1ull<<63;
            while(p2.cmp_abs(rhs) <= 0){
               p2 *= p2;
            }
            *this %= p2;
            if(empty()){
               return *this;
            }
            std::vector<int> vl = to_vec_in_base(2), vr = rhs.to_vec_in_base(2);
            size_t n = std::min(vl.size(), vr.size());
            vl.resize(n); vr.resize(n);
            for(size_t i=0; i<n; ++i){
               vl[i] &= vr[i];
            }
            return move_from_vec_in_base(std::move(vl), 2);
         }
         BigInt r2 = ~rhs, p2 = 1ull<<63;
         while(p2.cmp_abs(r2) <= 0){
            p2 *= p2;
         }
         BigInt rem = div(p2).rem;
         std::vector<int> vl = rem.to_vec_in_base(2), vr = r2.to_vec_in_base(2);
         size_t n = std::min(vl.size(), vr.size());
         for(size_t i=0; i<n; ++i){
            vl[i] &= ~vr[i];
         }
         *this -= rem;
         rem.move_from_vec_in_base(std::move(vl), 2);
         return *this += rem;
      }
      into_bitwise_not();
      if(rhs.sgn() == 1){
         BigInt p2 = 1ull<<63;
         while(p2.cmp_abs(rhs) <= 0){
            p2 *= p2;
         }
         *this %= p2;
         std::vector<int> vl = to_vec_in_base(2), vr = rhs.to_vec_in_base(2);
         vl.resize(vr.size());
         for(size_t i=0; i<vl.size(); ++i){
            vl[i] = ~vl[i] & vr[i];
         }
         return move_from_vec_in_base(std::move(vl), 2);
      }
      BigInt r2 = ~rhs, p2 = 1ull<<63;
      while(p2.cmp_abs(r2) <= 0){
         p2 *= p2;
      }
      BigInt rem = div(p2).rem;
      std::vector<int> vl = rem.to_vec_in_base(2), vr = r2.to_vec_in_base(2);
      if(vl.size() < vr.size()){
         vl.resize(vr.size());
      }
      for(size_t i=0; i<vr.size(); ++i){
         vl[i] |= vr[i];
      }
      *this -= rem;
      rem.move_from_vec_in_base(std::move(vl), 2);
      *this += rem;
      return into_bitwise_not();
   }
   BigInt &nontrivial_xor_eq(BigInt const &rhs){
      bool rev = (sgn()==-1) ^ (rhs.sgn()==-1);
      if(sgn() == -1){
         into_bitwise_not();
      }
      BigInt r2 = rhs;
      if(rhs.sgn() == -1){
         r2.into_bitwise_not();
      }
      BigInt p2 = 1ull<<63;
      while(p2.cmp_abs(r2) <= 0){
         p2 *= p2;
      }
      BigInt rem = div(p2).rem;
      std::vector<int> vl = to_vec_in_base(2), vr = r2.to_vec_in_base(2);
      if(vl.size() < vr.size()){
         vl.resize(vr.size());
      }
      for(size_t i=0; i<vr.size(); ++i){
         vl[i] ^= vr[i];
      }
      *this -= rem;
      rem.move_from_vec_in_base(std::move(vl), 2);
      *this += rem;
      if(rev){
         into_bitwise_not();
      }
      return *this;
   }
   friend bool operator==(BigInt const &lhs, BigInt const &rhs) noexcept;
   friend bool operator<(BigInt const &lhs, BigInt const &rhs) noexcept;
   friend BigInt from_vec(std::vector<int> const &a);
   friend BigInt from_vec(std::vector<int> &&a);
   friend BigInt stobi(char const *s, size_t *pos, int base);
   friend BigInt stobi(wchar_t const *s, size_t *pos, int base);
};

inline bool operator==(BigInt const &lhs, BigInt const &rhs) noexcept{
   return static_cast<std::vector<int> const&>(lhs) == static_cast<std::vector<int> const&>(rhs);
}

inline bool operator<(BigInt const &lhs, BigInt const &rhs) noexcept{
   if(lhs.sgn() == -1){
      if(rhs.sgn()>=0 || lhs.size()>rhs.size()){
         return true;
      }
      if(lhs.size()<rhs.size() || lhs.back()>rhs.back()){
         return false;
      }
      if(lhs.back() < rhs.back()){
         return true;
      }
      return std::lexicographical_compare(rhs.crbegin()+1, rhs.crend(), lhs.crbegin()+1, lhs.crend());
   }
   if(lhs.sgn() == 0){
      return rhs.sgn() == 1;
   }
   if(rhs.sgn()<=0 || lhs.size()>rhs.size()){
      return false;
   }
   if(lhs.size() < rhs.size()){
      return true;
   }
   return std::lexicographical_compare(lhs.crbegin(), lhs.crend(), rhs.crbegin(), rhs.crend());
}

inline BigInt from_vec(std::vector<int> const &a){
   if(a.empty()){
      return {};
   }
   auto ok = [BASE=BigInt::BASE](int d){
      return 0<=d && d<=BASE-1;
   };
   if(a.back()<=-BigInt::BASE || a.back()==0 || a.back()>=BigInt::BASE || !std::all_of(a.cbegin(), a.cend()-1, ok)){
      throw std::invalid_argument("no conversion");
   }
   return *static_cast<BigInt const*>(&a);
}

inline BigInt from_vec(std::vector<int> &&a){
   if(a.empty()){
      return {};
   }
   auto ok = [BASE=BigInt::BASE](int d){
      return 0<=d && d<=BASE-1;
   };
   if(a.back()<=-BigInt::BASE || a.back()==0 || a.back()>=BigInt::BASE || !std::all_of(a.cbegin(), a.cend()-1, ok)){
      throw std::invalid_argument("no conversion");
   }
   BigInt res;
   *static_cast<std::vector<int>*>(&res) = std::move(a);
   return res;
}

inline BigInt stobi(char const *s, size_t *pos=nullptr, int base=10){
   if(base<0 || base==1 || base>36){
      throw std::invalid_argument("expect base = 0, 2, 3, ..., 36");
   }
   BigInt res;
   char const *b = s;
   int base_bak = base;
   bool neg;
   while(std::isspace(*b)){
      ++b;
   }
   neg = false;
   for(; ; ++b){
      if(*b == '+'){
         continue;
      }
      if(*b == '-'){
         neg = !neg;
         continue;
      }
      break;
   }
   if(base == 0){
      if(!std::isdigit(*b)){
         goto err_no_conversion;
      }
      if(*b == '0'){
         if((b[1]=='X' || b[1]=='x') && std::isxdigit(b[2])){
            b += 2;
            base = 16;
         }else{
            base = 8;
         }
      }else{
         base = 10;
      }
   }
   if(base == 10){
      if(!std::isdigit(*b)){
         goto err_no_conversion;
      }
      while(*b == '0'){
         ++b;
      }
      char const *e = b;
      while(std::isdigit(*e)){
         ++e;
      }
      if(e > b){
         res.resize((e-b-1)/BigInt::LG10_BASE+1);
         int *di = res.data()+res.size()-1;
         for(int r=static_cast<int>((e-b-1)%BigInt::LG10_BASE)+1; r-->0; ++b){
            *di = 10**di + (*b-'0');
         }
         while(di-- > res.data()){
            for(int j=BigInt::LG10_BASE; j-->0; ++b){
               *di = 10**di + (*b-'0');
            }
         }
         if(neg){
            res.back() = -res.back();
         }
      }
      if(pos != nullptr){
         *pos = e-s;
      }
   }else{
      auto digit = [base](char c){
         if(base > 10){
            if(std::isdigit(c)){
               return c-'0';
            }
            if('A'<=c && c<'A'+(base-10)){
               return c-'A'+10;
            }
            if('a'<=c && c<'a'+(base-10)){
               return c-'a'+10;
            }
            return -1;
         }
         if('0'<=c && c<'0'+base){
            return c-'0';
         }
         return -1;
      };
      if(digit(*b) == -1){
         goto err_no_conversion;
      }
      if(base_bak==16 && b[0]=='0' && (b[1]=='X' || b[1]=='x') && std::isxdigit(b[2])){
         b += 2;
      }
      while(*b == '0'){
         ++b;
      }
      char const *e = b;
      while(digit(*e) != -1){
         ++e;
      }
      if(e > b){
         std::vector<int> a((e-b-1)/BigInt::LGD_BASE[base]+1);
         int *di = a.data()+a.size()-1;
         for(int r=static_cast<int>((e-b-1)%BigInt::LGD_BASE[base])+1; r-->0; ++b){
            *di = base**di + digit(*b);
         }
         while(di-- > a.data()){
            for(int j=BigInt::LGD_BASE[base]; j-->0; ++b){
               *di = base**di + digit(*b);
            }
         }
         res.move_from_vec_in_base(std::move(a), base);
      }
      if(pos != nullptr){
         *pos = e-s;
      }
   }
   return res;
err_no_conversion:
   throw std::invalid_argument("no conversion");
}

inline BigInt stobi(wchar_t const *s, size_t *pos=nullptr, int base=10){
   if(base<0 || base==1 || base>36){
      throw std::invalid_argument("expect base = 0, 2, 3, ..., 36");
   }
   BigInt res;
   wchar_t const *b = s;
   int base_bak = base;
   bool neg;
   while(std::iswspace(*b)){
      ++b;
   }
   neg = false;
   for(; ; ++b){
      if(*b == L'+'){
         continue;
      }
      if(*b == L'-'){
         neg = !neg;
         continue;
      }
      break;
   }
   if(base == 0){
      if(!std::iswdigit(*b)){
         goto err_no_conversion;
      }
      if(*b == L'0'){
         if((b[1]==L'X' || b[1]==L'x') && std::iswxdigit(b[2])){
            b += 2;
            base = 16;
         }else{
            base = 8;
         }
      }else{
         base = 10;
      }
   }
   if(base == 10){
      if(!std::iswdigit(*b)){
         goto err_no_conversion;
      }
      while(*b == '0'){
         ++b;
      }
      wchar_t const *e = b;
      while(std::iswdigit(*e)){
         ++e;
      }
      if(e > b){
         res.resize((e-b-1)/BigInt::LG10_BASE+1);
         int *di = res.data()+res.size()-1;
         for(int r=static_cast<int>((e-b-1)%BigInt::LG10_BASE)+1; r-->0; ++b){
            *di = 10**di + (*b-L'0');
         }
         while(di-- > res.data()){
            for(int j=BigInt::LG10_BASE; j-->0; ++b){
               *di = 10**di + (*b-L'0');
            }
         }
         if(neg){
            res.back() = -res.back();
         }
      }
      if(pos != nullptr){
         *pos = e-s;
      }
   }else{
      auto digit = [base](wchar_t c){
         if(base > 10){
            if(std::iswdigit(c)){
               return c-L'0';
            }
            if(L'A'<=c && c<L'A'+(base-10)){
               return c-L'A'+10;
            }
            if(L'a'<=c && c<L'a'+(base-10)){
               return c-L'a'+10;
            }
            return -1;
         }
         if(L'0'<=c && c<L'0'+base){
            return c-L'0';
         }
         return -1;
      };
      if(digit(*b) == -1){
         goto err_no_conversion;
      }
      if(base_bak==16 && b[0]==L'0' && (b[1]==L'X' || b[1]==L'x') && std::iswxdigit(b[2])){
         b += 2;
      }
      while(*b == L'0'){
         ++b;
      }
      wchar_t const *e = b;
      while(digit(*e) != -1){
         ++e;
      }
      if(e > b){
         std::vector<int> a((e-b-1)/BigInt::LGD_BASE[base]+1);
         int *di = a.data()+a.size()-1;
         for(int r=static_cast<int>((e-b-1)%BigInt::LGD_BASE[base])+1; r-->0; ++b){
            *di = base**di + digit(*b);
         }
         while(di-- > a.data()){
            for(int j=BigInt::LGD_BASE[base]; j-->0; ++b){
               *di = base**di + digit(*b);
            }
         }
         res.move_from_vec_in_base(std::move(a), base);
      }
      if(pos != nullptr){
         *pos = e-s;
      }
   }
   return res;
err_no_conversion:
   throw std::invalid_argument("no conversion");
}

namespace std{
inline string to_string(BigInt const &a){
   return static_cast<string>(a);
}
inline wstring to_wstring(BigInt const &a){
   return static_cast<std::wstring>(a);
}
}

inline BigInt sgn(BigInt const &a) noexcept{
   return a.sgn();
}

inline BigInt abs(BigInt const &a){
   BigInt res = a;
   return res.into_abs();
}

inline BigInt abs(BigInt &&a) noexcept{
   return a.into_abs();
}

inline DivT<BigInt> div(BigInt const &lhs, BigInt const &rhs){
   return lhs.div(rhs);
}

inline bool operator!=(BigInt const &lhs, BigInt const &rhs) noexcept{
   return !(lhs == rhs);
}

inline bool operator<=(BigInt const &lhs, BigInt const &rhs) noexcept{
   return !(rhs < lhs);
}

inline bool operator>(BigInt const &lhs, BigInt const &rhs) noexcept{
   return rhs < lhs;
}

inline bool operator>=(BigInt const &lhs, BigInt const &rhs) noexcept{
   return !(lhs < rhs);
}

#define DEF_BIOP(OP)\
inline BigInt operator OP(BigInt const &lhs, BigInt const &rhs){\
   BigInt res = lhs;\
   return res OP##= rhs;\
}\
inline BigInt operator OP(BigInt &&lhs, BigInt const &rhs){\
   lhs OP##= rhs;\
   return std::move(lhs);\
}
DEF_BIOP(+)
DEF_BIOP(-)
DEF_BIOP(*)
DEF_BIOP(/)
DEF_BIOP(%)
DEF_BIOP(&)
DEF_BIOP(^)
DEF_BIOP(|)
#undef DEF_BIOP
#define DEF_SA(OP)\
inline BigInt operator OP(BigInt const &lhs, size_t rhs){\
   BigInt res = lhs;\
   return res OP##= rhs;\
}\
inline BigInt operator OP(BigInt &&lhs, size_t rhs){\
   lhs OP##= rhs;\
   return std::move(lhs);\
}
DEF_SA(<<)
DEF_SA(>>)
#undef DEF_SA

#undef LOWBIT

#undef U64
#undef I64
#undef U32

#endif // ARBITRARY_PREC_HH
