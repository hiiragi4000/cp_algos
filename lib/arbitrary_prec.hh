#ifndef ARBITRARY_PREC_HH
#define ARBITRARY_PREC_HH

#include"convolution.hh"
#include<algorithm>
#include<iomanip>
#include<sstream>
#include<stdexcept>
#include<utility>
#include<vector>
#include<cctype>
#include<csignal>
#include<cstdlib>
#include<cstring>
#include<cwctype>

#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

template<typename T> struct DivT{
   T quot, rem;
};

struct BigInt: private std::vector<int>{
   constexpr static int BASE = 100'000'000;
   constexpr static int LG10_BASE = 8;
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
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>
   BigInt(INT a): BigInt(static_cast<I64>(a)){}
   explicit operator bool() const noexcept{
      return !empty();
   }
   explicit operator std::string() const{
      std::ostringstream oss;
      if(empty()){
         oss << '0';
      }else{
         oss << back() << std::setfill('0');
         for(size_t i=size()-1; i-->0; ){
            oss << std::setw(LG10_BASE) << data()[i];
         }
      }
      return oss.str();
   }
   explicit operator std::wstring() const{
      std::wostringstream oss;
      if(empty()){
         oss << L'0';
      }else{
         oss << back() << std::setfill(L'0');
         for(size_t i=size()-1; i-->0; ){
            oss << std::setw(LG10_BASE) << data()[i];
         }
      }
      return oss.str();
   }
   int sgn() const noexcept{
      return empty()? 0: back()>0? 1: -1;
   }
   BigInt operator+() const{
      return *this;
   }
   BigInt &into_inv() noexcept{
      if(!empty()){
         back() = -back();
      }
      return *this;
   }
   BigInt operator-() const{
      auto res = *this;
      return res.into_inv();
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
               into_inv();
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
            return sub_inv(rhs, false);
         }
         return sub_from_inv(rhs, false);
      }
      if(rhs.sgn() == -1){
         return into_inv().add_std(rhs, true).into_inv();
      }
      if(cmp_abs(rhs) >= 0){
         return into_inv().sub_inv(rhs, true).into_inv();
      }
      return into_inv().sub_from_inv(rhs, true).into_inv();
   }
   BigInt &operator-=(BigInt const &rhs){
      return (into_inv()+=rhs).into_inv();
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
         into_inv();
      }
      if(size()<=THRES_MUL || rhs.size()<=THRES_MUL){
         naive_mul(rhs);
      }else{
         fft_mul(rhs);
      }
      if(sign == -1){
         into_inv();
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
         qr.quot.into_inv();
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
   BigInt &add_std(BigInt const &rhs, bool this_inv){
      if(size() < rhs.size()){
         resize(rhs.size());
      }
      bool carry = false;
      for(size_t i=0; i<rhs.size(); ++i){
         data()[i] += (this_inv && i==rhs.size()-1? -rhs[i]: rhs[i]) + carry;
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
   BigInt &sub_inv(BigInt const &rhs, bool this_inv){
      for(size_t i=0; i<rhs.size()-1; ++i){
         data()[i] -= rhs[i];
         if(data()[i] < 0){
            data()[i] += BASE;
            --data()[i+1];
         }
      }
      data()[rhs.size()-1] += this_inv? -rhs.back(): rhs.back();
      for(size_t i=rhs.size()-1; data()[i]<0; ++i){
         data()[i] += BASE;
         --data()[i+1];
      }
      return trunc();
   }
   BigInt &sub_from_inv(BigInt const &rhs, bool this_inv){
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
      back() = (this_inv? rhs.back(): -rhs.back())-back()-borrow;
      return trunc().into_inv();
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
#define LOWBIT(N) ((N)&(~(N)+1))
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
         U64 q = static_cast<U64>(round((z[i].imag()/2-x[i])/NttU32::prime()));
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
#undef LOWBIT
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
            a.into_inv(); r = -r;
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
         while(a_cmp_at(b, i) >= 0){
            a_sub_at(b, i);
            ++qi;
         }
         q[i] = qi;
         h = a[i+b.size()-1];
      }
      q.trunc(); // directly checking q.back() doesn't work. (why?)
      a.trunc().base_div_eq(c);
      if(sgn() == -1){
         q.into_inv(); a.into_inv();
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
         std::copy_n(b.crbegin(), std::min(y.size(), b.size()), y.rbegin());
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
      if(r.size()>b.size() || (r.size()==b.size() && !std::lexicographical_compare(r.crbegin(), r.crend(), b.crbegin(), b.crend()))){
         ++q; r -= b;
      }
      r.base_div_eq(c);
      if(sgn() == -1){
         q.into_inv(); r.into_inv();
      }
      return {std::move(q), std::move(r)};
   }
   friend bool operator==(BigInt const &lhs, BigInt const &rhs) noexcept;
   friend bool operator<(BigInt const &lhs, BigInt const &rhs) noexcept;
   friend BigInt stob(char const *s, size_t *pos);
   friend BigInt stob(wchar_t const *s, size_t *pos);
};

inline bool operator==(BigInt const &lhs, BigInt const &rhs) noexcept{
   if(lhs.size() != rhs.size()){
      return false;
   }
   return std::memcmp(lhs.data(), rhs.data(), lhs.size()*sizeof*lhs.data()) == 0;
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
      return lexicographical_compare(rhs.crbegin()+1, rhs.crend(), lhs.crbegin()+1, lhs.crend());
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
   return lexicographical_compare(lhs.crbegin(), lhs.crend(), rhs.crbegin(), rhs.crend());
}

inline BigInt stob(char const *s, size_t *pos=nullptr){
   BigInt res;
   char const *b = s, *e;
   bool neg;
   int *di;
   while(std::isspace(*b)){
      ++b;
   }
   if(*b == '\0'){
      goto err_no_conversion;
   }
   neg = false;
   for(;; ++b){
      if(*b == '+'){
         continue;
      }
      if(*b == '-'){
         neg = !neg;
         continue;
      }
      break;
   }
   if(!std::isdigit(*b)){
      goto err_no_conversion;
   }
   while(*b == '0'){
      ++b;
   }
   e = b;
   while(std::isdigit(*e)){
      ++e;
   }
   if(e > b){
      res.resize((e-b-1)/BigInt::LG10_BASE+1);
      di = res.data()+res.size()-1;
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
   return res;
err_no_conversion:
   throw std::invalid_argument("no conversion");
}

inline BigInt stob(wchar_t const *s, size_t *pos=nullptr){
   BigInt res;
   wchar_t const *b = s, *e;
   bool neg;
   int *di;
   while(std::iswspace(*b)){
      ++b;
   }
   if(*b == L'\0'){
      goto err_no_conversion;
   }
   neg = false;
   for(;; ++b){
      if(*b == L'+'){
         continue;
      }
      if(*b == L'-'){
         neg = !neg;
         continue;
      }
      break;
   }
   if(!std::iswdigit(*b)){
      goto err_no_conversion;
   }
   while(*b == L'0'){
      ++b;
   }
   e = b;
   while(std::iswdigit(*e)){
      ++e;
   }
   if(e > b){
      res.resize((e-b-1)/BigInt::LG10_BASE+1);
      di = res.data()+res.size()-1;
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
   return res;
err_no_conversion:
   throw std::invalid_argument("no conversion");
}

inline BigInt sgn(BigInt const &a){
   return a.sgn();
}

inline BigInt abs(BigInt const &a){
   BigInt res = a;
   return res.into_abs();
}

inline BigInt abs(BigInt &&a){
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
   return lhs OP##= rhs;\
}
DEF_BIOP(+)
DEF_BIOP(-)
DEF_BIOP(*)
DEF_BIOP(/)
DEF_BIOP(%)
#undef DEF_BIOP

#undef U64
#undef I64
#undef U32

#endif // ARBITRARY_PREC_HH
