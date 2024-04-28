#ifndef ARBITRARY_PREC_HH
#define ARBITRARY_PREC_HH

#include<iomanip>
#include<sstream>
#include<stdexcept>
#include<utility>
#include<vector>
#include<cctype>
#include<cstdlib>
#include<cstring>
#include<cwctype>

#define I64 long long
#define U64 unsigned long long

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
   BigInt operator+() const noexcept{
      return *this;
   }
   BigInt &into_inv(){
      if(!empty()){
         data()[size()-1] = -data()[size()-1];
      }
      return *this;
   }
   BigInt operator-() const noexcept{
      auto res = *this;
      return res.into_inv();
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
            return subtract_inv(rhs, false);
         }
         return subtract_from_inv(rhs, false);
      }
      if(rhs.sgn() == -1){
         return into_inv().add_std(rhs, true).into_inv();
      }
      if(cmp_abs(rhs) >= 0){
         return into_inv().subtract_inv(rhs, true).into_inv();
      }
      return into_inv().subtract_from_inv(rhs, true).into_inv();
   }
   BigInt &operator-=(BigInt const &rhs){
      return (into_inv()+=rhs).into_inv();
   }
private:
   BigInt &trunc(){
      while(!empty() && back()==0){
         pop_back();
      }
      return *this;
   }
   int cmp_abs(BigInt const &rhs) const{
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
   BigInt &subtract_inv(BigInt const &rhs, bool this_inv){
      for(size_t i=0; i<rhs.size()-1; ++i){
         data()[i] -= rhs[i];
         if(data()[i] < 0){
            data()[i] += BASE;
            --data()[i+1];
         }
      }
      data()[rhs.size()-1] += this_inv? -rhs[rhs.size()-1]: rhs[rhs.size()-1];
      for(size_t i=rhs.size()-1; data()[i]<0; ++i){
         data()[i] += BASE;
         --data()[i+1];
      }
      return trunc();
   }
   BigInt &subtract_from_inv(BigInt const &rhs, bool this_inv){
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
      data()[size()-1] = (this_inv? rhs[size()-1]: -rhs[size()-1])-data()[size()-1]-borrow;
      return trunc().into_inv();
   }
   friend bool operator==(BigInt const &lhs, BigInt const &rhs);
   friend bool operator<(BigInt const &lhs, BigInt const &rhs);
   friend BigInt stob(char const *s, size_t *pos);
   friend BigInt stob(wchar_t const *s, size_t *pos);
};

inline bool operator==(BigInt const &lhs, BigInt const &rhs){
   if(lhs.size() != rhs.size()){
      return false;
   }
   return std::memcmp(lhs.data(), rhs.data(), lhs.size()*sizeof*lhs.data()) == 0;
}

inline bool operator<(BigInt const &lhs, BigInt const &rhs){
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
         res[res.size()-1] = -res[res.size()-1];
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
         res[res.size()-1] = -res[res.size()-1];
      }
   }
   if(pos != nullptr){
      *pos = e-s;
   }
   return res;
err_no_conversion:
   throw std::invalid_argument("no conversion");
}

inline bool operator!=(BigInt const &lhs, BigInt const &rhs){
   return !(lhs == rhs);
}

inline bool operator<=(BigInt const &lhs, BigInt const &rhs){
   return !(rhs < lhs);
}

inline bool operator>(BigInt const &lhs, BigInt const &rhs){
   return rhs < lhs;
}

inline bool operator>=(BigInt const &lhs, BigInt const &rhs){
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
#undef DEF_BIOP

#undef U64
#undef I64

#endif // ARBITRARY_PREC_HH
