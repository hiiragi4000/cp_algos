#ifndef POLYNOMIAL_HH
#define POLYNOMIAL_HH

#include"convolution.hh"
#include<initializer_list>

#define I64 long long

template<typename NttT> struct PolynomialNtt;
template<typename NttT>
bool operator==(PolynomialNtt<NttT> const&, PolynomialNtt<NttT> const&);

template<typename NttT> struct PolynomialNtt: private std::vector<int>{
   static constexpr int max_size() noexcept{
      return NttT::max_size();
   }
   static constexpr int mod() noexcept{
      return NttT::prime();
   }
   using Base = std::vector<int>;
   using Base::capacity;
   using Base::clear;
   using Base::get_allocator;
   using Base::reserve;
   using Base::shrink_to_fit;
   using Base::size;
   PolynomialNtt() = default;
   PolynomialNtt(int c): PolynomialNtt(Base{c}){}
   PolynomialNtt(std::initializer_list<int> coef): PolynomialNtt(Base(coef)){}
   template<typename It> PolynomialNtt(It b, It e): PolynomialNtt(Base(b, e)){}
   explicit PolynomialNtt(Base coef): Base(move(coef)){
      for(size_t i=0; i<size(); ++i){
         data()[i] = ((I64)data()[i]%mod() + mod()) % mod();
      }
   }
   PolynomialNtt &resize(size_t n){
      Base::resize(n);
      return *this;
   }
   int deg() const noexcept{
      auto it = find_if(crbegin(), crend(), [](int ci){return ci;});
      return it.base()-cbegin()-1;
   }
   explicit operator bool() const noexcept{
      return deg() >= 0;
   }
   int coef(int i) const noexcept{
      return i>=(int)size()? 0: data()[i];
   }
   template<typename T> PolynomialNtt &set_coef(int i, T ci){
      if(i >= (int)size()) resize(i+1);
      data()[i] = (int)RingZn(mod(), ci);
      return *this;
   }
   int operator()(int x) const noexcept{
      RingZn res(mod(), 0);
      for(size_t i=size(); i-->0; ){
         res = x*res + data()[i];
      }
      return (int)res;
   }
   PolynomialNtt &truncate() noexcept{
      return resize(deg()+1);
   }
   PolynomialNtt rev() const{
      int d = deg();
      vector<int> res(d+1);
      reverse_copy(cbegin(), cbegin()+(d+1), res.begin());
      return PolynomialNtt(move(res));
   }
   PolynomialNtt operator+() const{
      return *this;
   }
   PolynomialNtt operator-() const{
      auto res = *this;
      for(size_t i=0; i<size(); ++i) if(res[i]){
         res[i] = mod() - res[i];
      }
      return res;
   }
   PolynomialNtt &operator+=(PolynomialNtt const &rhs){
      int d2 = rhs.deg();
      if((int)size() < d2+1) resize(d2+1);
      for(int i=0; i<=d2; ++i){
         data()[i] = ((I64)data()[i] + rhs[i]) % mod();
      }
      return *this;
   }
   PolynomialNtt &operator-=(PolynomialNtt const &rhs){
      int d2 = rhs.deg();
      if((int)size() < d2+1) resize(d2+1);
      for(int i=0; i<=d2; ++i){
         data()[i] = ((I64)data()[i] + mod() - rhs[i]) % mod();
      }
      return *this;
   }
   PolynomialNtt &operator*=(PolynomialNtt const &rhs){
      int d1 = deg(), d2 = rhs.deg(), n = d1+d2+1;
      if(d1==-1 || d2==-1){
         clear(); return *this;
      }
      if(d1<32 || d2<32){
         return *this = PolynomialNtt(naive_multiply(data(), d1, rhs.data(), d2, 0, d1+d2));
      }
      while((n&-n) != n) n += n&-n;
      resize(n);
      NttT::transform_in_place(false, begin(), n);
      if(this == &rhs){
         for(int i=0; i<n; ++i){
            data()[i] = (I64)data()[i] * data()[i] % mod();
         }
      }else{
         auto b = NttT::transform(false, rhs.cbegin(), rhs.cend(), n);
         for(int i=0; i<n; ++i){
            data()[i] = (I64)data()[i] * b[i] % mod();
         }
      }
      NttT::transform_in_place(true, begin(), n);
      return *this;
   }
   PolynomialNtt &operator/=(PolynomialNtt const &rhs){
      auto rf = rev(), rg = rhs.rev();
      int n = (int)rf.size()-1, d = (int)rg.size()-1;
      if(n < d){
         clear(); return *this;
      }
      *this = rf.resize(n-d+1)*rg.reciprocal(n-d+1).resize(n-d+1);
      reverse(begin(), end());
      return *this;
   }
   PolynomialNtt &operator%=(PolynomialNtt const &rhs){
      int n = deg(), d = rhs.deg();
      if(n < d) return *this;
      auto q = *this;
      q /= rhs; q *= rhs;
      resize(d); q.resize(d);
      return *this -= q;
   }
   PolynomialNtt derivative() const{
      PolynomialNtt res;
      if(int d=deg(); d>0){
         res.resize(d);
         for(int i=1; i<=d; ++i){
            res[i-1] = (I64)i * data()[i] % mod();
         }
      }
      return res;
   }
   PolynomialNtt integral() const{
      PolynomialNtt res;
      if(int d=deg(); d>=0){
         res.resize(d+2);
         for(int i=0; i<=d; ++i){
            res[i+1] = inv_mod(i+1, mod()) * data()[i] % mod();
         }
      }
      return res;
   }
   PolynomialNtt reciprocal(int n_terms) const{
      if(n_terms <= 0) return {};
      int d = deg(), n = n_terms;
      while((n&-n) != n) n += n&-n;
      std::vector<int> res(n);
      res[0] = inv_mod(coef(0), mod());
      for(int l=1; l<n; l<<=1){
         std::vector<int> g2;
         if(l <= 32){
            PolynomialNtt h(naive_multiply(data(), std::min(d, l-1), res.data(), l-1, l, 2*l-1));
            if(d >= l){
               h += PolynomialNtt(naive_multiply(data()+l, std::min(d-l, l-1), res.data(), l-1, 0, l-1));
            }
            g2 = naive_multiply(res.data(), l-1, h.data(), l-1, 0, l-1);
            for(int i=0; i<l; ++i) if(g2[i]){
               g2[i] = mod()-g2[i];
            }
         }else{
            auto f1 = NttT::transform(false, cbegin(), cbegin()+std::min(d+1, l), 2*l);
            auto g1 = NttT::transform(false, res.cbegin(), res.cbegin()+l, 2*l);
            std::vector<int> h(2*l);
            for(int i=0; i<2*l; ++i){
               h[i] = ((I64)f1[i]*g1[i] + mod() - 1) % mod();
               if(i&1 && h[i]) h[i] = mod()-h[i];
            }
            if(d >= l){
               auto f2 = NttT::transform(false, cbegin()+l, cbegin()+std::min(d+1, 2*l), 2*l);
               for(int i=0; i<2*l; ++i){
                  f2[i] = (I64)f2[i] * g1[i] % mod();
               }
               NttT::transform_in_place(true, f2.begin(), 2*l);
               fill_n(f2.begin()+l, l, 0);
               NttT::transform_in_place(false, f2.begin(), 2*l);
               for(int i=0; i<2*l; ++i){
                  h[i] = ((I64)h[i] + f2[i]) % mod();
               }
            }
            g2.resize(2*l);
            for(int i=0; i<2*l; ++i){
               g2[i] = (I64)(mod()-g1[i]) * h[i] % mod();
            }
            NttT::transform_in_place(true, g2.begin(), 2*l);
         }
         copy_n(g2.cbegin(), l, res.begin()+l);
      }
      res.resize(n_terms);
      return PolynomialNtt(move(res));
   }
   friend bool operator==<>(PolynomialNtt const &lhs, PolynomialNtt const &rhs);
private:
   static std::vector<int> naive_multiply(int const *f, int deg_f, int const *g, int deg_g, int lo, int hi){
      std::vector<int> res(hi-lo+1);
      for(int i=std::max(lo-deg_g, 0); i<=std::min(deg_f, hi); ++i){
         for(int j=std::max(lo-i, 0); j<=std::min(deg_g, hi-i); ++j){
            res[i+j-lo] = (res[i+j-lo] + (I64)f[i]*g[j]) % mod();
         }
      }
      return res;
   }
};

#define DEF_BIOP(op) template<typename NttT>\
PolynomialNtt<NttT> operator op(PolynomialNtt<NttT> lhs, PolynomialNtt<NttT> const &rhs){\
   return lhs op##= rhs;\
}
DEF_BIOP(+)
DEF_BIOP(-)
DEF_BIOP(*)
DEF_BIOP(/)
DEF_BIOP(%)
#undef DEF_BIOP
template<typename NttT>
bool operator==(PolynomialNtt<NttT> const &lhs, PolynomialNtt<NttT> const &rhs){
   int d = lhs.deg();
   return rhs.deg()==d && equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
}
template<typename NttT>
bool operator!=(PolynomialNtt<NttT> const &lhs, PolynomialNtt<NttT> const &rhs){
   return !(lhs == rhs);
}

#undef I64

#endif