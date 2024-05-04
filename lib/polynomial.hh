#ifndef POLYNOMIAL_HH
#define POLYNOMIAL_HH

#include"convolution.hh"
#include"number_theory.hh"
#include<algorithm>
#include<functional>
#include<initializer_list>
#include<optional>
#include<type_traits>
#include<vector>

#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

namespace impl{
template<typename RingT>
int fps_deg(std::vector<RingT> const &f) noexcept{
   auto it = find_if(f.crbegin(), f.crend(), [](RingT const &ci){return ci != 0;});
   return it.base()-f.cbegin()-1;
}
template<typename RingT>
std::vector<RingT> &fps_add_assign(std::vector<RingT> &f, std::vector<RingT> const &g){
   int d2 = fps_deg(g);
   if(static_cast<int>(f.size()) < d2+1) f.resize(d2+1);
   for(int i=0; i<=d2; ++i){
      f[i] += g[i];
   }
   return f;
}
template<typename RingT>
std::vector<RingT> &fps_sub_assign(std::vector<RingT> &f, std::vector<RingT> const &g){
   int d2 = fps_deg(g);
   if(static_cast<int>(f.size()) < d2+1) f.resize(d2+1);
   for(int i=0; i<=d2; ++i){
      f[i] -= g[i];
   }
   return f;
}
template<typename RingT>
std::vector<RingT> fps_naive_mul(RingT const *f, int deg_f, RingT const *g, int deg_g, int lo, int hi){
   std::vector<RingT> res(hi-lo+1, 0);
   for(int i=std::max(lo-deg_g, 0); i<=std::min(deg_f, hi); ++i){
      for(int j=std::max(lo-i, 0); j<=std::min(deg_g, hi-i); ++j){
         res[i+j-lo] += f[i]*g[j];
      }
   }
   return res;
}
template<typename RingT> struct FpsMulAssign{
   std::vector<RingT> &operator()(std::vector<RingT> &f, std::vector<RingT> const &g){
      int d1 = fps_deg(f), d2 = fps_deg(g);
      if(d1==-1 || d2==-1){
         f.clear(); return f;
      }
      if(d1<32 || d2<32){
         return f = fps_naive_mul(f.data(), d1, g.data(), d2, 0, d1+d2);
      }
      int n = std::max(d1, d2)/2+1;
      auto mid = f.cbegin() + std::min(n, static_cast<int>(f.size()));
      std::vector<RingT> a(f.cbegin(), mid), b(mid, f.cend());
      mid = g.cbegin() + std::min(n, static_cast<int>(g.size()));
      std::vector<RingT> c(g.cbegin(), mid), d(mid, g.cend());
      auto ac = a; (*this)(ac, c);
      auto bd = b; (*this)(bd, d);
      (*this)(fps_add_assign(a, b), fps_add_assign(c, d));
      fps_sub_assign(fps_sub_assign(a, ac), bd);
      bd.insert(bd.cbegin(), 2*n, 0);
      int da = fps_deg(a);
      if(static_cast<int>(bd.size()) < n+da+1){
         bd.resize(n+da+1);
      }
      for(int i=0; i<=da; ++i){
         bd[n+i] += a[i];
      }
      return f = move(fps_add_assign(bd, ac));
   }
};
template<U32 N> struct FpsMulAssign<RingZn<N>>{
   std::vector<RingZn<N>> &operator()(std::vector<RingZn<N>> &f, std::vector<RingZn<N>> const &g){
      int d1 = fps_deg(f), d2 = fps_deg(g);
      if(d1==-1 || d2==-1){
         f.clear(); return f;
      }
      f.resize(d1+d2+1);
      convolution(f.cbegin(), f.cbegin()+d1+1, g.cbegin(), g.cbegin()+d2+1, d1+d2+1, f.begin());
      return f;
   }
};
#define RECIPROCAL_LIFTING(T) do{\
   std::vector<T> two{2};\
   std::vector<T> first_2i(f.cbegin(), f.cbegin()+std::min(2*i, static_cast<int>(f.size())));\
   FpsMulAssign<T> mul_assign;\
   mul_assign(res, fps_sub_assign(two, mul_assign(first_2i, res)));\
   if(static_cast<int>(res.size()) > 2*i) res.resize(2*i);\
}while(0)
template<typename RingT> struct FpsReciprocal{
   std::vector<RingT> operator()(std::vector<RingT> const &f, int n_terms){
      if(n_terms <= 0) return {};
      std::vector<RingT> res{1/f[0]};
      for(int i=1; i<n_terms; i*=2){
         RECIPROCAL_LIFTING(RingT);
      }
      if(static_cast<int>(res.size()) > n_terms) res.resize(n_terms);
      return res;
   }
};
template<U32 N> struct FpsReciprocal<RingZn<N>>{
   std::vector<RingZn<N>> operator()(std::vector<RingZn<N>> const &f, int n_terms){
      if(n_terms <= 0) return {};
      std::vector<RingZn<N>> res{1/f[0]};
      for(int i=1; i<n_terms; i*=2){
         if(primality_test(N) && 4*i<=Ntt<N>::max_size()){
            std::vector<RingZn<N>> g(f.cbegin(), f.cbegin()+std::min(2*i, static_cast<int>(f.size())));
            res.resize(4*i); g.resize(4*i);
            Ntt<N>::transform_in_place(false, res.begin(), 4*i);
            Ntt<N>::transform_in_place(false, g.begin(), 4*i);
            for(int j=0; j<4*i; ++j){
               res[j] *= 2-g[j]*res[j];
            }
            Ntt<N>::transform_in_place(true, res.begin(), 4*i);
            res.resize(2*i); continue;
         }
         RECIPROCAL_LIFTING(RingZn<N>);
      }
      return res;
   }
};
#undef RECIPROCAL_LIFTING
} // namespace impl

template<typename RingT> struct BasicFps: private std::vector<RingT>{
   using Base = std::vector<RingT>;
   using Base::capacity;
   using Base::clear;
   using Base::get_allocator;
   using Base::reserve;
   using Base::resize;
   using Base::shrink_to_fit;
   using Base::size;
   BasicFps() = default;
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>
   BasicFps(INT c): BasicFps(static_cast<RingT>(c)){}
   BasicFps(RingT const &c): BasicFps(c==0? Base{}: Base{c}){}
   BasicFps(std::initializer_list<RingT> coef): BasicFps(Base(coef)){}
   template<typename It> BasicFps(It b, It e): BasicFps(Base(b, e)){}
   explicit BasicFps(Base coef): Base(move(coef)){}
   Base const &as_vec() const noexcept{
      return *this;
   }
   Base into_vec() noexcept{
      Base res;
      res.swap(*this);
      return res;
   }
   int deg() const noexcept{
      return impl::fps_deg(*this);
   }
   explicit operator bool() const noexcept{
      return deg() >= 0;
   }
   RingT coef(int i) const noexcept{
      return i>=static_cast<int>(size())? 0: (*this)[i];
   }
   template<typename T> BasicFps &set_coef(int i, T &&ci){
      if(i >= static_cast<int>(size())) resize(i+1);
      (*this)[i] = std::forward<T>(ci);
      return *this;
   }
   RingT operator()(RingT const &x) const noexcept{
      RingT res = 0;
      for(size_t i=size(); i-->0; ){
         res = x*res + (*this)[i];
      }
      return res;
   }
   BasicFps &operator<<=(int n_terms){
      this->insert(this->cbegin(), n_terms, 0);
      return *this;
   }
   BasicFps operator<<(int n_terms) const{
      auto res = *this;
      return res <<= n_terms;
   }
   BasicFps &operator>>=(int n_terms){
      erase(this->cbegin(), this->cbegin()+std::min(n_terms, static_cast<int>(size())));
      return *this;
   }
   BasicFps operator>>(int n_terms) const{
      if(n_terms < static_cast<int>(size())){
         return BasicFps(this->cbegin()+n_terms, this->cend());
      }
      return {};
   }
   BasicFps &trunc(int n_terms) noexcept{
      if(n_terms < static_cast<int>(size())){
         resize(n_terms);
      }
      return *this;
   }
   BasicFps first_n_terms(int n) const{
      return BasicFps(this->cbegin(), this->cbegin()+std::min(n, static_cast<int>(size())));
   }
   BasicFps rev() const{
      int d = deg();
      std::vector<RingT> res;
      res.reserve(d+1);
      reverse_copy(this->cbegin(), this->cbegin()+(d+1), back_inserter(res));
      return BasicFps(move(res));
   }
   BasicFps operator+() const{
      return *this;
   }
   BasicFps operator-() const{
      auto res = *this;
      transform(res.begin(), res.end(), res.begin(), std::negate<RingT>());
      return res;
   }
   BasicFps &operator+=(BasicFps const &rhs){
      impl::fps_add_assign(*this, rhs);
      return *this;
   }
   BasicFps &operator-=(BasicFps const &rhs){
      impl::fps_sub_assign(*this, rhs);
      return *this;
   }
   BasicFps &operator*=(BasicFps const &rhs){
      impl::FpsMulAssign<RingT>()(*this, rhs);
      return *this;
   }
   BasicFps reciprocal(int n_terms) const{
      return BasicFps(impl::FpsReciprocal<RingT>()(*this, n_terms));
   }
   BasicFps &operator/=(BasicFps const &rhs){
      int n = deg(), d = rhs.deg();
      if(d<32 || n-d<32){
         resize(n+1);
         return *this = BasicFps(naive_division(std::move(*this), rhs.data(), d).first);
      }
      auto rf = rev(), rg = rhs.rev();
      *this = rf.trunc(n-d+1)*rg.reciprocal(n-d+1);
      resize(n-d+1);
      reverse(this->begin(), this->end());
      return *this;
   }
   BasicFps &operator%=(BasicFps const &rhs){
      int n = deg(), d = rhs.deg();
      if(d<32 || n-d<32){
         resize(n+1);
         return *this = BasicFps(naive_division(std::move(*this), rhs.data(), d).second);
      }
      auto q = *this;
      q /= rhs; q *= rhs;
      trunc(d); q.trunc(d);
      return *this -= q;
   }
   BasicFps derivative() const{
      BasicFps res;
      if(int d=deg(); d>0){
         res.resize(d);
         for(int i=1; i<=d; ++i){
            res[i-1] = i*(*this)[i];
         }
      }
      return res;
   }
   BasicFps integral() const{
      BasicFps res;
      if(int d=deg(); d>=0){
         res.resize(d+2);
         for(int i=0; i<=d; ++i){
            res[i+1] = (*this)[i]/(i+1);
         }
      }
      return res;
   }
#define DEF_BIOP(OP)\
   friend BasicFps operator OP(BasicFps lhs, BasicFps const &rhs){\
      return lhs OP##= rhs;\
   }
   DEF_BIOP(+)
   DEF_BIOP(-)
   DEF_BIOP(*)
   DEF_BIOP(/)
   DEF_BIOP(%)
#undef DEF_BIOP
   friend bool operator==(BasicFps const &lhs, BasicFps const &rhs){
      int d = lhs.deg();
      return rhs.deg()==d && equal(lhs.cbegin(), lhs.cbegin()+(d+1), rhs.cbegin());
   }
   friend bool operator!=(BasicFps const &lhs, BasicFps const &rhs){
      return !(lhs == rhs);
   }
private:
   static std::pair<std::vector<RingT>, std::vector<RingT>> naive_division(std::vector<RingT> f, RingT const *g, int deg_g){
      int deg_f = static_cast<int>(f.size())-1;
      if(deg_f < deg_g){
         return {{}, move(f)};
      }
      std::vector<RingT> q(deg_f-deg_g+1);
      for(int i=deg_f-deg_g; i>=0; --i){
         RingT qi = f[i+deg_g]/g[deg_g];
         f[i+deg_g] = 0;
         for(int j=0; j<deg_g; ++j){
            f[i+j] -= qi*g[j];
         }
         q[i] = qi;
      }
      f.resize(deg_g);
      return {move(q), move(f)};
   }
};

template<U32 N> using FpsMod = BasicFps<RingZn<N>>;

template<typename RingT>
BasicFps<RingT> log_fps(BasicFps<RingT> const &f, int n_terms){
   // assert(f.coef(0) == 1);
   if(n_terms <= 0) return {};
   return (f.derivative()*f.reciprocal(n_terms-1)).trunc(n_terms-1).integral();
}

template<typename RingT>
BasicFps<RingT> exp_fps(BasicFps<RingT> const &f, int n_terms){
   // assert(f.coef(0) == 0);
   if(n_terms <= 0) return {};
   if(n_terms == 1) return 1;
   auto g = exp_fps(f.first_n_terms((n_terms+1)/2), (n_terms+1)/2);
   return (g*(1-log_fps(g, n_terms)+f)).trunc(n_terms);
}

template<typename RingT>
BasicFps<RingT> pow_fps(BasicFps<RingT> const &f, I64 n, int n_terms){
   // assert(n >= 0);
   if(n_terms <= 0) return {};
   if(n == 0) return 1;
   int d = f.deg();
   if(d == -1) return 0;
   int low = 0;
   while(!f.coef(low)) ++low;
   if((n_terms-1)/n < low) return 0;
   int m = n_terms-n*low;
   auto res = f >> low;
   res /= f.coef(low);
   res = exp_fps(n*log_fps(res, m), m);
   res *= pow(f.coef(low), n);
   return res <<= n*low;
}

template<U32 P, typename = std::enable_if_t<primality_test(P) && P!=2>>
std::optional<FpsMod<P>> sqrt_fps(FpsMod<P> const &f, int n_terms){
   int d = f.deg();
   if(d == -1) return 0;
   int low = 0;
   while(!f.coef(low)) ++low;
   I64 c0;
   if(low % 2 || (c0 = sqrt_modp(static_cast<U32>(f.coef(low)), P)) == -1){
      return {};
   }
   FpsMod<P> g = f>>low, res = c0;
   for(int i=1; i<n_terms-low/2; i*=2){
      res = (res+g.first_n_terms(2*i)*res.reciprocal(2*i))/2;
   }
   return std::move(res<<=low/2);
}

namespace impl{
template<typename RingT>
BasicFps<RingT> naive_comp_fps_dac(BasicFps<RingT> const &f, int n_terms, BasicFps<RingT> const *p, int l, int r){
   if(l == r) return f.coef(l);
   int m = l+(r-l)/2;
   auto lo = naive_comp_fps_dac(f, n_terms, p+1, l, m);
   auto hi = naive_comp_fps_dac(f, n_terms, p+1, m+1, r);
   return (lo+*p*hi).trunc(n_terms);
}

template<typename RingT>
BasicFps<RingT> naive_comp_fps(BasicFps<RingT> const &f, BasicFps<RingT> const &g, int n_terms){
   if(n_terms <= 0) return 0;
   if(n_terms == 1) return f.coef(0);
   int n = n_terms;
   while(n != (n&-n)) n += n&-n;
   std::vector<BasicFps<RingT>> powg{g};
   for(int i=2; i<n; i*=2){
      powg.push_back((powg.back()*powg.back()).trunc(n_terms));
   }
   reverse(powg.begin(), powg.end());
   return naive_comp_fps_dac(f, n_terms, powg.data(), 0, n-1);
}
} // namespace impl

template<typename RingT>
BasicFps<RingT> composite_fps(BasicFps<RingT> const &f, BasicFps<RingT> const &g, int n_terms){
   // assert(g.coef(0));
   if(n_terms <= 32){
      return impl::naive_comp_fps(f, g, n_terms);
   }
   int m = 5*sqrt(n_terms/log2(n_terms));
   auto gl = g.first_n_terms(m), gh = g-gl;
   if(!gl){
      BasicFps<RingT> res = f.coef(0), p = 1;
      for(int i=1; i<=n_terms/m; ++i){
         (p *= gh).trunc(n_terms);
         res += f.coef(i)*p;
      }
      return res;
   }
   auto fogl = impl::naive_comp_fps(f, gl, n_terms), h = gl.derivative();
   int shift = 0;
   while(h.coef(shift) == 0) ++shift;
   h = (h>>shift).reciprocal(n_terms-shift);
   BasicFps<RingT> res = fogl, p = 1;
   for(int i=1; i<=n_terms/m; ++i){
      h.trunc(n_terms-i*shift);
      fogl = ((fogl.derivative()>>shift)*h).trunc(n_terms-i*shift);
      (p *= gh/i).trunc(n_terms);
      res += (fogl*p).trunc(n_terms);
   }
   return res;
}

#undef U64
#undef I64
#undef U32

#endif // POLYNOMIAL_HH
