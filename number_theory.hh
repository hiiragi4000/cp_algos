#ifndef NUMBER_THEORY_HH
#define NUMBER_THEORY_HH

#include<array>
#include<functional>
#include<numeric>
#include<unordered_map>
#include<utility>
#include<vector>
#include<cmath>

#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

struct PrimeEnumerator{
   static std::pair<std::vector<int>, std::function<bool(int)>> leq(int n){
      // assert(n >= 0);
      int b = n/30+1;
      std::vector<bool> isp(b<<3, true);
      std::vector<int> primes{2, 3, 5};
      primes.reserve(UB*30*b/std::log(30*b));
      isp[0] = false;
      for(int i=0; i<b; ++i){
         for(int j=0; j<8; ++j){
            I64 d = 30ll*i + r[j];
            if(d == 1) continue;
            if(isp[i<<3|j]) primes.push_back(d);
            for(int k=3; ; ++k){
               I64 pd = primes[k]*d;
               if(pd >= 30*b) break;
               isp[pd/30<<3|rinv[pd%30]] = false;
               if(d%primes[k] == 0) break;
            }
         }
      }
      while(!primes.empty() && primes.back()>n){
         primes.pop_back();
      }
      auto is_prime = [isp=move(isp)](int i){
         if(i <= 1) return false;
         if(i==2 || i==3 || i==5) return true;
         if(i%2==0 || i%3==0 || i%5==0) return false;
         return (bool)isp[i/30<<3|rinv[i%30]];
      };
      return {primes, is_prime};
   }
private:
   static constexpr int r[] = {1, 7, 11, 13, 17, 19, 23, 29};
   static constexpr signed char rinv[] = {
      -1,  0, -1, -1, -1, -1, -1,  1, -1, -1,
      -1,  2, -1,  3, -1, -1, -1,  4, -1,  5,
      -1, -1, -1,  6, -1, -1, -1, -1, -1,  7
   };
   static constexpr double UB = 1.25506;
};

inline std::vector<int> phi_table_leq(int n){
   // assert(n >= -1);
   std::vector<int> phi(n+1u);
   std::vector<int> primes;
   iota(phi.begin(), phi.end(), 0u);
   for(I64 i=2; i<=n; ++i){
      if(phi[i] == i){
         --phi[i];
         primes.push_back(i);
      }
      for(int p: primes){
         if(p*i > n) break;
         if(i%p){
            phi[p*i] = (p-1)*phi[i];
         }else{
            phi[p*i] = p*phi[i];
            break;
         }
      }
   }
   return phi;
}

inline std::vector<signed char> mu_table_leq(int n){
   // assert(n >= 1);
   std::vector<signed char> mu(n+1u, -2);
   std::vector<int> primes;
   mu[1] = 1;
   for(I64 i=2; i<=n; ++i){
      if(mu[i] == -2){
         mu[i] = -1;
         primes.push_back(i);
      }
      for(int p: primes){
         if(p*i > n) break;
         if(i%p){
            mu[p*i] = -mu[i];
         }else{
            mu[p*i] = 0;
            break;
         }
      }
   }
   return mu;
}

constexpr std::array<I64, 2> extgcd(I64 a, I64 b){
   // assert(a>=0 && b>=0);
   if(b == 0) return {1, 0};
   auto [x, y] = extgcd(b, a%b);
   return {y, x-a/b*y};
}

constexpr I64 inv_mod(I64 a, I64 mod){
   // assert(a>0 && mod>0 && std::gcd(a, mod)==1);
   auto [x, y] = extgcd(a, mod);
   return (x%=mod)<0? x+mod: x;
}

constexpr U32 pow_mod(I64 a, I64 n, U32 mod){
   // assert(mod > 0);
   if(mod == 1) return 0;
   a = (a%mod+mod)%mod;
   if(n < 0) a = inv_mod(a, mod);
   U64 res = 1, t = a;
   while(n){
      if(n%2) res = res*t%mod;
      if(n/=2) t = t*t%mod;
   }
   return res;
}

constexpr bool miller_rabin(U32 n, U32 witness){
   // assert(n>=3 && n%2==0);
   if((witness%=n) == 0) return true;
   U32 d = n-1;
   int r = 0;
   while(d%2 == 0){
      d>>=1, ++r;
   }
   U32 x = pow_mod(witness, d, n);
   if(x==1 || x==n-1) return true;
   for(int cnt=r-1; cnt-->0; ){
      x = (U64)x*x%n;
      if(x == n-1) return true;
   }
   return false;
}

constexpr bool primality_test(U32 n) noexcept{
   if(n <= 1) return false;
   if(n == 2) return true;
   if(n%2 == 0) return false;
   return miller_rabin(n, 2) && (n < 2047 || miller_rabin(n, 7) && miller_rabin(n, 61));
}

constexpr bool primality_test(int n) noexcept{
   return n<0? false: primality_test((U32)n);
}

constexpr int jacobi_symbol(I64 a, I64 b){
   // assert(b>=1 && b%2==1);
   if((a%=b) < 0) a += b;
   bool ok = true;
   while(a){
      bool n2 = false;
      while(a%2 == 0){
         a>>=1, n2=!n2;
      }
      if(n2 && (b%8==3 || b%8==5)){
         ok = !ok;
      }
      if(a%4==3 && b%4==3){
         ok = !ok;
      }
      int t = b%a; b = a; a = t;
   }
   return b>1? 0: ok? 1: -1;
}

constexpr int sqrt_modp(int a, int p){
   // assert(primality_test(p));
   a = ((I64)a%p+p)%p;
   if(p == 2) return a;
   if(a == 0) return 0;
   if(jacobi_symbol(a, p) == -1) return -1;
   if(p%4 == 3) return pow_mod(a, (p+1ll)/4, p);
   I64 k = 1;
   while(jacobi_symbol(k*k-a, p) == 1) ++k;
   if(k*k%p == a) return k;
   I64 w = ((k*k-a)%p+p)%p;
   I64 rc = 1, rx = 0, tc = k, tx = 1;
   for(I64 i=(p+1ll)/2; i>0; ){
      if(i%2){
         I64 rc2 = (rc*tc + rx*tx%p*w)%p;
         rx = (rc*tx + rx*tc)%p;
         rc = rc2;
      }
      if(i/=2){
         I64 tc2 = (tc*tc + tx*tx%p*w)%p;
         tx = 2*tc*tx%p;
         tc = tc2;
      }
   }
   if(rc > p-rc) rc = p-rc;
   return rc;
}

constexpr bool primitive_root_modp_check(int r, int p){
   // assert(primality_test(p));
   if(r%p == 0) return false;
   int d = p-1;
   for(I64 q=2; q*q<=d; ++q) if(d%q == 0){
      do d/=q; while(d%q==0);
      if(pow_mod(r, (p-1)/q, p) == 1) return false;
   }
   if(d>1 && pow_mod(r, (p-1)/d, p) == 1){
      return false;
   }
   return true;
}

constexpr int primitive_root_modp(int p){
   // assert(primality_test(p));
   for(int i=1; ; ++i){
      if(primitive_root_modp_check(i, p)){
         return i;
      }
   }
}

inline int discrete_log(int base, int power, int mod){
   // assert(mod > 0);
   if((base%=mod) < 0) base += mod;
   if((power%=mod) < 0) power += mod;
   if(mod==1 || power==1) return 0;
   int s = 1, d = 1, res = 0;
   while(1){
      s = (I64)base*s%mod;
      ++res;
      if(s == power) return res;
      int d2 = std::gcd(s, mod);
      if(d2 == d){
         if(power%d) return -1;
         s /= d; power /= d; mod /= d;
         break;
      }
      d = d2;
   }
   int m = std::ceil(std::sqrt(mod));
   std::unordered_map<int, int> bs(m);
   for(int i=0; i<m; ++i){
      bs.emplace(s, i);
      s = (I64)base*s%mod;
   }
   int gs = pow_mod(base, -m, mod);
   for(int i=0; i<m; ++i){
      if(auto it=bs.find(power); it!=bs.end()){
         return res + i*m + it->second;
      }
      power = (I64)gs*power%mod;
   }
   return -1;
}

template<U32 n> struct RingZn{
   RingZn() = default;
   template<typename INT> constexpr RingZn(INT a) noexcept: a((a%=(I64)n)<0? (I64)a+n: a){}
   constexpr U32 mod() const noexcept{
      return n;
   }
   template<typename INT> constexpr explicit operator INT() const noexcept{
      return a;
   }
   constexpr RingZn operator+() const noexcept{
      return *this;
   }
   constexpr RingZn operator-() const noexcept{
      return a? n-a: 0;
   }
   constexpr RingZn &operator+=(RingZn rhs) noexcept{
      a += a<n-rhs.a? rhs.a: n-rhs.a;
      return *this;
   }
   constexpr RingZn &operator-=(RingZn rhs) noexcept{
      a = a<rhs.a? a+(n-rhs.a): a-rhs.a;
      return *this;
   }
   constexpr RingZn &operator*=(RingZn rhs) noexcept{
      a = (U64)a * rhs.a % n;
      return *this;
   }
   constexpr RingZn &operator/=(RingZn rhs) noexcept{
      a = (U64)a * inv_mod(rhs.a, n) % n;
      return *this;
   }
   constexpr RingZn &operator++() noexcept{
      a = a==n-1? 0: a+1;
      return *this;
   }
   constexpr RingZn operator++(int) noexcept{
      RingZn res = *this; ++*this;
      return res;
   }
   constexpr RingZn &operator--() noexcept{
      a = a? a-1: n-1;
      return *this;
   }
   constexpr RingZn operator--(int) noexcept{
      RingZn res = *this; --*this;
      return res;
   }
private:
   U32 a = 0;
};
#define DEF_BIOP(op)\
template<U32 n> constexpr RingZn<n> operator op(RingZn<n> lhs, RingZn<n> rhs) noexcept{\
   return lhs op##= rhs;\
}
DEF_BIOP(+)
DEF_BIOP(-)
DEF_BIOP(*)
DEF_BIOP(/)
#undef DEF_BIOP
template<U32 n> constexpr bool operator==(RingZn<n> lhs, RingZn<n> rhs) noexcept{
   return (U32)lhs == (U32)rhs;
}
template<U32 n> constexpr bool operator!=(RingZn<n> lhs, RingZn<n> rhs) noexcept{
   return !(lhs == rhs);
}

struct RingZnDyn{
   RingZnDyn() = default;
   constexpr RingZnDyn(int t) noexcept: t(t){}
   template<typename INT> constexpr RingZnDyn(U32 n, INT a) noexcept: n(n), a((a%=(I64)n)<0? a+n: a){}
   constexpr U32 mod() const noexcept{
      return n;
   }
   template<typename INT> constexpr explicit operator INT() const noexcept{
      return n? (I64)a: (I64)t;
   }
   constexpr RingZnDyn operator+() const noexcept{
      return *this;
   }
   constexpr RingZnDyn operator-() const noexcept{
      RingZnDyn res = *this;
      if(n){
         res.a = a? n-a: 0;
      }else res.t = -t;
      return res;
   }
   constexpr RingZnDyn &operator+=(RingZnDyn rhs){
      if(!n){
         if(!rhs.n){
            t += rhs.t; return *this;
         }
         *this = RingZnDyn(rhs.n, t);
      }else if(!rhs.n){
         rhs = RingZnDyn(n, rhs.t);
      } // else assert(n == rhs.n);
      a += a<n-rhs.a? rhs.a: n-rhs.a;
      return *this;
   }
   constexpr RingZnDyn &operator-=(RingZnDyn rhs){
      if(!n){
         if(!rhs.n){
            t -= rhs.t; return *this;
         }
         *this = RingZnDyn(rhs.n, t);
      }else if(!rhs.n){
         rhs = RingZnDyn(n, rhs.t);
      } // else assert(n == rhs.n);
      a = a<rhs.a? a+(n-rhs.a): a-rhs.a;
      return *this;
   }
   constexpr RingZnDyn &operator*=(RingZnDyn rhs){
      if(!n){
         if(!rhs.n){
            t *= rhs.t; return *this;
         }
         *this = RingZnDyn(rhs.n, t);
      }else if(!rhs.n){
         rhs = RingZnDyn(n, rhs.t);
      } // else assert(n == rhs.n);
      a = (U64)a * rhs.a % n;
      return *this;
   }
   constexpr RingZnDyn &operator/=(RingZnDyn rhs){
      if(!n){
         if(!rhs.n){
            t /= rhs.t; return *this;
         }
         *this = RingZnDyn(rhs.n, t);
      }else if(!rhs.n){
         rhs = RingZnDyn(n, rhs.t);
      } // else assert(n == rhs.n);
      a = (U64)a * inv_mod(rhs.a, n) % n;
      return *this;
   }
   constexpr RingZnDyn &operator++() noexcept{
      if(n) a = a==n-1? 0: a+1;
      else ++t;
      return *this;
   }
   constexpr RingZnDyn operator++(int) noexcept{
      RingZnDyn res = *this; ++*this;
      return res;
   }
   constexpr RingZnDyn &operator--() noexcept{
      if(n) a = a==0? n-1: a-1;
      else --t;
      return *this;
   }
   constexpr RingZnDyn operator--(int) noexcept{
      RingZnDyn res = *this; --*this;
      return res;
   }
private:
   U32 n = 0;
   union{
      int t = 0;
      U32 a;
   };
};

#define DEF_BIOP(op)\
constexpr RingZnDyn operator op(RingZnDyn lhs, RingZnDyn rhs){\
   return lhs op##= rhs;\
}
DEF_BIOP(+)
DEF_BIOP(-)
DEF_BIOP(*)
DEF_BIOP(/)
#undef DEF_BIOP

constexpr bool operator==(RingZnDyn lhs, RingZnDyn rhs) noexcept{
   if(!lhs.mod()){
      if(!rhs.mod()){
         return (int)lhs == (int)rhs;
      }
      lhs = RingZnDyn(rhs.mod(), (int)lhs);
   }else if(!rhs.mod()){
      rhs = RingZnDyn(lhs.mod(), (int)rhs);
   }else if(lhs.mod() != rhs.mod()){
      return false;
   }
   return (U32)lhs == (U32)rhs;
}
constexpr bool operator!=(RingZnDyn lhs, RingZnDyn rhs) noexcept{
   return !(lhs == rhs);
}

#undef U64
#undef I64
#undef U32

#endif
