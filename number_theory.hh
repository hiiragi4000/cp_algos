#ifndef NUMBER_THEORY_HH
#define NUMBER_THEORY_HH

#include<array>
#include<functional>
#include<limits>
#include<numeric>
#include<ostream>
#include<type_traits>
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
   mu[0] = 0; mu[1] = 1;
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
   // assert(mod>1 && std::gcd(mod, a%mod)==1);
   if((a%=mod) < 0) a += mod;
   auto [_, x] = extgcd(mod, a);
   return x<0? x+mod: x;
}

constexpr U32 pow_mod(I64 a, I64 n, U32 mod){
   // assert(mod > 0);
   if(mod == 1) return 0;
   if((a%=mod) < 0) a += mod;
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
   return miller_rabin(n, 2) && (n < 2047 || (miller_rabin(n, 7) && miller_rabin(n, 61)));
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

constexpr I64 sqrt_modp(I64 a, U32 p){
   // assert(primality_test(p));
   if((a%=p) < 0) a += p;
   if(a <= 1) return a;
   if(jacobi_symbol(a, p) == -1) return -1;
   if(p%4 == 3) return pow_mod(a, (p+1ll)/4, p);
   I64 k = 1;
   while(jacobi_symbol(k*k-a, p) == 1) ++k;
   if(k*k%p == a) return k;
   U64 w = ((k*k-a)%p+p)%p;
   U64 rc = 1, rx = 0, tc = k, tx = 1;
   for(U64 i=(p+1ll)/2; i>0; ){
      if(i%2){
         U64 rc2 = (rc*tc%p + rx*tx%p*w)%p;
         rx = (rc*tx%p + rx*tc)%p;
         rc = rc2;
      }
      if(i/=2){
         U64 tc2 = (tc*tc%p + tx*tx%p*w)%p;
         tx = (tc<p-tc? 2*tc: tc-(p-tc))*tx%p;
         tc = tc2;
      }
   }
   if(rc > p-rc) rc = p-rc;
   return rc;
}

constexpr bool primitive_root_modp_check(U32 r, U32 p){
   // assert(primality_test(p));
   if(r%p == 0) return false;
   U32 d = p-1;
   for(I64 q=2; q*q<=d; ++q) if(d%q == 0){
      do d/=q; while(d%q==0);
      if(pow_mod(r, (p-1)/q, p) == 1) return false;
   }
   if(d>1 && pow_mod(r, (p-1)/d, p) == 1){
      return false;
   }
   return true;
}

constexpr U32 primitive_root_modp(U32 p){
   // assert(primality_test(p));
   for(U32 i=1; ; ++i){
      if(primitive_root_modp_check(i, p)){
         return i;
      }
   }
}

inline I64 discrete_log(I64 base, I64 power, U32 mod){
   // assert(mod > 0);
   if((base%=mod) < 0) base += mod;
   if((power%=mod) < 0) power += mod;
   if(mod==1 || power==1) return 0;
   U32 s = 1, d = 1, res = 0;
   while(1){
      s = (U64)base*s%mod;
      ++res;
      if(s == power) return res;
      U32 d2 = std::gcd(s, mod);
      if(d2 == d){
         if(power%d) return -1;
         s /= d; power /= d; mod /= d;
         break;
      }
      d = d2;
   }
   U32 m = std::ceil(std::sqrt(mod));
   std::unordered_map<U32, U32> bs(m);
   for(U32 i=0; i<m; ++i){
      bs.emplace(s, i);
      s = (U64)base*s%mod;
   }
   U32 gs = pow_mod(base, -(I64)m, mod);
   for(U32 i=0; i<m; ++i){
      if(auto it=bs.find(power); it!=bs.end()){
         return res + i*m + it->second;
      }
      power = (U64)gs*power%mod;
   }
   return -1;
}

template<U32 N> struct RingZn{
   using Base = std::conditional_t<(N<=(U32)std::numeric_limits<int>::max()), int, U32>;
   static constexpr U32 mod() noexcept{
      return N;
   }
   RingZn() = default;
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>
   constexpr RingZn(INT a) noexcept: a((a%=(I64)N)<0? (I64)a+N: a){}
   template<typename INT, typename = std::enable_if_t<std::is_arithmetic_v<INT>>>
   constexpr explicit operator INT() const noexcept{
      return a;
   }
   constexpr RingZn operator+() const noexcept{
      return *this;
   }
   constexpr RingZn operator-() const noexcept{
      return a? (Base)N-a: 0;
   }
   constexpr RingZn &operator+=(RingZn rhs) noexcept{
      a = a<(Base)N-rhs.a? a+rhs.a: a-((Base)N-rhs.a);
      return *this;
   }
   constexpr RingZn &operator-=(RingZn rhs) noexcept{
      a = a<rhs.a? a+((Base)N-rhs.a): a-rhs.a;
      return *this;
   }
   constexpr RingZn &operator*=(RingZn rhs) noexcept{
      a = (U64)a * rhs.a % N;
      return *this;
   }
   constexpr RingZn &operator/=(RingZn rhs) noexcept{
      a = (U64)a * inv_mod(rhs.a, N) % N;
      return *this;
   }
   constexpr RingZn &operator++() noexcept{
      a = a==(Base)N-1? 0: a+1;
      return *this;
   }
   constexpr RingZn operator++(int) noexcept{
      RingZn res = *this; ++*this;
      return res;
   }
   constexpr RingZn &operator--() noexcept{
      a = a? a-1: (Base)N-1;
      return *this;
   }
   constexpr RingZn operator--(int) noexcept{
      RingZn res = *this; --*this;
      return res;
   }
#define DEF_BIOP(op)\
   friend constexpr RingZn operator op(RingZn lhs, RingZn rhs) noexcept{\
      return lhs op##= rhs;\
   }
   DEF_BIOP(+)
   DEF_BIOP(-)
   DEF_BIOP(*)
   DEF_BIOP(/)
#undef DEF_BIOP
   friend constexpr bool operator==(RingZn lhs, RingZn rhs) noexcept{
      return (U32)lhs == (U32)rhs;
   }
   friend constexpr bool operator!=(RingZn lhs, RingZn rhs) noexcept{
      return !(lhs == rhs);
   }
private:
   Base a = 0;
};
template<U32 N> RingZn<N> pow(RingZn<N> a, I64 b){
   return pow_mod((I64)a, b, N);
}
template<U32 N> std::ostream &operator<<(std::ostream &os, RingZn<N> a){
   return os << (U32)a;
}

template<typename> struct RingZnMod: std::integral_constant<U32, 0>{};
template<U32 N> struct RingZnMod<RingZn<N>>: std::integral_constant<U32, N>{};

#undef U64
#undef I64
#undef U32

#endif
