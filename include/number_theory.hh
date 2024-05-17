#ifndef NUMBER_THEORY_HH
#define NUMBER_THEORY_HH

#include"assert_platform.hh"
#include<array>
#include<limits>
#include<numeric>
#include<ostream>
#include<type_traits>
#include<unordered_map>
#include<utility>
#include<vector>
#include<climits>
#include<cmath>

#define I8 signed char
#define U8 unsigned char
#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

struct IsprimePrimePair{
   std::vector<U8> isprime;
   std::vector<int> prime;
};

inline IsprimePrimePair isprime_prime_pair(int n){
   if(n <= 0){
      return {};
   }
   if(n == 1){
      return {{false}, {}};
   }
   std::vector<U8> inp(n, true);
   inp[0] = inp[1] = false;
   std::vector<int> p;
   for(int i=2; i<=n-1; ++i){
      if(inp[i]){
         p.push_back(i);
      }
      for(int pj: p){
         I64 m = static_cast<I64>(i)*pj;
         if(m >= n){
            break;
         }
         inp[m] = false;
         if(i%pj == 0){
            break;
         }
      }
   }
   return {std::move(inp), std::move(p)};
}

inline std::vector<U8> isprime_table(int n){
   return isprime_prime_pair(n).isprime;
}

inline std::vector<int> prime_list(int n){
   return isprime_prime_pair(n).prime;
}

struct PhiPrimePair{
   std::vector<int> phi, prime;
};

inline PhiPrimePair phi_prime_pair(int n){
   if(n <= 0){
      return {};
   }
   std::vector<int> phi(n);
   std::iota(phi.begin(), phi.end(), 0);
   std::vector<int> p;
   for(int i=2; i<=n-1; ++i){
      if(phi[i] == i){
         --phi[i];
         p.push_back(i);
      }
      for(int pj: p){
         I64 m = static_cast<I64>(i)*pj;
         if(m >= n){
            break;
         }
         if(i%pj == 0){
            phi[m] = pj*phi[i];
            break;
         }
         phi[m] = (pj-1)*phi[i];
      }
   }
   return {std::move(phi), std::move(p)};
}

inline std::vector<int> phi_table(int n){
   return phi_prime_pair(n).phi;
}

constexpr U64 euler_phi(U64 n) noexcept{
   U64 res = n;
   for(U64 i=2; i<=UINT_MAX && i*i<=n; ++i){
      if(n%i == 0){
         res = res/i*(i-1);
         do{
            n /= i;
         }while(n%i == 0);
      }
   }
   if(n > 1){
      res = res/n*(n-1);
   }
   return res;
}

struct MuPrimePair{
   std::vector<I8> mu;
   std::vector<int> prime;
};

inline MuPrimePair mu_prime_pair(int n){
   if(n <= 0){
      return {};
   }
   if(n == 1){
      return {{0}, {}};
   }
   std::vector<I8> mu(n, -2);
   mu[1] = 1;
   std::vector<int> p;
   for(int i=2; i<=n-1; ++i){
      if(mu[i] == -2){
         mu[i] = -1;
         p.push_back(i);
      }
      for(int pj: p){
         I64 m = static_cast<I64>(i)*pj;
         if(m >= n){
            break;
         }
         if(i%pj == 0){
            mu[m] = 0;
            break;
         }
         mu[m] = -mu[i];
      }
   }
   return {std::move(mu), std::move(p)};
}

inline std::vector<I8> mu_table(int n){
   return mu_prime_pair(n).mu;
}

constexpr int moebius_mu(U64 n) noexcept{
   int res = 1;
   for(U64 i=2; i<=UINT_MAX && i*i<=n; ++i){
      if(n%i == 0){
         n /= i;
         if(n%i == 0){
            return 0;
         }
         res = -res;
      }
   }
   if(n > 1){
      res = -res;
   }
   return res;
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
   return static_cast<U32>(res);
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
      x = static_cast<U64>(x)*x%n;
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
   return n<0? false: primality_test(static_cast<U32>(n));
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
      I64 t = b%a; b = a; a = t;
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
      do d/=static_cast<U32>(q); while(d%q==0);
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
      s = static_cast<U64>(base)*s%mod;
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
   U32 m = static_cast<U32>(std::ceil(std::sqrt(mod)));
   std::unordered_map<U32, U32> bs(m);
   for(U32 i=0; i<m; ++i){
      bs.emplace(s, i);
      s = static_cast<U64>(base)*s%mod;
   }
   U32 gs = pow_mod(base, -static_cast<I64>(m), mod);
   for(U32 i=0; i<m; ++i){
      if(auto it=bs.find(static_cast<U32>(power)); it!=bs.end()){
         return res + i*m + it->second;
      }
      power = static_cast<U64>(gs)*power%mod;
   }
   return -1;
}

template<U32 N> struct RingZn{
   using Base = std::conditional_t<(N<=static_cast<U32>(std::numeric_limits<int>::max())), int, U32>;
   static constexpr U32 mod() noexcept{
      return N;
   }
   RingZn() = default;
   template<typename INT, typename = std::enable_if_t<std::is_integral_v<INT>>>
   constexpr RingZn(INT a) noexcept: a((a%=static_cast<I64>(N))<0? static_cast<I64>(a)+N: a){}
   template<typename INT, typename = std::enable_if_t<std::is_arithmetic_v<INT>>>
   constexpr explicit operator INT() const noexcept{
      return a;
   }
   constexpr RingZn operator+() const noexcept{
      return *this;
   }
   constexpr RingZn operator-() const noexcept{
      return a? static_cast<Base>(N)-a: 0;
   }
   constexpr RingZn &operator+=(RingZn rhs) noexcept{
      a = a<static_cast<Base>(N)-rhs.a? a+rhs.a: a-(static_cast<Base>(N)-rhs.a);
      return *this;
   }
   constexpr RingZn &operator-=(RingZn rhs) noexcept{
      a = a<rhs.a? a+(static_cast<Base>(N)-rhs.a): a-rhs.a;
      return *this;
   }
   constexpr RingZn &operator*=(RingZn rhs) noexcept{
      a = static_cast<U64>(a) * rhs.a % N;
      return *this;
   }
   constexpr RingZn &operator/=(RingZn rhs) noexcept{
      a = static_cast<U64>(a) * inv_mod(rhs.a, N) % N;
      return *this;
   }
   constexpr RingZn &operator++() noexcept{
      a = a==static_cast<Base>(N)-1? 0: a+1;
      return *this;
   }
   constexpr RingZn operator++(int) noexcept{
      RingZn res = *this; ++*this;
      return res;
   }
   constexpr RingZn &operator--() noexcept{
      a = a? a-1: static_cast<Base>(N)-1;
      return *this;
   }
   constexpr RingZn operator--(int) noexcept{
      RingZn res = *this; --*this;
      return res;
   }
#define DEF_BIOP(OP)\
   friend constexpr RingZn operator OP(RingZn lhs, RingZn rhs) noexcept{\
      return lhs OP##= rhs;\
   }
   DEF_BIOP(+)
   DEF_BIOP(-)
   DEF_BIOP(*)
   DEF_BIOP(/)
#undef DEF_BIOP
   friend constexpr bool operator==(RingZn lhs, RingZn rhs) noexcept{
      return static_cast<U32>(lhs) == static_cast<U32>(rhs);
   }
   friend constexpr bool operator!=(RingZn lhs, RingZn rhs) noexcept{
      return !(lhs == rhs);
   }
private:
   Base a = 0;
};
template<U32 N> constexpr RingZn<N> pow(RingZn<N> a, I64 b){
   return pow_mod(static_cast<I64>(a), b, N);
}
template<U32 N> std::ostream &operator<<(std::ostream &os, RingZn<N> a){
   return os << static_cast<U32>(a);
}

template<typename> struct RingZnMod: std::integral_constant<U32, 0>{};
template<U32 N> struct RingZnMod<RingZn<N>>: std::integral_constant<U32, N>{};

#undef U64
#undef I64
#undef U32
#undef U8
#undef I8

#endif // NUMBER_THEORY_HH
