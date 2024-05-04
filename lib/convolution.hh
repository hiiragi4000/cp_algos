#ifndef CONVOLUTION_HH
#define CONVOLUTION_HH

#include"number_theory.hh"
#include"util.hh"
#include<algorithm>
#include<complex>
#include<type_traits>
#include<utility>
#include<vector>
#include<cmath>

#define U32 unsigned
#define I64 long long
#define U64 unsigned long long

#define LOWBIT(NUM) ((NUM)&(~(NUM)+1))

template<U32 P, U32 R> struct BasicNtt{
   static constexpr U32 prime() noexcept{
      return P;
   }
   static constexpr U32 primitive_root() noexcept{
      return R;
   }
   static constexpr int lg_max_size() noexcept{
      if(primality_test(P) && primitive_root_modp_check(R, P)){
         U32 d = P-1, res = 0;
         if((d&65535) == 0) d>>=16, res|=16;
         if((d&255) == 0) d>>=8, res|=8;
         if((d&15) == 0) d>>=4, res|=4;
         if((d&3) == 0) d>>=2, res|=2;
         if((d&1) == 0) res|=1;
         return res;
      }
      return -1;
   }
   static constexpr size_t max_size() noexcept{
      int lg = lg_max_size();
      return lg>=0? 1ull<<lg: 0;
   }
   template<typename It>
   static void transform_in_place(bool inverse, It x, size_t n){
      // assert(n>=0 && n==LOWBIT(n) && n<=max_size());
      using T = DerefType<It>;
      for(size_t i=1, j=0; i<n; ++i){
         for(size_t k=n>>1; !((j^=k)&k); k>>=1);
         if(i < j) std::swap(x[i], x[j]);
      }
      RingZn<P> z2[lg_max_size()+1] = {pow_mod(R, inverse? (P-1)/n: (P-1)-(P-1)/n, P)};
      int lg = 0;
      for(; z2[lg]!=1; ++lg){
         z2[lg+1] = z2[lg]*z2[lg];
      }
      std::reverse(z2, z2+lg+1);
      for(size_t i=1; 1ull<<i<=n; ++i){
         for(size_t j=0; j<n; j+=1ull<<i){
            RingZn<P> t = 1;
            for(size_t k=0; k<1ull<<(i-1); ++k){
               RingZn<P> u = x[j|k], v = t*x[j|1ull<<(i-1)|k];
               x[j|k] = static_cast<T>(u+v);
               x[j|1ull<<(i-1)|k] = static_cast<T>(u-v);
               t *= z2[i];
            }
         }
      }
      if(inverse){
         RingZn<P> factor = static_cast<U32>(P-(P-1)/n);
         for(size_t i=0; i<n; ++i){
            x[i] = static_cast<T>(factor*x[i]);
         }
      }
   }
   template<typename It>
   static std::vector<RingZn<P>> transform(bool inverse, It b, It e, size_t n){
      std::vector<RingZn<P>> x(n);
      size_t d = e-b<n? e-b: n;
      copy(b, b+d, x.begin());
      transform_in_place(inverse, x.begin(), n);
      return x;
   }
};
template<U32 P> using Ntt = BasicNtt<P, primitive_root_modp(P)>;
using NttI32 = Ntt<(15<<27|1)>;
using NttU32 = Ntt<(3u<<30|1)>;

template<typename T> struct BasicFft{
   BasicFft() = default;
   template<typename It> void transform_in_place(bool inverse, It x, size_t n){
      // assert(n>=0 && n==LOWBIT(n));
      for(size_t i=1, j=0; i<n; ++i){
         for(size_t k=n>>1; !((j^=k)&k); k>>=1);
         if(i < j) std::swap(x[i], x[j]);
      }
      if(table.size() < n) calc_exp(n);
      for(size_t i=1; 1ull<<i<=n; ++i){
         size_t d = inverse? table.size()>>i: table.size()-(table.size()>>i);
         for(size_t j=0; j<n; j+=1ull<<i){
            for(size_t k=0, a=0; k<1ull<<(i-1); ++k, a=(a+d)&(table.size()-1)){
               std::complex<T> u = x[j|k], v = table[a]*x[j|1ull<<(i-1)|k];
               x[j|k] = u+v;
               x[j|1ull<<(i-1)|k] = u-v;
            }
         }
      }
      if(inverse){
         for(size_t i=0; i<n; ++i){
            x[i] /= static_cast<T>(n);
         }
      }
   }
   template<typename It>
   std::vector<std::complex<T>> transform(bool inverse, It b, It e, size_t n){
      std::vector<std::complex<T>> x(n);
      copy_n(b, e-b<n? e-b: n, x.begin());
      transform_in_place(inverse, x.begin(), n);
      return x;
   }
private:
   static constexpr T PI = 3.14159265358979323846l;
   void calc_exp(size_t n){
      if(n == 1){
         table = {1};
      }else if(n == 2){
         table = {1, -1};
      }else if(n == 4){
         table = {1, {0, 1}, -1, {0, -1}};
      }else if(n >= 8){
         table.resize(n);
         table[0] = 1;
         for(size_t i=1; i<n/8; ++i){
            table[i] = {std::cos(2*PI*i/n), std::sin(2*PI*i/n)};
         }
         T cos45 = std::sqrt((T).5);
         table[n/8] = {cos45, cos45};
         for(size_t i=n/8+1; i<=n/4; ++i){
            table[i] = {table[n/4-i].imag(), table[n/4-i].real()};
         }
         for(size_t i=n/4+1; i<=n/2; ++i){
            table[i] = {-table[n/2-i].real(), table[n/2-i].imag()};
         }
         for(size_t i=n/2+1; i<n; ++i){
            table[i] = -table[i-n/2];
         }
      }
   }
   std::vector<std::complex<T>> table;
};
using Fft = BasicFft<double>;
using Fftf = BasicFft<float>;
using Fftl = BasicFft<long double>;

namespace impl{
template<U32 N>
std::vector<RingZn<N>> convolution_pow2(std::vector<RingZn<N>> x, std::vector<RingZn<N>> y){
   size_t n = x.size();
   // assert(x.size()==y.size() && n==LOWBIT(n));
   if(primality_test(N) && n <= Ntt<N>::max_size()){
      if(x == y){
         Ntt<N>::transform_in_place(false, x.begin(), n);
         for(size_t i=0; i<n; ++i){
            x[i] *= x[i];
         }
      }else{
         Ntt<N>::transform_in_place(false, x.begin(), n);
         Ntt<N>::transform_in_place(false, y.begin(), n);
         for(size_t i=0; i<n; ++i){
            x[i] *= y[i];
         }
      }
      Ntt<N>::transform_in_place(true, x.begin(), n);
      return x;
   }
   std::vector<std::complex<double>> z(n);
   for(size_t i=0; i<n; ++i) z[i] = {static_cast<double>(x[i]), static_cast<double>(y[i])};
   Fft fft;
   fft.transform_in_place(false, z.begin(), n);
   for(size_t i=0; i<n; ++i){
      z[i] *= z[i];
   }
   fft.transform_in_place(true, z.begin(), n);
   if(static_cast<double>(n)*N*N < 2e15){
      for(size_t i=0; i<n; ++i){
         x[i] = static_cast<U64>(round(z[i].imag()/2));
      }
      return x;
   }
   std::vector<RingZn<NttU32::prime()>> a(n);
   for(size_t i=0; i<n; ++i){
      a[i] = static_cast<U32>(x[i]);
   }
   NttU32::transform_in_place(false, a.begin(), n);
   if(x == y){
      for(size_t i=0; i<n; ++i){
         a[i] *= a[i];
      }
   }else{
      std::vector<RingZn<NttU32::prime()>> b(n);
      for(size_t i=0; i<n; ++i){
         b[i] = static_cast<U32>(y[i]);
      }
      NttU32::transform_in_place(false, b.begin(), n);
      for(size_t i=0; i<n; ++i){
         a[i] *= b[i];
      }
   }
   NttU32::transform_in_place(true, a.begin(), n);
   for(size_t i=0; i<n; ++i){
      U64 q = round((z[i].imag()/2-static_cast<U32>(a[i]))/NttU32::prime());
      x[i] = (q%N)*NttU32::prime()+static_cast<U32>(a[i]);
   }
   return x;
}
} // namespace impl

template<
   typename InputIt1,
   typename InputIt2,
   typename OutputIt,
   typename = std::enable_if_t<
      RingZnMod<DerefType<InputIt1>>::value == RingZnMod<DerefType<InputIt2>>::value
      && RingZnMod<DerefType<InputIt2>>::value == RingZnMod<DerefType<OutputIt>>::value
      && RingZnMod<DerefType<OutputIt>>::value != 0
   >
> void convolution(InputIt1 b1, InputIt1 e1, InputIt2 b2, InputIt2 e2, size_t n, OutputIt res){
   constexpr U32 N = RingZnMod<DerefType<InputIt1>>::value;
   if(n <= 0) return;
   size_t n1 = std::min(static_cast<size_t>(e1-b1), n);
   size_t n2 = std::min(static_cast<size_t>(e2-b2), n);
   if(n1<=32 || n2<=32){
      std::vector<RingZn<N>> t(n);
      for(size_t i=0; i<n1; ++i) for(size_t j=0; j<n2; ++j){
         t[(i+j)%n] += b1[i]*b2[j];
      }
      copy(t.cbegin(), t.cend(), res);
      return;
   }
   size_t nf = n;
   while(nf != LOWBIT(nf)) nf += LOWBIT(nf);
   if(nf>n && n<n1+n2-1) nf *= 2;
   std::vector<RingZn<N>> x(nf), y(nf);
   copy(b1, b1+n1, x.begin());
   copy(b2, b2+n2, y.begin());
   if(nf > 2*n){
      copy(b1+1, b1+n1, x.end()-n+1);
   }
   x = impl::convolution_pow2(move(x), move(y));
   copy(x.cbegin(), x.cbegin()+n, res);
}

template<typename InputIt1, typename InputIt2, typename OutputIt>
void convolution_int(InputIt1 b1, InputIt1 e1, InputIt2 b2, InputIt2 e2, size_t n, OutputIt res){
   if(n <= 0) return;
   using OutT = DerefType<OutputIt>;
   size_t n1 = e1-b1, n2 = e2-b2;
   n1 = std::min(n1, n); n2 = std::min(n2, n);
   if(n1<=32 || n2<=32){
      std::vector<OutT> t(n);
      for(size_t i=0; i<n1; ++i) for(size_t j=0; j<n2; ++j){
         t[(i+j)%n] += static_cast<OutT>(b1[i])*static_cast<OutT>(b2[j]);
      }
      std::copy(t.cbegin(), t.cend(), res);
      return;
   }
   size_t n_fft = n;
   while((n_fft & -n_fft) != n_fft){
      n_fft += n_fft & -n_fft;
   }
   if(n_fft>n && n<n1+n2-1){
      n_fft <<= 1;
   }
   std::vector<std::complex<double>> x(n_fft);
   for(size_t i=0; i<n1; ++i){
      x[i] = b1[i];
   }
   for(size_t i=0; i<n2; ++i){
      x[i] += std::complex<double>(0, b2[i]);
   }
   if(n_fft > 2*n){
      for(size_t i=1; i<n1; ++i){
         x[n_fft-n+i] = b1[i];
      }
   }
   Fft fft;
   fft.transform_in_place(false, x.begin(), n_fft);
   for(size_t i=0; i<n_fft; ++i){
      x[i] *= x[i];
   }
   fft.transform_in_place(true, x.begin(), n_fft);
   for(size_t i=0; i<n; ++i){
      res[i] = static_cast<OutT>(round(x[i].imag()/2));
   }
}

#undef LOWBIT

#undef U64
#undef I64
#undef U32

#endif // CONVOLUTION_HH
