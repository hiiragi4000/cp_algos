#ifndef CONVOLUTION_HH
#define CONVOLUTION_HH

#include"number_theory.hh"
#include<algorithm>
#include<complex>

#define I64 long long

template<int P, int R=primitive_root_modp(P)> struct BasicNtt{
   static_assert(primality_test(P));
   static_assert(primitive_root_modp_check(R, P));
   static constexpr int prime() noexcept{
      return P;
   }
   static constexpr int primitive_root() noexcept{
      return R;
   }
   static constexpr int lg_max_size() noexcept{
      int d = P-1, res = 0;
      if((d&65535) == 0) d>>=16, res|=16;
      if((d&255) == 0) d>>=8, res|=8;
      if((d&15) == 0) d>>=4, res|=4;
      if((d&3) == 0) d>>=2, res|=2;
      if((d&1) == 0) res|=1;
      return res;
   }
   static constexpr int max_size() noexcept{
      return 1<<lg_max_size();
   }
   template<typename It>
   static void transform_in_place(bool inverse, It x, int n){
      // assert(n>=0 && (n&-n)==n && n<=max_size());
      // assert(*std::min_element(x, x+n) >= 0);
      for(int i=1, j=0; i<n; ++i){
         for(int k=n>>1; !((j^=k)&k); k>>=1);
         if(i < j) std::swap(x[i], x[j]);
      }
      int z2[lg_max_size()+1] = {(int)pow_mod(R, inverse? (P-1)/n: -(P-1)/n, P)}, lg = 1;
      for(; z2[lg-1]!=1; ++lg){
         z2[lg] = (I64)z2[lg-1]*z2[lg-1]%P;
      }
      std::reverse(z2, z2+lg);
      for(int i=1; 1<<i<=n; ++i){
         for(int j=0; j<n; j+=1<<i){
            for(int k=0, t=1; k<1<<(i-1); ++k, t=(I64)t*z2[i]%P){
               I64 u = x[j|k], v = (I64)t*x[j|1<<(i-1)|k]%P;
               x[j|k] = (u+v)%P;
               x[j|1<<(i-1)|k] = (u+P-v)%P;
            }
         }
      }
      if(inverse){
         for(int i=0; i<n; ++i){
            x[i] = x[i]*(P-(P-1ll)/n)%P;
         }
      }
   }
   template<typename It>
   static std::vector<int> transform(bool inverse, It b, It e, int n){
      std::vector<int> x(n);
      int d = e-b<n? e-b: n;
      for(int i=0; i<d; ++i){
         x[i] = (b[i]%(I64)P+P)%P;
      }
      transform_in_place(inverse, x.begin(), n);
      return x;
   }
};
using Ntt = BasicNtt<(15<<27|1)>;

template<typename T> struct BasicFft{
   template<typename It>
   static void transform_in_place(bool inverse, It x, int n){
      // assert(n>=0 && (n&-n)==n);
      for(int i=1, j=0; i<n; ++i){
         for(int k=n>>1; !((j^=k)&k); k>>=1);
         if(i < j) std::swap(x[i], x[j]);
      }
      auto exp_table = calc_exp(n);
      for(int i=1; 1<<i<=n; ++i){
         int d = inverse? n>>i: n-(n>>i);
         for(int j=0; j<n; j+=1<<i){
            for(int k=0, a=0; k<1<<(i-1); ++k, a=a+d&n-1){
               std::complex<T> u = x[j|k], v = exp_table[a]*x[j|1<<(i-1)|k];
               x[j|k] = u+v;
               x[j|1<<(i-1)|k] = u-v;
            }
         }
      }
      if(inverse){
         for(int i=0; i<n; ++i){
            x[i] /= n;
         }
      }
   }
   template<typename It>
   static std::vector<std::complex<T>> transform(bool inverse, It b, It e, int n){
      std::vector<std::complex<T>> x(n);
      copy_n(b, e-b<n? e-b: n, x.begin());
      transform_in_place(inverse, x.begin(), n);
      return x;
   }
private:
   static constexpr T PI = 3.14159265358979323846l;
   static std::vector<std::complex<T>> calc_exp(int n){
      std::vector<std::complex<T>> table;
      if(n == 1){
         table = {1};
      }else if(n == 2){
         table = {1, -1};
      }else if(n == 4){
         table = {1, {0, 1}, -1, {0, -1}};
      }else if(n >= 8){
         table.resize(n);
         table[0] = {1, 0};
         for(int i=1; i<n/8; ++i){
            table[i] = {std::cos(2*PI*i/n), std::sin(2*PI*i/n)};
         }
         T cos45 = std::sqrt((T).5);
         table[n/8] = {cos45, cos45};
         for(int i=n/8+1; i<=n/4; ++i){
            table[i] = {table[n/4-i].imag(), table[n/4-i].real()};
         }
         for(int i=n/4+1; i<=n/2; ++i){
            table[i] = {-table[n/2-i].real(), table[n/2-i].imag()};
         }
         for(int i=n/2+1; i<n; ++i){
            table[i] = -table[i-n/2];
         }
      }
      return table;
   }
};
using Fft = BasicFft<double>;
using Fftf = BasicFft<float>;
using Fftl = BasicFft<long double>;

#undef I64

#endif
