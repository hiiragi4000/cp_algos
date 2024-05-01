// Verified in the problem https://judge.yosupo.jp/problem/division_of_big_integers
#include"arbitrary_prec.hh"
#include<iostream>
#include<string>
#include<cstdio>
using namespace std;

int main(){
   int T;
   cin >> T;
   while(T-- > 0){
      string sa, sb;
      cin >> sa >> sb;
      BigInt a = stob(sa.data()), b = stob(sb.data());
      // DivT<BigInt> div(BigInt const &a, BigInt const &b):
      // Computes both the quotient `quot` and the remainder `rem` of the division of the numerator `a` by the denominator `b`.
      // If b = 0, raises SIGFPE.
      // The quotient is the algebraic quotient with any fractional part discarded (truncated towards zero).
      // The remainder is such that quot*y + rem = x.
      // The struct template DivT is defined as follows:
      // template<typename T> struct DivT<T> { T quot, rem; };
      auto [q, r] = div(a, b);
      printf("%s %s\n", string(q).data(), string(r).data());
   }
}
