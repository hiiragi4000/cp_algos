// Verified in the problem https://judge.yosupo.jp/problem/multiplication_of_big_integers
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
      puts(string(a*b).data());
   }
}
