// Verified in the problem https://judge.yosupo.jp/problem/many_aplusb_128bit
#include"int128.hh"
#include<cstdio>
using namespace std;

int main(){
   int T;
   scanf("%d", &T);
   while(T-- > 0){
      char sa[48], sb[48];
      scanf("%47s%47s", sa, sb);
      i128 a = stoi128(sa), b = stoi128(sb);
      puts(to_string(a+b).data());
   }
}
