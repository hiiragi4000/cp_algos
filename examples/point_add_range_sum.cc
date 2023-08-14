// Verified in the problem https://judge.yosupo.jp/problem/point_add_range_sum
#include"range_query.hh"
#include<cstdio>
using namespace std;
using i64 = long long;

int main(){
   int n, q;
   scanf("%d%d", &n, &q);
   // Set a Fenwick tree (aka binary indexed tree) from a[0..n-1], where a[0] = a[1] = ... = a[n-1] = 0.
   Fenwick<i64> bit(n);
   for(int i=0; i<=n-1; ++i){
      int ai;
      scanf("%d", &ai);
      // Fenwick::add(pos, x): add x to a[pos].
      bit.add(i, ai);
   }
   while(q-- > 0){
      int cmd[3];
      scanf("%d%d%d", cmd, cmd+1, cmd+2);
      if(cmd[0] == 0){
         bit.add(cmd[1], cmd[2]);
      }else{
         // Fenwick::sum(i): returns a[0]+a[1]+...+a[i].
         // If i = -1, returns 0.
         printf("%lld\n", bit.sum(cmd[2]-1)-bit.sum(cmd[1]-1));
      }
   }
}
