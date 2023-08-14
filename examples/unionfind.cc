// Verified in the problem https://judge.yosupo.jp/problem/unionfind
#include"disjoint_sets.hh"
#include<cstdio>
using namespace std;

int main(){
   int n, q;
   scanf("%d%d", &n, &q);
   // Set the disjoint sets as {0}, {1}, ..., {n-1}.
   DisjointSets djs(n);
   while(q-- > 0){
      int t, u, v;
      scanf("%d%d%d", &t, &u, &v);
      if(t == 0){
         // Suppose u \in U and v \in V.
         // If U = V, the method returns false.
         // Otherwise, the method returns true and merge U and V.
         djs.merge(u, v);
      }else{
         // Suppose u \in U.
         // |djs.find(u)| returns the representative of U.
         printf("%d\n", djs.find(u)==djs.find(v));
      }
   }
}
