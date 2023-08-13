#include"disjoint_sets.hh"
#include<cstdio>
using namespace std;

int main(){
   int n, q;
   scanf("%d%d", &n, &q);
   DisjointSets djs(n);
   while(q-- > 0){
      int t, u, v;
      scanf("%d%d%d", &t, &u, &v);
      if(t == 0){
         djs.merge(u, v);
      }else{
         printf("%d\n", djs.find(u)==djs.find(v));
      }
   }
}
