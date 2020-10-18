#ifndef LCA_HH
#define LCA_HH

#include"rmq.hh"
#include<utility>

struct Lca{
   Lca() = default;
   Lca(std::vector<int> const *g, int n): fi(n), dep(2*n), vis(2*n){
      int n_step = 0;
      dfs(g, 0, 0, 0, n_step);
      rmq = Rmq(dep.cbegin(), dep.size());
   }
   int query(int u, int v) const{
      auto l = fi[u], r = fi[v];
      if(l > r) std::swap(l, r);
      auto i4 = rmq.query(l, r);
      int i = i4[0];
      if(dep[i4[1]] < dep[i]) i = i4[1];
      if(dep[i4[2]] < dep[i]) i = i4[2];
      if(dep[i4[3]] < dep[i]) i = i4[3];
      return vis[i];
   }
private:
   std::vector<int> fi, dep, vis;
   Rmq rmq;
   void dfs(std::vector<int> const *g, int u, int p, int d, int &n_step){
      fi[u] = n_step;
      dep[n_step] = d;
      vis[n_step++] = u;
      for(int v: g[u]) if(v != p){
         dfs(g, v, u, d+1, n_step);
         dep[n_step] = d;
         vis[n_step++] = u;
      }
   }
};

#endif
