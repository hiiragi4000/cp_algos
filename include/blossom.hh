#ifndef BLOSSOM_HH
#define BLOSSOM_HH

#include<algorithm>
#include<queue>
#include<utility>
#include<vector>

struct Blossom{
   Blossom() = default;
   explicit Blossom(int n): g(n), mat(n, -1){}
   explicit Blossom(std::vector<std::vector<int>> g0): g(move(g0)), mat(g.size(), -1){}
   void insert_edge(int u, int v){
      g[u].push_back(v);
      g[v].push_back(u);
   }
   bool find_augmenting_path(int s){
      if(mat[s] >= 0) return false;
      std::vector<int> b(g.size(), -1), p(g.size()), ce1(g.size()), ce2(g.size()), tag(g.size());
      std::queue<int> bfs;
      std::vector<signed char> colo(g.size(), -1);
      int n_bl = 0;
      auto basis = [&](int i){
         int r = i;
         while(b[r]!=-1) r = b[r];
         for(int j; i!=r; ) j = b[i], b[i] = r, i = j;
         return r;
      };
      auto lca = [&](int u, int v){
         u = basis(u); v = basis(v);
         ++n_bl;
         while(tag[u] != n_bl){
            tag[u] = n_bl;
            if(u != s) u = basis(p[u]);
            std::swap(u, v);
         }
         return u;
      };
      auto contract = [&](int u, int v, int a){
         for(int i=basis(u); i!=a; i=basis(p[i])){
            b[i] = a;
            if(colo[i] == 1){
               ce1[i] = u; ce2[i] = v;
               bfs.push(i);
            }
         }
      };
      bfs.push(s); colo[s] = 0;
      while(!bfs.empty()){
         int u = bfs.front(); bfs.pop();
         for(int v: g[u]) if(mat[v]!=u && basis(v)!=basis(u)){
            if(colo[v] == -1){
               if(mat[v] == -1){
                  alternate(s, u, p.data(), ce1.data(), ce2.data(), colo.data());
                  mat[u] = v; mat[v] = u;
                  ++n_mat;
                  return true;
               }
               p[v] = u; colo[v] = 1;
               p[mat[v]] = v; colo[mat[v]] = 0;
               bfs.push(mat[v]);
            }else if(colo[basis(v)] == 0){
               int a = lca(u, v);
               contract(u, v, a);
               contract(v, u, a);
            }
         }
      }
      return false;
   }
   int maximum_matching(){
      for(int u=0; u<(int)g.size(); ++u) if(mat[u] == -1){
         int mu = -1;
         for(int v: g[u]) if(mat[v] == -1){
            mu = v; break;
         }
         if(mu >= 0){
            mat[u] = mu; mat[mu] = u;
            ++n_mat;
         }
      }
      for(int u=0; u<(int)g.size(); ++u){
         find_augmenting_path(u);
      }
      return n_mat;
   }
   int match(int u) noexcept{
      return mat[u];
   }
private:
   std::vector<std::vector<int>> g;
   std::vector<int> mat;
   int n_mat = 0;
   void alternate(int s, int t, int const *p, int const *ce1, int const *ce2, signed char *colo){
      if(s == t) return;
      if(colo[t] == 0){
         int pt = p[t], ppt = p[pt];
         alternate(s, ppt, p, ce1, ce2, colo);
         mat[pt] = ppt; mat[ppt] = pt;
      }else{
         int c1 = ce1[t], c2 = ce2[t];
         alternate(mat[t], c1, p, ce1, ce2, colo);
         alternate(s, c2, p, ce1, ce2, colo);
         mat[c1] = c2; mat[c2] = c1;
      }
   }
};

#endif // BLOSSOM_HH
