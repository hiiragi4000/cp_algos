#ifndef HOPCROFT_KARP_HH
#define HOPCROFT_KARP_HH

#include<algorithm>
#include<queue>
#include<utility>
#include<vector>
#include<climits>

struct HopcroftKarp{
   HopcroftKarp() = default;
   HopcroftKarp(int l, int r): g(l), x_mat(l, -1), y_mat(r, -1){}
   void insert_edge(int x, int y){
      g[x].push_back(y);
   }
   int maximum_matching(){
      while(bfs());
      return n_mat;
   }
   int left_match(int x) noexcept{
      return x_mat[x];
   }
   int right_match(int y) noexcept{
      return y_mat[y];
   }
private:
   std::vector<std::vector<int>> g;
   std::vector<int> x_mat, y_mat;
   int n_mat = 0;
   bool bfs(){
      std::vector<int> lv(g.size(), INT_MAX);
      std::queue<int> q;
      for(int x=0; x<(int)g.size(); ++x) if(x_mat[x] == -1){
         lv[x] = 0; q.push(x);
      }
      int d_min = INT_MAX;
      while(!q.empty()){
         int x = q.front(); q.pop();
         if(lv[x] >= d_min) continue;
         for(int y: g[x]) if(int x2 = y_mat[y]; x2 == -1){
            d_min = std::min(d_min, lv[x]+1);
         }else if(lv[x2] == INT_MAX){
            lv[x2] = lv[x]+1; q.push(x2);
         }
      }
      if(d_min == INT_MAX) return false;
      for(int x=0; x<(int)g.size(); ++x){
         n_mat += x_mat[x]==-1 && dfs(x, d_min, lv.data());
      }
      return true;
   }
   bool dfs(int x, int d_max, int *lv){
      if(lv[x] == d_max) return false;
      for(int y: g[x]) if(int x2=y_mat[y]; x2==-1 || lv[x2]==lv[x]+1 && dfs(x2, d_max, lv)){
         x_mat[x] = y; y_mat[y] = x;
         return true;
      }
      lv[x] = INT_MAX;
      return false;
   }
};

#endif
