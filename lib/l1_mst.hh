#ifndef L1_MST_HH
#define L1_MST_HH

#include"disjoint_sets.hh"
#include"point2.hh"
#include<algorithm>
#include<map>
#include<numeric>
#include<tuple>
#include<utility>
#include<vector>

template<typename T>
std::vector<std::pair<int, int>> l1_mst(Point2<T> const *pts, int n){
   static constexpr std::pair<Mat2<int>, Point2<int>> CONES[]{
      {{1, 0, -1, 1}, {2, 1}},
      {{1, -1, 0, 1}, {1, 2}},
      {{0, -1, 1, 1}, {2, 1}},
      {{-1, -1, 1, 0}, {1, 2}}
   };
   std::vector<Point2<T>> tp(n);
   std::vector<int> si(n);
   iota(si.begin(), si.end(), 0);
   std::vector<std::tuple<T, int, int>> g;
   g.reserve(4*n);
   for(auto [A, coef]: CONES){
      for(int i=0; i<n; ++i){
         tp[i] = A*pts[i];
      }
      sort(si.begin(), si.end(), [&tp](int i, int j){
         return tp[j] < tp[i];
      });
      std::map<T, int> sl;
      for(int i0=0; i0<n; ++i0){
         int i = si[i0];
         if(i0>0 && tp[si[i0-1]]==tp[i]){
            g.emplace_back(0, i, si[i0-1]);
            continue;
         }
         if(auto it=sl.upper_bound(tp[i].y); it!=sl.end()){
            g.emplace_back(coef*tp[it->second]-coef*tp[i], i, it->second);
         }
         auto it = sl.insert_or_assign(tp[i].y, i).first;
         while(it!=sl.begin() && coef*tp[i]<=coef*tp[prev(it)->second]){
            sl.erase(prev(it));
         }
      }
   }
   DisjointSets djs(n);
   std::vector<std::pair<int, int>> res;
   res.reserve(n-1);
   sort(g.begin(), g.end());
   for(auto [w, u, v]: g) if(djs.find(u) != djs.find(v)){
      djs.merge(u, v);
      res.emplace_back(u, v);
   }
   return res;
}

#endif // L1_MST_HH
