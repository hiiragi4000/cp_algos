#ifndef DISJOINT_SETS_HH
#define DISJOINT_SETS_HH

#include"cow.hh"

struct DisjointSets{
   DisjointSets() = default;
   explicit DisjointSets(int n): p(n, -1), s(n, 1){}
   int find(int i){
      int r = i;
      while(p[r]!=-1) r = p[r];
      for(int j; i!=r; ) j = p[i], p[i] = r, i = j;
      return r;
   }
   bool merge(int i, int j){
      i = find(i); j = find(j);
      if(i == j) return false;
      if(s[i] < s[j]) std::swap(i, j);
      p[j] = i; s[i] += s[j];
      return true;
   }
   int component_size(int i){
      return s[find(i)];
   }
private:
   std::vector<int> p, s;
};

struct DjsCow: private VlaCow<std::pair<int, int>>{
   using Base = VlaCow<std::pair<int, int>>;
   using Base::size;
   using Base::load;
   using Base::save;
   DjsCow() = default;
   explicit DjsCow(int n): Base(n, {-1, 0}){}
   int find(int u) const noexcept{
      int r = u;
      for(int p; (p=std::get<PARENT>((*this)[r]))!=-1; ){
         r = p;
      }
      return r;
   }
   bool merge(int u, int v){
      u = find(u); v = find(v);
      if(u == v) return false;
      int ru = std::get<RANK>((*this)[u]), rv = std::get<RANK>((*this)[v]);
      if(ru < rv){
         std::swap(u, v);
         std::swap(ru, rv);
      }
      write(v, {u, rv});
      if(ru == rv){
         write(u, {-1, ru+1});
      }
      return true;
   }
private:
   enum {PARENT, RANK};
};

#endif
