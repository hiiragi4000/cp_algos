#ifndef DISJOINT_SETS_HH
#define DISJOINT_SETS_HH

#include"cow.hh"

struct DisjointSets{
   DisjointSets() = default;
   explicit DisjointSets(int n): p(n, -1), s(n, 1){}
   int find(int i) noexcept{
      int r = static_cast<DisjointSets const*>(this)->find(i);
      for(int j; i!=r; ) j = p[i], p[i] = r, i = j;
      return r;
   }
   int find(int i) const noexcept{
      int r = i;
      while(p[r] != -1) r = p[r];
      return r;
   }
   bool merge(int i, int j) noexcept{
      i = find(i); j = find(j);
      if(i == j) return false;
      if(s[i] < s[j]) std::swap(i, j);
      p[j] = i; s[i] += s[j];
      return true;
   }
   int component_size(int i) noexcept{
      return s[find(i)];
   }
   int component_size(int i) const noexcept{
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
   explicit DjsCow(int n): Base(n, {-1, 1}){}
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
      int su = std::get<SIZE>((*this)[u]), sv = std::get<SIZE>((*this)[v]);
      if(su < sv){
         std::swap(u, v);
         std::swap(su, sv);
      }
      write(v, {u, sv});
      write(u, {-1, su+sv});
      return true;
   }
   int component_size(int u) const noexcept{
      return std::get<SIZE>((*this)[find(u)]);
   }
private:
   enum {PARENT, SIZE};
};

#endif // DISJOINT_SETS_HH
