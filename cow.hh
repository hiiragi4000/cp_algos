#ifndef COW_HH
#define COW_HH

#include<array>
#include<utility>
#include<vector>

template<typename T>
struct CowVla{
   CowVla() = default;
   explicit CowVla(int n, T const &x = T()): n(n){
      auto f = [&x](int) -> decltype(auto){
         return x;
      };
      r.push_back(cur = ver0_build(f, 0, n-1));
   }
   template<typename It> CowVla(It b, int n): n(n){
      auto f = [&b](int i) -> decltype(auto){
         return b[i];
      };
      r.push_back(cur = ver0_build(f, 0, n-1));
   }
   constexpr int size() const noexcept{
      return n;
   }
   void load(int ver) noexcept{
      cur = r[ver];
   }
   int save(){
      int ver = r.size();
      r.push_back(cur);
      return ver;
   }
   T const &operator[](int i) const noexcept{
      int res = cur, l = 0, r = n-1, m;
      while((m=l+(r-l)/2) != i){
         if(i < m){
            res = a[res].ch[0];
            r = m-1;
         }else{
            res = a[res].ch[1];
            l = m+1;
         }
      }
      return a[res].val;
   }
   void write(int i, T const &x){
      cur = write_impl(i, x, cur, 0, n-1);
   }
private:
   struct Node{
      std::array<int, 2> ch{};
      T val;
      Node() = default;
      explicit Node(T const &x): val(x){}
      explicit Node(T &&x): val(std::move(x)){}
   };
   int n=0, cur=0;
   std::vector<int> r;
   std::vector<Node> a;
   template<typename FnT>
   int ver0_build(FnT f, int l, int r){
      int m = l+(r-l)/2, res = a.size();
      a.emplace_back(f(m));
      if(l < m) a[res].ch[0] = ver0_build(f, l, m-1);
      if(m < r) a[res].ch[1] = ver0_build(f, m+1, r);
      return res;
   }
   int write_impl(int i, T const &x, int p, int l, int r){
      int m = l+(r-l)/2, res = a.size();
      if(i == m){
         a.emplace_back(x);
         a[res].ch = a[p].ch;
      }else{
         a.push_back(a[p]);
         if(i < m) a[res].ch[0] = write_impl(i, x, a[p].ch[0], l, m-1);
         else a[res].ch[1] = write_impl(i, x, a[p].ch[1], m+1, r);
      }
      return res;
   }
};

#endif
