#ifndef COW_RC_HH
#define COW_RC_HH

#include<array>
#include<memory>
#include<utility>
#include<vector>

template<typename T>
struct VlaCow{
   VlaCow() = default;
   explicit VlaCow(int n, T const &x = T()): n(n){
      auto f = [&x](int) -> decltype(auto){
         return x;
      };
      r.push_back(cur = ver0_build(f, 0, n-1));
   }
   template<typename It> VlaCow(It b, int n): n(n){
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
      Node *p = cur.get();
      int l = 0, r = n-1, m;
      while((m=l+(r-l)/2) != i){
         if(i < m){
            p = p->ch[0].get();
            r = m-1;
         }else{
            p = p->ch[1].get();
            l = m+1;
         }
      }
      return p->val;
   }
   void write(int i, T const &x){
      cur = write_impl(i, x, cur.get(), 0, n-1);
   }
private:
   struct Node{
      std::array<std::shared_ptr<Node>, 2> ch{};
      T val;
      Node() = default;
      explicit Node(T const &x): val(x){}
      explicit Node(T &&x): val(std::move(x)){}
   };
   int n = 0;
   std::shared_ptr<Node> cur;
   std::vector<std::shared_ptr<Node>> r;
   template<typename FnT>
   std::shared_ptr<Node> ver0_build(FnT f, int l, int r){
      int m = l+(r-l)/2;
      auto p = std::make_shared<Node>(f(m));
      if(l < m) p->ch[0] = ver0_build(f, l, m-1);
      if(m < r) p->ch[1] = ver0_build(f, m+1, r);
      return p;
   }
   std::shared_ptr<Node> write_impl(int i, T const &x, Node const *p, int l, int r){
      int m = l+(r-l)/2;
      std::shared_ptr<Node> q;
      if(i == m){
         q = std::make_shared<Node>(x);
         q->ch = p->ch;
      }else{
         q = std::make_shared<Node>(*p);
         if(i < m) q->ch[0] = write_impl(i, x, p->ch[0].get(), l, m-1);
         else q->ch[1] = write_impl(i, x, p->ch[1].get(), m+1, r);
      }
      return q;
   }
};

#endif // COW_RC_HH
