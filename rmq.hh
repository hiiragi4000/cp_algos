#ifndef RMQ_HH
#define RMQ_HH

#include<array>
#include<functional>
#include<iterator>
#include<vector>

#define LB 6
#define P 67

struct Rmq{
   Rmq() = default;
   template<typename It, typename Lt=std::less<typename std::iterator_traits<It>::value_type>>
   Rmq(It a, int n, Lt const &lt={}){
      if(n <= 0) return;
      int n2 = ((n-1)>>LB)+1;
      lg.resize(n2+1);
      for(int i=2; i<=n2; ++i){
         lg[i] = lg[i>>1]+1;
      }
      b.resize(n);
      st.resize((lg[n2]+1)*n2);
      for(int i=0; i<n; i+=1<<LB){
         i8 stk[1<<LB];
         u64 mask = 0;
         int top = 0, min_id = i;
         for(int j=0; (i|j)<n && j<1<<LB; ++j){
            while(top && !lt(a[i|stk[top-1]], a[i|j])){
               mask &= ~(1ull<<stk[--top]);
            }
            b[i|j] = mask |= 1ull<<(stk[top++]=j);
            if(lt(a[i|j], a[min_id])){
               min_id = i|j;
            }
         }
         st[i>>LB] = min_id;
      }
      for(int i=1; i<=lg[n2]; ++i){
         for(int j=0; j+(1<<i)<=n2; ++j){
            int k = st[(i-1)*n2+j], l = st[(i-1)*n2+j+(1<<(i-1))];
            st[i*n2+j] = lt(a[l], a[k])? l: k;
         }
      }
   }
   std::array<int, 4> query(int l, int r) const{
      int l1 = l>>LB, l2 = l&(1<<LB)-1, r1 = r>>LB;
      if(l1 == r1){
         u64 mask = b[r] & ~((1ull<<l2)-1);
         u64 bit = mask & ~mask+1;
         int k = l1<<LB | LG_LOWBIT[bit%P];
         return {k, k, k, k};
      }
      u64 mask = b[((l1+1)<<LB)-1] & ~((1ull<<l2)-1);
      u64 bit = mask & ~mask+1;
      int k1 = l1<<LB | LG_LOWBIT[bit%P];
      int k4 = r1<<LB | LG_LOWBIT[(b[r]&~b[r]+1)%P];
      if(l1+1 == r1){
         return {k1, k4, k4, k4};
      }
      int m = lg[r1-l1-1], n2 = ((b.size()-1)>>LB)+1;
      int k2 = st[m*n2+l1+1];
      int k3 = st[m*n2+r1-(1<<m)];
      return {k1, k2, k3, k4};
   }
private:
   // LG_LOWBIT[2^n mod P] = n for n in [0, P-1)
   // TODO: initialize LG_LOWBIT from a constexpr fn instead
   static constexpr signed char LG_LOWBIT[P] = {0, 0, 1, 39, 2, 15, 40, 23, 3, 12, 16, 59, 41, 19, 24, 54, 4, 64, 13, 10, 17, 62, 60, 28, 42, 30, 20, 51, 25, 44, 55, 47, 5, 32, 65, 38, 14, 22, 11, 58, 18, 53, 63, 9, 61, 27, 29, 50, 43, 46, 31, 37, 21, 57, 52, 8, 26, 49, 45, 36, 56, 7, 48, 35, 6, 34, 33};
   using u64 = unsigned long long;
   using i8 = signed char;
   std::vector<u64> b;
   std::vector<i8> lg;
   std::vector<int> st;
};

#undef P
#undef LB

#endif
