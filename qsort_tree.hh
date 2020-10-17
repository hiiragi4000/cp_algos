#ifndef QSORT_TREE_HH
#define QSORT_TREE_HH

#include<algorithm>
#include<functional>
#include<iterator>
#include<numeric>
#include<vector>

struct QsortTree{
public:
   QsortTree() = default;
   template<typename It, typename Lt=std::less<typename std::iterator_traits<It>::value_type>>
   QsortTree(It a, int n, Lt const &lt={}): si(n), val(buf_size(n)), sl(buf_size(n)), pos(4*n){
      iota(si.begin(), si.end(), 0);
      sort(si.begin(), si.end(), [&](int i, int j){
         return lt(a[i], a[j]);
      });
      for(int i=0; i<n; ++i) val[si[i]+1] = i;
      int top = n+1; pos[1] = 0;
      build(1, 0, n-1, top);
   }
   int query(int l, int r, int k) const{
      return si[m_query(l+1, r+1, k, 1)];
   }
private:
   std::vector<int> si, val, sl, pos;
   static constexpr int buf_size(int n) noexcept{
      int lg = 0;
      for(int t=n; t>1; t>>=1) ++lg;
      if(n > 1<<lg) ++lg;
      return n*(lg+3);
   }
   void build(int i, int l, int r, int &top){
      if(l == r) return;
      int m = l+(r-l)/2;
      pos[2*i] = top; top += m-l+2;
      pos[2*i+1] = top; top += r-m+1;
      sl[pos[i]] = 0;
      for(int j=1, k1=0, k2=0; j<=r-l+1; ++j){
         if(val[pos[i]+j] <= m){
            sl[pos[i]+j] = sl[pos[i]+j-1]+1;
            val[pos[2*i]+ ++k1] = val[pos[i]+j];
         }else{
            sl[pos[i]+j] = sl[pos[i]+j-1];
            val[pos[2*i+1]+ ++k2] = val[pos[i]+j];
         }
      }
      build(2*i, l, m, top);
      build(2*i+1, m+1, r, top);
   }
   int m_query(int l, int r, int k, int i) const{
      if(l == r) return val[pos[i]+l];
      if(k <= sl[pos[i]+r]-sl[pos[i]+l-1]){
         return m_query(sl[pos[i]+l-1]+1, sl[pos[i]+r], k, 2*i);
      }
      return m_query(l-sl[pos[i]+l-1], r-sl[pos[i]+r], k-sl[pos[i]+r]+sl[pos[i]+l-1], 2*i+1);
   }
};

#endif
