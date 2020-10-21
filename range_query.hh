#ifndef RANGE_QUERY_HH
#define RANGE_QUERY_HH

#include<algorithm>
#include<array>
#include<functional>
#include<iterator>
#include<numeric>
#include<utility>
#include<vector>
#include<cmath>

#define I8 signed char
#define I64 long long
#define U64 unsigned long long

#define LB 6
#define P 67

struct RangeMinimum{
   RangeMinimum() = default;
   template<typename It, typename Lt=std::less<typename std::iterator_traits<It>::value_type>>
   RangeMinimum(It a, int n, Lt const &lt={}){
      if(n <= 0) return;
      int n2 = ((n-1)>>LB)+1;
      lg.resize(n2+1);
      for(int i=2; i<=n2; ++i){
         lg[i] = lg[i>>1]+1;
      }
      b.resize(n);
      st.resize((lg[n2]+1)*n2);
      for(int i=0; i<n; i+=1<<LB){
         I8 stk[1<<LB];
         U64 mask = 0;
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
         U64 mask = b[r] & ~((1ull<<l2)-1);
         U64 bit = mask & ~mask+1;
         int k = l1<<LB | LG_LOWBIT[bit%P];
         return {k, k, k, k};
      }
      U64 mask = b[((l1+1)<<LB)-1] & ~((1ull<<l2)-1);
      U64 bit = mask & ~mask+1;
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
   static constexpr std::array<I8, P> lg_lowbit_build(){
      std::array<I8, P> res{0};
      int p = 1;
      for(int i=0; i<P-1; ++i){
         res[p] = i;
         p = 2*p%P;
      }
      return res;
   }
   static std::array<I8, P> LG_LOWBIT;
   std::vector<U64> b;
   std::vector<I8> lg;
   std::vector<int> st;
};

// LG_LOWBIT[2^n mod P] = n for n in [0, P-1)
inline std::array<I8, P> RangeMinimum::LG_LOWBIT = RangeMinimum::lg_lowbit_build();

#undef P
#undef LB

struct LowestCommonAncestor{
   LowestCommonAncestor() = default;
   LowestCommonAncestor(std::vector<int> const *g, int n): fi(n), dep(2*n), vis(2*n){
      int n_step = 0;
      dfs(g, 0, 0, 0, n_step);
      rmq = RangeMinimum(dep.cbegin(), dep.size());
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
   RangeMinimum rmq;
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

struct RangeInversion{
   RangeInversion() = default;
   template<typename It, typename Lt = std::less<typename std::iterator_traits<It>::value_type>>
   RangeInversion(It first, int n, Lt const &lt={}){
      std::vector<int> si(n);
      iota(si.begin(), si.end(), 0);
      sort(si.begin(), si.end(), [&first, &lt](int i, int j){
         if(lt(first[i], first[j])) return true;
         if(lt(first[j], first[i])) return false;
         return i < j;
      });
      int bsize = std::round(std::sqrt(n));
      int len = (n-1)/bsize + 1;
      std::vector<int> a(len*bsize);
      for(int i=0; i<n; ++i) a[si[i]] = i;
      iota(a.begin()+n, a.end(), n);
      b.resize(len);
      macro_inv.resize(len, std::vector<I64>(len));
      for(int i=0; i<len; ++i){
         b[i] = Bucket(bsize);
         std::vector<int> isort;
         isort.reserve(bsize);
         for(int j=0; j<bsize; ++j){
            auto it = upper_bound(isort.cbegin(), isort.cend(), a[i*bsize+j]);
            b[i].pinv[j] = isort.cend() - it;
            isort.insert(it, a[i*bsize+j]);
            b[i].sp[j] = isort;
         }
         partial_sum(b[i].pinv.cbegin(), b[i].pinv.cend(), b[i].pinv.begin());
         isort.clear();
         for(int j=bsize-1; j>=0; --j){
            auto it = lower_bound(isort.cbegin(), isort.cend(), a[i*bsize+j]);
            b[i].sinv[j] = it - isort.cbegin();
            isort.insert(it, a[i*bsize+j]);
            b[i].ss[j] = isort;
         }
         partial_sum(b[i].sinv.crbegin(), b[i].sinv.crend(), b[i].sinv.rbegin());
         for(int j=0; j<bsize; ++j){
            b[i].si[lower_bound(isort.cbegin(), isort.cend(), a[i*bsize+j])-isort.cbegin()] = j;
         }
         macro_inv[i][i] = b[i].sinv[0];
      }
      dp.resize(len*bsize, std::vector<int>(len));
      for(int i=0; i<len-1; ++i){
         macro_inv[i][i+1] = b[i].sinv[0] + b[i+1].sinv[0] + b_inv(i, i+1);
      }
      for(int i=2; i<len; ++i) for(int j=0; j<len-i; ++j){
         macro_inv[j][j+i] = macro_inv[j][j+i-1] + macro_inv[j+1][j+i] + b_inv(j, j+i) - macro_inv[j+1][j+i-1];
      }
      for(int i=0; i<len; ++i){
         for(int j=0; j<bsize; ++j){
            for(int k=i-2; k>=0; --k){
               dp[i*bsize+j][k] += dp[i*bsize+j][k+1];
            }
            for(int k=i+2; k<len; ++k){
               dp[i*bsize+j][k] += dp[i*bsize+j][k-1];
            }
         }
         for(int j=1; j<bsize; ++j){
            for(int k=0; k<=i-1; ++k){
               dp[i*bsize+j][k] += dp[i*bsize+j-1][k];
            }
         }
         for(int j=bsize-2; j>=0; --j){
            for(int k=i+1; k<len; ++k){
               dp[i*bsize+j][k] += dp[i*bsize+j+1][k];
            }
         }
      }
   }
   I64 query(int l, int r) const{
      int bsize = bucket_size();
      int k1 = l/bsize, r1 = l%bsize, k2 = r/bsize, r2 = r%bsize;
      if(k1 == k2){
         if(r1 == 0) return b[k1].pinv[r2];
         if(r2 == bsize-1) return b[k1].sinv[r1];
         return b[k1].pinv[r2] + b[k1].sinv[r1] + n_inv(b[k1].sp[r1-1].data(), r1, b[k1].ss[r2+1].data(), bsize-r2-1) - macro_inv[k1][k1];
      }
      if(k1+1 == k2){
         return b[k1].sinv[r1] + b[k2].pinv[r2] + n_inv(b[k1].ss[r1].data(), bsize-r1, b[k2].sp[r2].data(), r2+1);
      }
      return b[k1].sinv[r1] + b[k2].pinv[r2] + n_inv(b[k1].ss[r1].data(), bsize-r1, b[k2].sp[r2].data(), r2+1) + macro_inv[k1+1][k2-1] + dp[l][k2-1] + dp[r][k1+1];
   }
private:
   struct Bucket{
      std::vector<int> si, pinv, sinv;
      std::vector<std::vector<int>> sp, ss;
      Bucket() = default;
      explicit Bucket(int bsize): si(bsize), pinv(bsize), sinv(bsize), sp(bsize), ss(bsize){}
   };
   std::vector<Bucket> b;
   std::vector<std::vector<I64>> macro_inv;
   std::vector<std::vector<int>> dp;
   int bucket_size() const noexcept{
      if(b.empty()) return 0;
      return b[0].si.size();
   }
   int b_inv(int i, int j){
      int bsize = bucket_size(), res = 0;
      for(int k1=0, k2=0; k1<bsize || k2<bsize; ){
         if(k2==bsize || k1<bsize && b[i].ss[0][k1]<=b[j].ss[0][k2]){
            dp[i*bsize+b[i].si[k1]][j] = k2;
            res += k2; ++k1;
         }else{
            dp[j*bsize+b[j].si[k2]][i] = bsize-k1;
            ++k2;
         }
      }
      return res;
   }
   static int n_inv(int const *a, int n, int const *b, int m) noexcept{
      int res = 0;
      for(int i=0, j=0; i<n || j<m; ){
         if(j==m || i<n && a[i]<=b[j]){
            res += j; ++i;
         }else ++j;
      }
      return res;
   }
};

#undef U64
#undef I64
#undef I8

#endif
