#ifndef RANGE_QUERY_HH
#define RANGE_QUERY_HH

#include"util.hh"
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

template<typename T, typename CompT=std::plus<T>> struct Fenwick{
   Fenwick() = default;
   explicit Fenwick(int n, T const &e={}, CompT const &op={}): b(n+1u, e), op(op){}
   T sum(int i) const noexcept{
      T res = b[0];
      for(size_t j=i+1; j; j-=j&~j+1){
         res = op(res, b[j]);
      }
      return res;
   }
   void add(int i, T const &dx) noexcept{
      for(size_t j=i+1; j<b.size(); j+=j&~j+1){
         b[j] = op(b[j], dx);
      }
   }
private:
   std::vector<T> b;
   CompT op;
};

struct RangeMinimum{
   RangeMinimum() = default;
   template<typename It, typename Lt=std::less<DerefType<It>>>
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
   static constexpr int LB = 6, P = 67;
   static constexpr std::array<I8, P> lg_lowbit_build(){
      std::array<I8, P> res{0};
      int p = 1;
      for(int i=0; i<P-1; ++i){
         res[p] = i;
         p = 2*p%P;
      }
      return res;
   }
   // LG_LOWBIT[2^n mod 67] = n for n in [0, 65]
   static std::array<I8, P> LG_LOWBIT;
   std::vector<U64> b;
   std::vector<I8> lg;
   std::vector<int> st;
};
inline std::array<I8, RangeMinimum::P> RangeMinimum::LG_LOWBIT = RangeMinimum::lg_lowbit_build();

template<typename T, typename MonT, typename CompT, typename ActT>
struct VlaSegMonAct{
   VlaSegMonAct() = default;
   VlaSegMonAct(std::vector<T> a0, MonT const &id, CompT const &comp, ActT const &act): id(id), op(comp), act(act), a(move(a0)), tag(a.size(), id){}
   void seg_act(int l, int r, MonT const &f){
      seg_act_impl(l, r, f, 0, (int)a.size()-1);
   }
   T operator[](int i) const{
      return bracket_impl(i, 0, (int)a.size()-1);
   }
private:
   MonT id;
   CompT op;
   ActT act;
   std::vector<T> a;
   std::vector<MonT> tag;
   void seg_act_impl(int l0, int r0, MonT const &f, int l, int r){
      if(r<l0 || r0<l) return;
      int m = l+(r-l)/2;
      if(l0<=l && r<=r0){
         tag[m] = op(f, tag[m]);
         return;
      }
      a[m] = act(tag[m], a[m]);
      if(l0<=m && m<=r0){
         a[m] = act(f, a[m]);
      }
      if(l <= m-1){
         int lm = l+(m-l-1)/2;
         tag[lm] = op(tag[m], tag[lm]);
         seg_act_impl(l0, r0, f, l, m-1);
      }
      if(m+1 <= r){
         int rm = (m+1)+(r-m-1)/2;
         tag[rm] = op(tag[m], tag[rm]);
         seg_act_impl(l0, r0, f, m+1, r);
      }
      tag[m] = id;
   }
   T bracket_impl(int i0, int l, int r) const{
      int m = l+(r-l)/2;
      return act(tag[m], i0<m? bracket_impl(i0, l, m-1): i0==m? a[m]: bracket_impl(i0, m+1, r));
   }
};

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

template<typename GrpT, typename CompT, typename InvT>
struct TreePathVertexComp{
   TreePathVertexComp() = default;
   template<typename It>
   TreePathVertexComp(std::vector<int> const *g, int n, It x, GrpT const &e, CompT const &comp={}, InvT const &inv={}): op(comp), inv(inv), lca(g, n), pos(n+1u), pa(n), sz(n, 1){
      // assert(n >= 1);
      using namespace std::placeholders;
      tr.reserve(n);
      pa[0] = n; pos[n] = n;
      std::vector<GrpT> lp(n+1u, e), rp(n+1u, e);
      dfs(g, x, 0, lp.data(), rp.data());
      auto rev_op = bind(op, _2, _1);
      left = {move(lp), e, rev_op, rev_op};
      right = {move(rp), e, op, op};
   }
   void set_vertex(int u, GrpT const &g){
      auto a = op(op(inv(left[pos[u]]), g), left[pos[pa[u]]]);
      left.seg_act(pos[u], pos[u]+sz[u]-1, a);
      auto b = op(op(right[pos[pa[u]]], g), inv(right[pos[u]]));
      right.seg_act(pos[u], pos[u]+sz[u]-1, b);
   }
   GrpT query(int u, int v) const{
      int a = lca.query(u, v);
      return op(op(op(left[pos[u]], inv(left[pos[pa[a]]])), inv(right[pos[a]])), right[pos[v]]);
   }
private:
   using RevCompT = std::function<GrpT(GrpT const&, GrpT const&)>;
   CompT op;
   InvT inv;
   LowestCommonAncestor lca;
   std::vector<int> tr, pos, pa, sz;
   VlaSegMonAct<GrpT, GrpT, RevCompT, RevCompT> left;
   VlaSegMonAct<GrpT, GrpT, CompT, CompT> right;
   template<typename It> void dfs(std::vector<int> const *g, It x, int u, GrpT *lp, GrpT *rp){
      pos[u] = tr.size();
      lp[tr.size()] = op(x[u], lp[pos[pa[u]]]);
      rp[tr.size()] = op(rp[pos[pa[u]]], x[u]);
      tr.push_back(u);
      for(int v: g[u]) if(v != pa[u]){
         pa[v] = u;
         dfs(g, x, v, lp, rp);
         sz[u] += sz[v];
      }
   }
};

template<typename AbT, typename CompT=std::plus<AbT>, typename InvT=std::negate<AbT>>
struct TreePathVertexSum{
   TreePathVertexSum() = default;
   template<typename It>
   TreePathVertexSum(std::vector<int> const *g, int n, It x, AbT const &e=0, CompT const &comp={}, InvT const &inv={}): op(comp), inv(inv), lca(g, n), pos(n+1u), pa(n), sz(n, 1){
      // assert(n >= 1);
      tr.reserve(n);
      pa[0] = n; pos[n] = n;
      std::vector<AbT> s(n+1u, e);
      dfs(g, x, 0, s.data());
      sum = {move(s), e, op, op};
   }
   void add_vertex(int u, AbT const &a){
      sum.seg_act(pos[u], pos[u]+sz[u]-1, a);
   }
   void set_vertex(int u, AbT const &a){
      add_vertex(op(op(inv(sum[pos[u]]), a), sum[pos[pa[u]]]));
   }
   AbT query(int u, int v) const{
      int a = lca.query(u, v);
      return op(op(op(sum[pos[u]], inv(sum[pos[pa[a]]])), inv(sum[pos[a]])), sum[pos[v]]);
   }
private:
   CompT op;
   InvT inv;
   LowestCommonAncestor lca;
   std::vector<int> tr, pos, pa, sz;
   VlaSegMonAct<AbT, AbT, CompT, CompT> sum;
   template<typename It> void dfs(std::vector<int> const *g, It x, int u, AbT *s){
      pos[u] = tr.size();
      s[tr.size()] = op(x[u], s[pos[pa[u]]]);
      tr.push_back(u);
      for(int v: g[u]) if(v != pa[u]){
         pa[v] = u;
         dfs(g, x, v, s);
         sz[u] += sz[v];
      }
   }
};

struct QsortTree{
public:
   QsortTree() = default;
   template<typename It, typename Lt=std::less<DerefType<It>>>
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

struct RangeInversion{
   RangeInversion() = default;
   template<typename It, typename Lt = std::less<DerefType<It>>>
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
      I64 res = b[k1].sinv[r1] + b[k2].pinv[r2] + n_inv(b[k1].ss[r1].data(), bsize-r1, b[k2].sp[r2].data(), r2+1);
      if(k1+1 < k2){
         res += macro_inv[k1+1][k2-1] + dp[l][k2-1] + dp[r][k1+1];
      }
      return res;
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
