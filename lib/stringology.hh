#ifndef STRINGOLOGY_HH
#define STRINGOLOGY_HH

#include<algorithm>
#include<numeric>
#include<tuple>
#include<vector>

struct SuffixArray{
   template<typename It>
   static std::vector<int> build(It s, int n){
      if(n <= 0) return {};
      // assert((*std::min_element(s, s+n) >= 0));
      // assert((*std::max_element(s, s+n) < INT_MAX));
      int m = *std::max_element(s, s+n)+1;
      std::vector<int> sa(n), buf(2*n+m+std::max(n, m));
      std::vector<bool> t_buf(2*n);
      build_impl(s, n, m, sa.data(), buf.data(), t_buf, 0);
      return sa;
   }
private:
   template<typename It>
   static void build_impl(It s, int n, int m, int *sa, int *p, std::vector<bool> &t_buf, int t_offset){
      int *b = p+(n+1)/2, *b2 = b+m, *r = b2+std::max(m, n), *s2 = r+n;
      auto t = [&t_buf, t_offset](int i){
         return t_buf[t_offset+i];
      };
      std::fill(b, b+m, 0);
      for(int i=0; i<n; ++i){
         ++b[s[i]];
      }
      std::partial_sum(b, b+m, b);
      t(n-1) = false;
      for(int i=n-2; i>=0; --i){
         t(i) = s[i]==s[i+1]? t(i+1): s[i]<s[i+1];
      }
      auto isort = [s, n, m, sa, b, b2, &t](int *q, int n2){
         std::fill(sa, sa+n, 0);
         std::copy(b, b+m, b2);
         for(int i=n2-1; i>=0; --i){
            sa[--b2[s[q[i]]]] = q[i];
         }
         b2[0] = 0;
         std::copy(b, b+m-1, b2+1);
         sa[b2[s[n-1]]++] = n-1;
         for(int i=0; i<n; ++i) if(sa[i] && !t(sa[i]-1)){
            sa[b2[s[sa[i]-1]]++] = sa[i]-1;
         }
         std::copy(b, b+m, b2);
         for(int i=n-1; i>=0; --i) if(sa[i] && t(sa[i]-1)){
            sa[--b2[s[sa[i]-1]]] = sa[i]-1;
         }
      };
      int n2 = 0;
      for(int i=0; i<n; ++i){
         if(i>0 && t(i) && !t(i-1)){
            r[i] = n2;
            p[n2++] = i;
         }else r[i] = -1;
      }
      p[n2] = n;
      isort(p, n2);
      auto diff = [s, n](int i, int j, int len){
         while(len-- > 0){
            if(i==n || j==n || s[i++]!=s[j++]) return true;
         }
         return false;
      };
      int m2 = 0;
      for(int i=0, j=-1; i<n; ++i) if(int k=r[sa[i]]; k>=0){
         m2 += j==-1 || diff(p[j], p[k], p[j+1]-p[j]+1);
         s2[k] = m2-1; j = k;
      }
      if(m2 < n2){
         build_impl(s2, n2, m2, sa, b2, t_buf, t_offset+n);
      }else{
         for(int i=0; i<n2; ++i){
            sa[s2[i]] = i;
         }
      }
      for(int i=0; i<n2; ++i){
         s2[i] = p[sa[i]];
      }
      isort(s2, n2);
   }
};

struct LcpArray{
   template<typename It>
   static std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> build_sa_rank_lcp(It s, int n){
      if(n <= 0) return {};
      auto sa = SuffixArray::build(s, n);
      std::vector<int> ra(n);
      for(int i=0; i<n; ++i){
         ra[sa[i]] = i;
      }
      std::vector<int> lcp(n-1);
      for(int i=0, k=0; i<n; ++i){
         if(ra[i] == n-1){
            k = 0; continue;
         }
         int j = sa[ra[i]+1];
         while(i+k<n && j+k<n && s[i+k]==s[j+k]){
            ++k;
         }
         lcp[ra[i]] = k;
         if(k > 0) --k;
      }
      return {sa, ra, lcp};
   }
};

template<typename TransFnT> struct SuffixAutomaton{
   struct State{
      int lnk=-1, len{};
      TransFnT next{};
   };
   template<typename It>
   static std::vector<State> build(It s, int n){
      std::vector<State> a;
      a.reserve(n>=2? 2*n-1: n+1);
      a.emplace_back();
      int cur = 0;
      for(int i=0; i<n; ++i){
         cur = extend(a, cur, *s++);
      }
      return a;
   }
private:
   static int extend(std::vector<State> &a, int last, int c){
      int cur = a.size();
      a.emplace_back();
      a[cur].len = a[last].len+1;
      int p = last;
      while(p>=0 && !a[p].next[c]){
         a[p].next[c] = cur;
         p = a[p].lnk;
      }
      if(p == -1){
         a[cur].lnk = 0;
      }else{
         int q = a[p].next[c];
         if(a[p].len+1 == a[q].len){
            a[cur].lnk = q;
         }else{
            int q2 = a.size();
            a.push_back(a[q]);
            a[q2].len = a[p].len+1;
            while(p>=0 && a[p].next[c]==q){
               a[p].next[c] = q2;
               p = a[p].lnk;
            }
            a[cur].lnk = a[q].lnk = q2;
         }
      }
      return cur;
   }
};

#endif
