// Verified in the problem https://judge.yosupo.jp/problem/range_kth_smallest
#include"range_query.hh"
#include<vector>
#include<cstdio>
using namespace std;

int main(){
   int n, q;
   scanf("%d%d", &n, &q);
   vector<int> a(n);
   for(int i=0; i<=n-1; ++i){
      scanf("%d", a.data()+i);
   }
   // QsortTree::QsortTree(it, n, lt): build a quicksort tree from [it, it+n-1] with the comparator lt.
   // By default std::less<>{} (aka operator<) is used.
   QsortTree qst(a.cbegin(), n);
   while(q-- > 0){
      int l, r, k;
      scanf("%d%d%d", &l, &r, &k);
      // QsortTree::query(l, r, k): returns the index that is the kth-smallest between *(it+l), *(it+l+1), ..., *(it+r).
      // If there are multiple elements that are the kth-smallest, returns any of them.
      // Note that a QsortTree only extracts the size relations between *it, *(it+1), ..., *(it+n-1).
      // To get the element itself, the user should refer to the arguments that initialized the QsortTree.
      printf("%d\n", a[qst.query(l, r-1, k+1)]);
   }
}
