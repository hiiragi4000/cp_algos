// Verified in the problem https://judge.yosupo.jp/problem/addition_of_big_integers
#include"arbitrary_prec.hh"
#include<iostream>
#include<string>
#include<utility>
#include<cstdio>
using namespace std;

int main(){
   int T;
   cin >> T;
   while(T-- > 0){
      string sa, sb;
      cin >> sa >> sb;
      // stob(char const *s, size_t *pos = nullptr):
      // Interprets a BigInt in the string |s|.
      // First it discards whitespaces (identified by std::isspace) until the first non-whitespace character is found.
      // Then it reads a BigInt consists of the following parts:
      // 1. (optional) plus and minus signs
      // 2. a (non-empty) sequence of digits
      // If |pos| is not nullptr, the number of characters consumed is stored in |*pos|.
      // If no conversion can be performed, an |std::invalid_argument| is thrown.
      BigInt a = stob(sa.data()), b = stob(sb.data());
      // A BigInt can be cast into an std::string for output.
      puts(string(a+b).data());
   }
}
