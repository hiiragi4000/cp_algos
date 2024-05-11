// Verified in the problem https://judge.yosupo.jp/problem/addition_of_big_integers
#include"arbitrary_prec.hh"
#include<iostream>
#include<string>
#include<cstdio>
using namespace std;

int main(){
   int T;
   cin >> T;
   while(T-- > 0){
      string sa, sb;
      cin >> sa >> sb;
      // BigInt stob(char const *s, size_t *pos=nullptr, int base=10):
      // Interprets a BigInt from the string |s| in base |base|.
      // The base must be in {0, 2, 3, ..., 36}; otherwise, an |std::invalid_argument| is thrown.
      // First the function discards whitespaces (identified by std::isspace) until the first non-whitespace character is found.
      // Then it reads a BigInt consists of the following parts:
      // 1. (optional) plus and minus signs
      // 2. (optional) prefix 0 or 0x/0X, indicating octal or hexadecimal base (applies only when the base is 0)
      // 3. a (non-empty) sequence of digits
      // The set of valid digits for base-2 integers is {0, 1}, for base-3 integers is {0, 1, 2}, and so on.
      // For bases larger than 10, the set of valid digits contains alphabets with the cases ignored.
      // If |pos| is not nullptr, the number of characters consumed is stored in |*pos|.
      // If no conversion can be performed, an |std::invalid_argument| is thrown.
      BigInt a = stob(sa.data()), b = stob(sb.data());
      // A BigInt can be cast into an std::string for output.
      puts(string(a+b).data());
   }
}
