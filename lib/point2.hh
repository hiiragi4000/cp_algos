#ifndef POINT2_HH
#define POINT2_HH

template<typename T> struct Point2{
   T x, y;
   Point2() = default;
   constexpr Point2(T const &x, T const &y): x(x), y(y){}
   template<typename T2> constexpr explicit operator Point2<T2>() const{
      return {x, y};
   }
};
template<typename T1, typename T2>
constexpr bool operator<(Point2<T1> const &p, Point2<T2> const &q){
   return p.x==q.x? p.y<q.y: p.x<q.x;
}
template<typename T1, typename T2>
constexpr bool operator==(Point2<T1> const &p, Point2<T2> const &q){
   return p.x==q.x && p.y==q.y;
}
template<typename T1, typename T2>
constexpr auto operator*(Point2<T1> const &p, Point2<T2> const &q){
   return p.x*q.x + p.y*q.y;
}
template<typename T1, typename T2>
constexpr auto operator%(Point2<T1> const &p, Point2<T2> const &q){
   return p.x*q.y - p.y*q.x;
}

template<typename T>
struct Mat2{
   T e11, e12, e21, e22;
   Mat2() = default;
   constexpr Mat2(T const &e11, T const &e12, T const &e21, T const &e22): e11(e11), e12(e12), e21(e21), e22(e22){}
   template<typename T2> constexpr explicit operator Mat2<T2>() const{
      return {e11, e12, e21, e22};
   }
};
template<typename T1, typename T2>
constexpr auto operator*(Mat2<T1> const &A, Point2<T2> const &v) -> Point2<decltype(A.e11*v.x)>{
   return {A.e11*v.x+A.e12*v.y, A.e21*v.x+A.e22*v.y};
}

#endif // POINT2_HH
