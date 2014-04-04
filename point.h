#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>

template <class T>
struct Point2d {
	T x,y;
  Point2d(T X, T Y);
  Point2d();
};

template <class T>
Point2d<T>::Point2d() {
  x = y = 0;
}

template <class T>
Point2d<T>::Point2d(T X, T Y) {
  x = X;
  y = Y;
}

template <class T>
Point2d<T> operator+(Point2d<T> a, Point2d<T> b) {
  return Point2d<T>(a.x+b.x, a.y+b.y);
}

template <class T>
Point2d<T> operator-(Point2d<T> a, Point2d<T> b) {
  return Point2d<T>(a.x-b.x, a.y-b.y);
}

template <class T>
Point2d<T> operator*(T sf, Point2d<T> a) {
  return Point2d<T>(sf*a.x, sf*a.y);
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Point2d<T>& p) {
  return os << "(" << p.x << "," << p.y << ")";
}

template <class T>
T dot(const Point2d<T>& a, const Point2d<T>& b) {
  return (a.x*b.x) + (a.y*b.y);
}

template <class T>
T cross_z(const Point2d<T>& a, const Point2d<T>& b) {
  return (a.x*b.y) - (a.y*b.x);
}

template <class T>
void segment_intersection(const Point2d<T>& p11, const Point2d<T>& p12, 
                          const Point2d<T>& p21, const Point2d<T>& p22,
                          bool& do_cross, T& t1, T& t2, Point2d<T>& cross_coords) {
  Point2d<T> delta_vec_1 = p12 - p11;
  Point2d<T> delta_vec_2 = p22 - p21;
  
  //solve the vector equation
  //p11 + t1*delta_vec_1 = p21 + t2*delta_vec_2
  Point2d<T> RHS = p21 - p11;
  T det = delta_vec_1.x*(-delta_vec_2.y) - (-delta_vec_2.x)*delta_vec_1.y;
  if (det == 0) {
    do_cross = false;
    return;
  }
  T a11 = (-delta_vec_2.y);
  T a12 = -(-delta_vec_2.x);
  T a21 = -delta_vec_1.y;
  T a22 = delta_vec_1.x;
  Point2d<T> ans(a11*RHS.x + a12*RHS.y, a21*RHS.x + a22*RHS.y);
  ans =  (1/det) * ans;
  
  t1 = ans.x;
  t2 = ans.y;
  if (t1 > 1 || t1 < 0 || t2 > 1 || t2 < 0) {
    do_cross = false;
    return;
  }
  cross_coords = p11 + t1*delta_vec_1;
  do_cross = true;

}


template <class T>
T max(T a, T b) {
  return (a > b ? a : b);
}

#endif