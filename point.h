#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>
#include <algorithm>

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
Point2d<T> operator/(const Point2d<T>& a, T f) {
  return Point2d<T>(a.x/f, a.y/f);
}

template <class T>
bool operator==(const Point2d<T>& a, const Point2d<T>& b) {
  return a.x == b.x && a.y == b.y;
}

template <class T>
bool operator!=(const Point2d<T>& a, const Point2d<T>& b) {
  return a.x != b.x || a.y != b.y;
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



//3d point

template <class T>
struct Point3d {
        T x,y,z;
  Point3d(T X, T Y, T Z);
  Point3d();
};

template <class T>
Point3d<T>::Point3d() {
  x = y = z = 0;
}

template <class T>
Point3d<T>::Point3d(T X, T Y, T Z) {
  x = X;
  y = Y;
  z = Z;
}

template <class T>
Point3d<T> operator+(Point3d<T> a, Point3d<T> b) {
  return Point3d<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template <class T>
Point3d<T> operator-(Point3d<T> a, Point3d<T> b) {
  return Point2d<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template <class T>
Point3d<T> operator*(T sf, Point3d<T> a) {
  return Point3d<T>(sf*a.x, sf*a.y, sf*a.z);
}

template <class T>
bool operator==(const Point3d<T>& a, const Point3d<T>& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

template <class T>
bool operator!=(const Point3d<T>& a, const Point3d<T>& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z;
}


template <class T>
std::ostream& operator<<(std::ostream& os, const Point3d<T>& p) {
  return os << "(" << p.x << "," << p.y << "," << p.z << ")";
}

//std::ostream& operator<<(std::ostream& os, const unsigned char c) {
//  return os << (unsigned int)c;
//}

template <class T>
T dot(const Point3d<T>& a, const Point3d<T>& b) {
  return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
}


//4d I guess why not
template <class T>
struct Point4d {
        T x,y,z,w;
  Point4d(T X, T Y, T Z, T W);
  Point4d();
};

template <class T>
Point4d<T>::Point4d() {
  x = y = z = w = 0;
}

template <class T>
Point4d<T>::Point4d(T X, T Y, T Z, T W) {
  x = X;
  y = Y;
  z = Z;
  w = W;
}

template <class T>
Point4d<T> operator+(const Point4d<T>& a, const Point4d<T>& b) {
  return Point4d<T>(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

template <class T>
Point4d<T> operator-(const Point4d<T>& a, const Point4d<T>& b) {
  return Point2d<T>(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}

template <class T>
Point4d<T> operator*(T sf, const Point4d<T>& a) {
  return Point4d<T>(sf*a.x, sf*a.y, sf*a.z, sf*a.w);
}

template <class T>
bool operator==(const Point4d<T>& a, const Point4d<T>& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

template <class T>
bool operator!=(const Point4d<T>& a, const Point4d<T>& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}


template <class T>
std::ostream& operator<<(std::ostream& os, const Point4d<T>& p) {
  return os << "(" << p.x << "," << p.y << "," << p.z << "," << p.w << ")";
}

//std::ostream& operator<<(std::ostream& os, const unsigned char c) {
//  return os << (unsigned int)c;
//}

template <class T>
T dot(const Point4d<T>& a, const Point4d<T>& b) {
  return (a.x*b.x) + (a.y*b.y) + (a.z*b.z) + (a.w*b.w);
}


//Nd -- it's made to generalize 4d, so it has x,y,z,w, then the rest
//it's an error if n<
template <int n, class T>
struct PointNd {
  T x,y,z,w;
  std::vector<T> the_rest;
        
  PointNd(const T& val);
  PointNd(const std::vector<T>& v);
  PointNd();
  T& operator[](int i);
};

template <int n, class T>
PointNd<n,T>::PointNd() {
  the_rest.resize(n-4);
  x = y = z = w = 0;
  for (int i=0; i<n-4; ++i) {
    the_rest[i] = 0;
  }
}

template <int n, class T>
PointNd<n,T>::PointNd(const T& val) {
  the_rest.resize(n-4);
  x = y = z = w = val;
  for (int i=0; i<n-4; ++i) {
    the_rest[i] = val;
  }
}

template <int n, class T>
PointNd<n,T>::PointNd(const std::vector<T>& v) {
  the_rest.resize(n-4);
  x = v[0]; y = v[1]; z = v[2]; w = v[3];
  std::copy(v.begin()+4, v.end(), the_rest.begin());
}


template <int n, class T>
std::ostream& operator<<(std::ostream& os, const PointNd<n,T>& p) {
  os << "(" << p.x << "," << p.y << "," << p.z << "," << p.w;
  for (int i=0; i<n-4; ++i) {
    os << "," << p.the_rest[i];
  }
  return os << ")";
}

template <int n, class T>
T& PointNd<n,T>::operator[](int i) {
  if (i<4) {
    return (i==0 ? x : (i==1 ? y : (i==2 ? z : w)));
  } else {
    return the_rest[i-4];
  }
}























#endif