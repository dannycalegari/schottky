
LinearMap::LinearMap() {
  a=b=c=d=0;
}

LinearMap::LinearMap(double A, double B, double C, double D) {
  a = A; b = B; c = C; d = D;
}

LinearMap LinearMap::operator*(const LinearMap& other) const {
  return LinearMap( a*other.a + b*other.c,  a*other.b + b*other.d,
                    c*other.a + d*other.c,  c*other.b + d*other.d );
}

Point2d<double> LinearMap::operator()(const Point2d<double>& X) const {
  return Point2d<double>( a*X.x + b*X.y, c*X.x + d*X.y );
}



AffineMap::AffineMap() {
  A = LinearMap();
  t = Point2d<double>(0,0);
}

AffineMap::AffineMap(double a, double b, double c, double d, double x, double y) {
  A = LinearMap(a,b,c,d);
  t = Point2d<double>(x,y);
}


AffineMap::AffineMap(const LinearMap& a, const Point2d<double>& T) {
  A = a;
  t = T;
}

AffineMap AffineMap::operator*(const AffineMap& other) const {
  return AffineMap( A*other.A, A(other.t) + t );
}

Point2d<double> AffineMap::operator()(const Point2d<double>& X) const {
  return A(X) + t;
}


AffineMap ifs2d::semigroup_element(const std::vector<int>& gen_word) {
  AffineMap A = gens[gen_word[0]];
  for (int i=1; i<(int)gen_word.size(); ++i) {
    A = A*gens[gen_word[i]];
  }
  return A;
}

AffineMap ifs2d::semigroup_element(int list, int n) {
  if ((int)gens.size() > 2 || n==0) return AffineMap();
  
  AffineMap A = gens[list&1];
  for (int i=1; i<n; ++i) {
    A = gens[(list>>i)&1]*A;
  }
  return A;
}

Point4d<double> point_as_weighted_average_in_box(const Point2d<double>& point_target, 
                                                 const Point2d<double>& center, 
                                                 double radius) {
  return Point4d<double>();
}







