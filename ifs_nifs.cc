nBall::nBall() {
  center = 0;
  radius = 0;
  one = 1;
}

nBall::nBall(cpx c, double r, cpx o) {
  center = c;
  radius = r;
  one = o;
}

bool nBall::is_disjoint(const cpx& ll, const cpx& ur) const {
  return Ball(center, radius).is_disjoint(ll, ur);
}

bool nBall::is_contained(const cpx& ll, const cpx& ur)  const {
  return Ball(center, radius).is_contained(ll,ur);
}



nIFS::nIFS() {
  c = 0;
  centers = std::vector<cpx>();
}

nIFS::nIFS(int n, cpx C) {
  c = C;
  set_std_centers(n);
}
  
void nIFS::set_param(cpx C) {
  c = C;
}

void nIFS::set_centers(const std::vector<cpx>& new_centers) {
  centers = new_centers;
}

void nIFS::set_std_centers(int n) {
  centers.resize(n);
  double angle = 2*3.14159265358979*(1.0/double(n));
  for (int i=0; i<(int)n; ++i) {
    centers[i] = exp(cpx(0,1.0)*angle*(double)i) / (1.0-c);
  }
}
  
double nIFS::minimal_initial_radius() {
  double lb = -1;
  for (int i=0; i<(int)centers.size(); ++i) {
    double lbi = abs(centers[i]*(1.0-c)/(1.0-abs(c)));
    if (lb < 0 || lbi < lb) {
      lb = lbi;
    }
  }
  return lb;
}
  
nBall nIFS::act_on_left(int i, const nBall& b) {
  cpx cent = centers[i];
  return nBall( c*(b.center - cent) + cent, abs(c*b.radius), c*b.one);
}

nBall nIFS::act_on_right(int i, const nBall& b) {
  cpx sub_center = b.center + centers[i]*b.one;
  return nBall( sub_center, abs(c*b.radius), c*b.one);
}










