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


bool nBall::is_disjoint(const nBall& other) const {
  return abs( center - other.center ) > radius + other.radius;
}

bool nBall::is_disjoint(const cpx& ll, const cpx& ur) const {
  return Ball(center, radius).is_disjoint(ll, ur);
}

bool nBall::is_contained(const cpx& ll, const cpx& ur)  const {
  return Ball(center, radius).is_contained(ll,ur);
}


std::ostream& operator<<(std::ostream& os, const nBall& b) {
  return os << "nBall(" << b.center << "," << b.radius << "," << b.one << ")";
}

std::ostream& operator<<(std::ostream& os, const nBall_stuff& b) {
  return os << "(" << b.contained << "," << b.last_gen << "," << b.depth << "," << b.ball << ")";
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
  
  //std::cout << "Finding min initial radius for nifs with param " << c
  //          << " and centers: \n";
  //for (int i=0; i<(int)centers.size(); ++i) {
  //  std::cout << centers[i] << "\n";
  //}
  //std::cout << "(1-|c|) = " << fabs(1.0-abs(c)) << "\n";
  
  if (abs(c) > 0.99) return -1;
  for (int i=0; i<(int)centers.size(); ++i) {
    double lbi = abs(centers[i]*(1.0-c))/(1.0-abs(c));
    if (lb < 0 || lbi > lb) {
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
  cpx sub_center = b.center + b.one * (1.0-c)*centers[i];
  //std::cout << "Mapping " << b << " to " << nBall( sub_center, abs(c*b.radius), c*b.one) << "\n";
  return nBall( sub_center, abs(c*b.radius), c*b.one);
}



bool nIFS::is_connected(int depth, int& difficulty) {
  double min_r = minimal_initial_radius();
  if (min_r < 0) return true;
  
  nBall initial_ball(0,min_r,1);
  std::vector<std::pair<nBall_stuff,nBall_stuff> > stack(0); 
  stack.push_back( std::make_pair( nBall_stuff(false, 0, 1, act_on_right(0, initial_ball)),
                                   nBall_stuff(false, 1, 1, act_on_right(1, initial_ball)) ));
  difficulty = 0;
  while (stack.size() > 0) {
    
    //std::cout << "Stack: " << "\n";
    //for (int i=0; i<(int)stack.size(); ++i) {
    //  std::cout << stack[i].first << "  -  " << stack[i].second << "\n";
    //}
  
    std::pair<nBall_stuff,nBall_stuff> bp = stack.back();
    stack.pop_back();
    if (bp.first.ball.is_disjoint(bp.second.ball)) {
      continue;
    } else if (bp.first.depth >= depth) {
      //std::cout << "Yep connected\n";
      return true;
    }
    for (int i=0; i<(int)centers.size(); ++i) {
      nBall new_b1 = act_on_right(i, bp.first.ball);
      for (int j=0; j<(int)centers.size(); ++j) {
        nBall new_b2 = act_on_right(j, bp.second.ball);
        stack.push_back(std::make_pair(nBall_stuff(false, 0, bp.first.depth+1, new_b1),
                                       nBall_stuff(false, 1, bp.second.depth+1, new_b2)) );
      }
    }
    ++difficulty;
  }
  //std::cout << "Nope not\n";
  return false;
}








