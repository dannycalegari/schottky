gBall::gBall() {
  center = 0;
  radius = 0;
  ifs_centers = std::vector<cpx>();
}

gBall::gBall(cpx c, double r, const std::vector<cpx>& i_c) {
  center = c;
  radius = r;
  ifs_centers = i_c;
}


bool gBall::is_disjoint(const gBall& other) const {
  return abs( center - other.center ) > radius + other.radius;
}

bool gBall::is_disjoint(const cpx& ll, const cpx& ur) const {
  return Ball(center, radius).is_disjoint(ll, ur);
}

bool gBall::is_contained(const cpx& ll, const cpx& ur)  const {
  return Ball(center, radius).is_contained(ll,ur);
}


std::ostream& operator<<(std::ostream& os, const gBall& b) {
  return os << "gBall(" << b.center << "," << b.radius  << ")";
}

std::ostream& operator<<(std::ostream& os, const gBall_stuff& b) {
  return os << "(" << b.contained << "," << b.last_gen << "," << b.depth << "," << b.ball << ")";
}





gIFS::gIFS() {
  centers = std::vector<cpx>();
  factors = std::vector<cpx>();
}

gIFS::gIFS(const std::vector<cpx>& facs, const std::vector<cpx>& cents) {
  centers = cents;
  factors = facs;
}
  
void gIFS::set_params(const std::vector<cpx>& new_factors) {
  factors = new_factors;
}

void gIFS::set_centers(const std::vector<cpx>& new_centers) {
  centers = new_centers;
}
  
bool gIFS::minimal_ball(gBall& b) {
  
  //average the centers to find where the ball will be centered
  cpx av = 0;
  for (int i=0; i<(int)centers.size(); ++i) {
    av += centers[i];
  }
  cpx ball_center = av / (cpx)centers.size();
  
  //now figure out how big the radius needs to be
  double lb = -1;
  for (int i=0; i<(int)factors.size(); ++i) {
    if (abs(factors[i]) > 0.99) return false;
    cpx ball_center_im = factors[i]*(ball_center - centers[i]) + centers[i];
    double lbi = abs( ball_center - ball_center_im ) / (1-abs(factors[i]));
    if (lb < 0 || lbi > lb) {
      lb = lbi;
    }
  }
  double ball_radius = lb;
  std::vector<cpx> ifs_center_displacements(centers.size());
  for (int i=0; i<(int)centers.size(); ++i) {
    ifs_center_displacements[i] = (factors[i]*(ball_center - centers[i]) + centers[i])  - ball_center;
  }
  
  b = gBall(ball_center, ball_radius, ifs_center_displacements);
  
  return true;
}
  
  
gBall gIFS::act_on_left(int i, const gBall& b) {
  cpx new_center = factors[i]*(b.center - centers[i]) + centers[i];
  double new_radius = abs( factors[i]*b.radius );
  std::vector<cpx> new_displacements(centers.size());
  for (int j=0; j<(int)centers.size(); ++j) {
    new_displacements[j] = factors[i]*b.ifs_centers[j];
  }
  return gBall( new_center, new_radius, new_displacements );
}

gBall gIFS::act_on_right(int i, const gBall& b) {
  cpx new_center = b.center + b.ifs_centers[i];
  double new_radius = abs(factors[i]*b.radius);
  std::vector<cpx> new_displacements(centers.size());
  for (int j=0; j<(int)centers.size(); ++j) {
    new_displacements[j] = factors[i]*b.ifs_centers[j];
  }
  return gBall( new_center, new_radius, new_displacements );
}



bool gIFS::is_connected(int depth, int& difficulty) {
  gBall initial_ball;
  if (!minimal_ball(initial_ball)) return true;
  
  std::vector<std::pair<gBall_stuff,gBall_stuff> > stack(0); 
  stack.push_back( std::make_pair( gBall_stuff(false, 0, 1, act_on_right(0, initial_ball)),
                                   gBall_stuff(false, centers.size()-1, 1, act_on_right(centers.size()-1, initial_ball)) ));
  
  //for (int i=0; i<(int)centers.size(); ++i) {
  //  for (int j=i+1; j<(int)centers.size(); ++j) {
  //    stack.push_back( std::make_pair( gBall_stuff(false, i, 1, act_on_right(i, initial_ball)),
  //                                     gBall_stuff(false, j, 1, act_on_right(j, initial_ball)) ));
  //  }
  //}
  difficulty = 0;
  while (stack.size() > 0) {
    
    //std::cout << "Stack: " << "\n";
    //for (int i=0; i<(int)stack.size(); ++i) {
    //  std::cout << stack[i].first << "  -  " << stack[i].second << "\n";
    //}
  
    std::pair<gBall_stuff,gBall_stuff> bp = stack.back();
    stack.pop_back();
    if (bp.first.ball.is_disjoint(bp.second.ball)) {
      continue;
    } else if (bp.first.depth >= depth) {
      //std::cout << "Yep connected\n";
      return true;
    }
    for (int i=0; i<(int)centers.size(); ++i) {
      gBall new_b1 = act_on_right(i, bp.first.ball);
      for (int j=0; j<(int)centers.size(); ++j) {
        gBall new_b2 = act_on_right(j, bp.second.ball);
        stack.push_back(std::make_pair(gBall_stuff(false, bp.first.last_gen, bp.first.depth+1, new_b1),
                                       gBall_stuff(false, bp.second.last_gen, bp.second.depth+1, new_b2)) );
      }
    }
    ++difficulty;
  }
  //std::cout << "Nope not\n";
  return false;
}




