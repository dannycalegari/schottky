#include <algorithm>

double delta(double C, double Cp, double epsilon) {
  double s = sqrt( ( (C*C + C*Cp + Cp*Cp)*epsilon*epsilon) / pow(Cp,4) );
  double d = C*epsilon + Cp*(epsilon + Cp*s);
  return Cp*epsilon / d;
}

//find the maximum and minimum possible values of the derivative
//around a ball
void ifs::deriv_bounds_around_ball(const Bitword& u, 
                                    cpx z0, 
                                    double r, 
                                    double& max,
                                    double& min) {
  int N = 100;
  cpx current_z;
  double E = 2.718281828459045;
  double PI = 3.1415926535897932;
  double small_r = PI*r / (double)N;
  cpx angle;
  max = -1;
  min = -1;
  cpx current_deriv;
  double err;
  double this_max, this_min;
  double deriv2_bound; 
  //std::cout << "Finding derivative bounds on the disk: " << z0 << "," << r << " (small_r: " << small_r << "\n";
  for (int i=0; i<N; ++i) {
    angle = cpx( cos(i*2*PI/(double)N), sin(i*2*PI/(double)N) );
    current_z = z0 + r*angle;
    word_deriv(u, current_z, current_deriv, err);
    deriv2_bound = 2/pow(1-(abs(current_z)+small_r), 3);
    //std::cout << "i=" << i << " angle=" << angle << " current_z=" << current_z
    //          << " current_deriv=" << current_deriv << " deriv2_bound=" << deriv2_bound << "\n";
    this_max = abs(current_deriv) + err + small_r*deriv2_bound;
    if (this_max > max || max < 0) max = this_max;
    this_min = abs(current_deriv) - err - small_r*deriv2_bound;
    if (this_min < min || min < 0) min = this_min;
    //std::cout << "this_max=" << this_max << " this_min=" << this_min << "\n";
  }
}



bool ifs::certify_set_B_point(const Bitword& u, bool certify_all, double& within) {
  //find out how far z is from taking u^inf to 1/2
  cpx hz0 = apply_bitword(u, 0.5);
  double epsilon = abs(hz0 - 0.5);
  double ball_rad = abs(0.5*z-0.5)/(1.0-az);
  if (certify_all) epsilon += ball_rad*pow(az,u.len);
  
  //get the deriv at this point
  cpx deriv;
  double err;
  word_deriv(u, z, deriv, err);
  
  //a first attempt at r will be something times epsilon/deriv
  double initial_r = (epsilon / (abs(deriv) - err));
  double r = initial_r;
  double C,Cp,d,RHS;
  
  //std::cout << "Certifying set B within r=" << r << "\n";
  //if (certify_all) {
  //  std::cout << "epsilon=" << abs(hz0 - 0.5) << "+" << ball_rad*pow(az,u.len) << " = " << epsilon << "\n";
  //} else {
  //  std::cout << "epsilon=" << epsilon << "\n";
  //}
  
  do {
    r *= 2;
    deriv_bounds_around_ball(u, z, r, C, Cp);
    d = delta(C, Cp, epsilon);
    RHS = epsilon*(1.0-d+d*d)/(d*((1.0-d)*Cp - d*C));
    //std::cout << "r=" << r << " RHS=" << RHS << " C=" << C << " Cp=" << Cp << " d=" << d << "\n";
  } while (r < RHS && r < 257*initial_r && C > 0 && Cp > 0);
  if (r > 257*initial_r || C < 0 || Cp < 0) {
    //std::cout << "No good\n";
    return false;
  }
  
  //std::cout << "Good\n";
  
  within = epsilon*(1.0-d+d*d)/((1.0-d)*Cp - d*C);
  return true;
}


//returns whether a<=c<=d<=b cyclically
bool contains_cyclic_range(int a, int b, int c, int d, int L) {
  int A=a;
  int B=b;
  int C=c;
  int D=d;
  if (B<A) B += L;
  if (C<A) C += L;
  if (D<A) D += L;
  return (A <= C && C <= D && D <= B);
}
  


std::vector<Bitword> ifs::get_half_balls_along_path(const std::vector<cpx>& path,
                                                    int d,
                                                    int verbose) {
  cpx old_z = z;
  cpx old_w = w;
  std::vector<std::pair<cpx, std::set<Bitword> > > ball_list(path.size());
  std::vector<Bitword> current_balls;
  for (int i=0; i<(int)path.size(); ++i) {
    set_params(path[i], path[i]);
    half_balls(current_balls, d, 0);
    ball_list[i] = std::make_pair(path[i], 
                                  std::set<Bitword>(current_balls.begin(), 
                                                    current_balls.end()));
  }
  
  if (verbose>0) {
    std::cout << "Got initial ball list along path:\n";
    for (int i=0; i<(int)ball_list.size(); ++i) {
      std::cout << i << " " << ball_list[i].first << ":\n";
      std::set<Bitword>::iterator it = ball_list[i].second.begin();
      while (it != ball_list[i].second.end()) {
        std::cout << *it << "\n";
        it++;
      }
    }
  }
  
  int i=0;
  //subdivide the path until every step is a strict subset or superset
  while (i < ball_list.size()) {
    int ip1 = (i+1)%ball_list.size();
    //check if subdivision is necessary
    if (std::includes(ball_list[i].second.begin(), ball_list[i].second.end(),
                      ball_list[ip1].second.begin(), ball_list[ip1].second.end()) 
        ||
        std::includes(ball_list[ip1].second.begin(), ball_list[ip1].second.end(),
                      ball_list[i].second.begin(), ball_list[i].second.end())) {
      ++i;
      continue;
    }
    //subdivide
    cpx new_c = 0.5*(ball_list[i].first + ball_list[ip1].first);
    set_params(new_c,new_c);
    half_balls(current_balls, d, 0);
    ball_list.insert((ip1 == 0 ? ball_list.end() : ball_list.begin()+ip1),
                     std::make_pair(new_c, std::set<Bitword>(current_balls.begin(), 
                                                             current_balls.end())));
  }

  if (verbose>0) {
    std::cout << "After subdividing:\n";
    for (int i=0; i<(int)ball_list.size(); ++i) {
      std::cout << i << " " << ball_list[i].first << ":\n";
      std::set<Bitword>::iterator it = ball_list[i].second.begin();
      while (it != ball_list[i].second.end()) {
        std::cout << *it << "\n";
        it++;
      }
    }
  }

  //this list records the indices on which the balls starts and stops
  std::map<Bitword, std::pair<int,int> > ball_indices;
  for (int i=0; i<ball_list.size(); ++i) {
    int ip1 = (i+1)%ball_list.size();
    std::vector<Bitword> in_first_not_second(ball_list[i].second.size() + 
                                             ball_list[ip1].second.size());
    std::vector<Bitword>::iterator diff_it;
    diff_it = std::set_difference(ball_list[i].second.begin(), ball_list[i].second.end(), 
                                  ball_list[ip1].second.begin(), ball_list[ip1].second.end(), 
                                  in_first_not_second.begin());
    in_first_not_second.resize(diff_it - in_first_not_second.begin());
    for (int j=0; j<(int)in_first_not_second.size(); ++j) {
      ball_indices[in_first_not_second[j]].second = i;
    }
    std::vector<Bitword> in_second_not_first(ball_list[i].second.size() + 
                                             ball_list[ip1].second.size());
    diff_it = std::set_difference(ball_list[ip1].second.begin(), ball_list[ip1].second.end(), 
                                  ball_list[i].second.begin(), ball_list[i].second.end(), 
                                  in_second_not_first.begin());
    in_second_not_first.resize(diff_it - in_second_not_first.begin());
    for (int j=0; j<(int)in_second_not_first.size(); ++j) {
      ball_indices[in_second_not_first[j]].first = ip1;
    }
  }
  
  if (verbose>0) {
    std::cout << "Start/stop indices::\n";
    std::map<Bitword, std::pair<int,int> >::iterator it = ball_indices.begin();
    while (it != ball_indices.end()) {
      std::cout << it->first << ": " << it->second.first << " -> " << it->second.second << "\n";
      it++;
    }
  }
  
  //this removes the balls whose index sets are subsets of other balls
  //this is naive, but hopefully it won't take very long
  std::map<Bitword, std::pair<int,int> >::iterator it = ball_indices.begin();
  std::map<Bitword, std::pair<int,int> >::iterator it2;
  std::map<Bitword, std::pair<int,int> >::iterator it3;
  while (false) { //it != ball_indices.end()) {
    it2 = ball_indices.begin();
    while (it2 != ball_indices.end()) {
      if (it2 == it) {
        it2++;
        continue;
      }
      it3 = it2;
      it3++;
      if (contains_cyclic_range(it->second.first, it->second.second,
                                it2->second.first, it2->second.second,
                                ball_list.size())) {
        if (verbose>0) {
          std::cout << it2->first << ": " << it2->second.first << "->" << it2->second.second << 
                       " is made redundant by " << it->first << ": " << it->second.first << "->" << it->second.second << "\n";
        }
        ball_indices.erase(it2);
      }
      it2 = it3;
    }
    it++;
  }
  //read out the balls that remain, in order of when they start
  std::vector<Bitword> balls(0);
  for (int i=0; i<(int)ball_list.size(); ++i) {
    std::set<Bitword>::iterator it4 = ball_list[i].second.begin();
    std::map<Bitword, std::pair<int,int> >::iterator it5;
    while (it4 != ball_list[i].second.end()) {
      it5 = ball_indices.find( *it4 );
      if (it5 != ball_indices.end() && it5->second.first == i) {
        balls.push_back( *it4 );
      }
      it4++;
    }
  }
  set_params(old_z, old_w);
  return balls; 
}


//find a parameter taking the infinite word u^inf to 1/2
cpx ifs::solve_for_half(const Bitword& u, cpx start, double tol) {
  cpx current_z = start;
  Bitword upow = u.pow(64/u.len);
  cpx current_output = upow.apply(current_z, 0.5);
  cpx deriv, step;
  double err;
  int i=0;
  while (abs(current_output - 0.5) > tol) {
    word_deriv(upow, current_z, deriv, err);
    step = (0.5-current_output)/deriv;
    current_z += step;
    current_output =  upow.apply(current_z, 0.5);
    //++i;
    //if (i%1000 == 0 && i>0) {
    //  std::cout << "Took " << i << " steps\n";
    //}
  }
  return current_z;
}


//open a separate window with a picture of the set B contained inside the 
//balls
void ifs::draw_set_B_balls(const std::vector<Bitword>& balls, 
                           cpx initial_point, 
                           int d, 
                           int verbose) {
  cpx old_z = z;
  cpx old_w = w;
  
  //first, get a guess about z for each of the balls
  double approx_radius = abs(0.5*initial_point-0.5)/(1.0-abs(initial_point));
  std::vector<cpx> ball_zs(balls.size());
  std::vector<double> ball_rads(balls.size());
  for (int i=0; i<(int)balls.size(); ++i) {
    ball_zs[i] = solve_for_half(balls[i], 
                                initial_point, 
                                0.01*approx_radius*pow(abs(initial_point), balls[i].len));
    if (verbose>0) {
      std::cout << "Placed ball " << balls[i] << " at " << ball_zs[i] << "\n";
    }
  }
  //certify all of the balls within some radius
  for (int i=0; i<(int)balls.size(); ++i) {
    set_params(ball_zs[i], ball_zs[i]);
    if (!certify_set_B_point(balls[i], true, ball_rads[i])) {
      std::cout << "Couldn't certify a ball; aborting\n";
      return;
    }
    if (verbose>0) {
      std::cout << "Certified " << balls[i] << " within " << ball_rads[i] << "\n";
    }
  }
  
  //find the ball extents and everything
  cpx ll = ball_zs[0] - cpx(ball_rads[0], ball_rads[0]);
  cpx ur = ball_zs[0] + cpx(ball_rads[0], ball_rads[0]);
  for (int i=1; i<(int)balls.size(); ++i) {
    double a = ball_zs[i].real() - ball_rads[i];
    double b = ball_zs[i].real() + ball_rads[i];
    double c = ball_zs[i].imag() - ball_rads[i];
    double d = ball_zs[i].imag() + ball_rads[i];
    if (a < ll.real()) ll = cpx(a, ll.imag());
    if (c < ll.imag()) ll = cpx(ll.real(), c);
    if (b > ur.real()) ur = cpx(b, ur.imag());
    if (d > ur.imag()) ur = cpx(ur.real(), d);
  }
  cpx center = 0.5*(ur+ll);
  double height_res = (ur-center).imag();
  double width_res = (ur-center).real();
  double box_rad = (height_res > width_res ? height_res : width_res);
  ll = center - cpx(box_rad,box_rad);
  ur = center + cpx(box_rad,box_rad);
  
  //get the pixel data
  int num_pix = 800;
  double pixel_width = 2.0*box_rad/(double)num_pix;
  
  if (verbose>0) {
    std::cout << "ll: " << ll << "\nur: " << ur << "\n";
    std::cout << "Pixel width: " << pixel_width << "\n";
  }
  
  //draw the circles
  XGraphics X2(num_pix, num_pix, 1, Point2d<float>(0,0));
  int bcol = X2.get_rgb_color(0,0,1);
  for (int i=0; i<(int)balls.size(); ++i) {
    Point2d<int> c( (ball_zs[i].real() - ll.real())/pixel_width,
                    (ball_zs[i].imag() - ll.imag())/pixel_width );
    double r = ball_rads[i]/pixel_width;
    X2.draw_circle(c, r, bcol);
    if (verbose>0) {
      std::cout << "Drew circle " << c << " radius " << r << "\n";
    }
  }
  
  //now draw the little disks inside
  for (int i=0; i<(int)balls.size(); ++i) {
    int blcol = X2.get_rgb_color(0,
                                 (double)i/(double)balls.size(), 
                                 (double)(balls.size()-i)/(double)balls.size());
    if (d==0) break;
    for (int j=0; j<1<<d; ++j) {
      Bitword u_added = balls[i].append(j, d);
      cpx little_z = solve_for_half(u_added, ball_zs[i], 0.01*approx_radius*pow(abs(ball_zs[i]), u_added.len));
      set_params(little_z, little_z);
      double little_radius;
      if (!certify_set_B_point(u_added, true, little_radius)) {
        std::cout << "Couldn't certify disk\n";
        return;
      }
      if (i==0 && verbose>0) {
        std::cout << "Did " << u_added << " at " << little_z << " radius " << little_radius << "\n";
      }
      Point2d<int> c( (little_z.real() - ll.real())/pixel_width,
                      (little_z.imag() - ll.imag())/pixel_width );
      double r = little_radius/pixel_width;
      X2.draw_disk(c, r, blcol);
    }
  }
  
  
  X2.wait_for_key();
  set_params(old_z, old_w);
}



std::vector<Bitword> ifs::get_certified_half_balls_along_path(const std::vector<cpx>& path,
                                                              int d,
                                                              int verbose) {
  std::vector<Bitword> balls = get_half_balls_along_path(path, d, verbose);
  cpx old_z = z;
  cpx old_w = w;
  cpx initial_point = path[0];
  double approx_radius = abs(0.5*initial_point-0.5)/(1.0-abs(initial_point));
  std::vector<cpx> ball_zs(balls.size());
  std::vector<double> ball_rads(balls.size());
  for (int i=0; i<(int)balls.size(); ++i) {
    ball_zs[i] = solve_for_half(balls[i], 
                                initial_point, 
                                0.005*approx_radius*pow(abs(initial_point), balls[i].len));
    if (verbose>0) {
      std::cout << "Placed ball " << balls[i] << " at " << ball_zs[i] << "\n";
    }
  }
  //certify all of the balls within some radius
  for (int i=0; i<(int)balls.size(); ++i) {
    set_params(ball_zs[i], ball_zs[i]);
    if (!certify_set_B_point(balls[i], true, ball_rads[i])) {
      std::cout << "Couldn't certify a ball; aborting\n";
      return std::vector<Bitword>(0);
    }
    if (verbose>0) {
      std::cout << "Certified " << balls[i] << " within " << ball_rads[i] << "\n";
    }
  }
  set_params(old_z, old_w);
  return balls;
}


bool ifs::certify_set_B_path(const std::vector<cpx>& path, int initial_depth, int verbose) {
  //first, get a list of the balls at that depth
  std::vector<Bitword> initial_balls = get_half_balls_along_path(path, initial_depth, verbose);
  if (verbose>0) {
    std::cout << "Half balls along path:\n";
    for (int i=0; i<(int)initial_balls.size(); ++i) {
      std::cout << i << ": " << initial_balls[i] << "\n";
    }
  }
  
  draw_set_B_balls(initial_balls, path[0], 10, verbose);
  
  
  return true;
}


std::vector<Ball> ifs::subdivide_half_prefix(const Bitword& u, 
                                             cpx start_z,
                                             int d, cpx ll, cpx ur) {
  cpx old_z = z;
  cpx old_w = w;
  double approx_radius = abs(0.5*start_z-0.5)/(1.0-abs(start_z));
  //set up the initial ball
  std::deque<Ball> stack(1);
  stack[0].word = u.w;
  stack[0].word_len = u.len;
  stack[0].center = solve_for_half(u, start_z, 0.005*approx_radius*pow(abs(start_z), u.len));
  set_params(stack[0].center, stack[0].center);
  if (!certify_set_B_point(u, true, stack[0].radius)) {
    std::cout << "Couldn't certify disk\n";
    return std::vector<Ball>(0);
  }
  
  std::vector<Ball> ans(0);
  
  if (d == 0) {
    ans.resize(1);
    ans[0] = stack[0];
    return ans;
  }
  
  double max_radius = 0;
  std::cout << "Starting to subdivide half prefix " << u << "\n";
  
  while (stack.size() > 0) {
    Ball b = stack.back();
    stack.pop_back();
    
    //std::cout << "Word: " << b.word << " word len: " << b.word_len << " center: " << b.center << " radius: " << b.radius << "\n";
    
    //find the next step balls
    Ball b0(b.center, 0, 0, b.radius, b.word << 1, b.word_len+1);
    approx_radius = abs(0.5*b0.center-0.5)/(1.0-abs(b0.center));
    Bitword b0word = Bitword(b0.word, b0.word_len);
    //std::cout << "Solving for new center within radius " << 0.005*approx_radius*pow(abs(b0.center), b0.word_len) << "\n";
    b0.center = solve_for_half(b0word, 
                               b0.center, 
                               0.005*approx_radius*pow(abs(b0.center), b0.word_len));
    //std::cout << "New center: " << b0.center << "\n";
    set_params(b0.center, b0.center);
    if (!certify_set_B_point(b0word, true, b0.radius)) {
      std::cout << "Couldn't certify disk\n";
      return std::vector<Ball>(0);
    }
    //std::cout << "New radius: " << b0.radius << "\n";
    
    Ball b1(b.center, 0, 0, b.radius, (b.word << 1).flip(0), b.word_len+1);
    approx_radius = abs(0.5*b1.center-0.5)/(1.0-abs(b1.center));
    Bitword b1word = Bitword(b1.word, b1.word_len);
    b1.center = solve_for_half(b1word,
                               b1.center, 
                               0.005*approx_radius*pow(abs(b1.center), b1.word_len));
    set_params(b1.center, b1.center);
    if (!certify_set_B_point(b1word, true, b1.radius)) {
      std::cout << "Couldn't certify disk\n";
      return std::vector<Ball>(0);
    }
    
    //decide what to do with them
    if (!b0.is_disjoint(ll, ur)) {
      if (b0.word_len - u.len >= d) {
        ans.push_back(b0);
        if (b0.radius > max_radius) max_radius = b0.radius;
      } else {
        stack.push_front(b0);
      }
    }
    if (!b1.is_disjoint(ll, ur)) {
      if (b1.word_len - u.len >= d) {
        ans.push_back(b1);
        if (b1.radius > max_radius) max_radius = b1.radius;
      } else {
        stack.push_front(b1);
      }
    }
  }
  std::cout << "Done subdividing; max radius: " << max_radius << "\n";
  set_params(old_z,old_w);
  return ans;
}












