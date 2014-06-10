#include <deque>

#include "ifs.h"


/**************************************************************************
 * There are a bunch of big functions, so the ifs class is in 
 * multiple files.  This file contains only basic computational functions
 * ************************************************************************/

//these are the other files
#include "ifs_draw.cc"         //everything about drawing it to the xwindow
#include "ifs_interface.cc"    //everything about interacting
#include "ifs_connected.cc"    //the function to detect if the IFS is connected 
#include "ifs_trap.cc"         //functions to build a trap
#include "ifs_trap_like.cc"    //functions to find trap balls via u,v words
#include "ifs_set_A.cc"        //function to find the boundary of holes in set A

//first some ball functions
Ball::Ball() { 
  center = 0.5;
  to_z = to_w = 0;
  radius = 1.0;
  word = std::bitset<64>(0);
  word_len = 0;
}

Ball::Ball(cpx c, double r) {
  center = c;
  radius = r;
  to_z = to_w = 0;
  word = std::bitset<64>(0);
  word_len = 0;
}

Ball::Ball(cpx c, cpx tz, cpx tw, double r) {
  center = c;
  to_z = tz;
  to_w = tw;
  radius = r;
  word = std::bitset<64>(0);
  word_len = 0; //just a ball, no words, to start
}

Ball::Ball(cpx c, cpx tz, cpx tw, double r, const std::bitset<64>& w, int wl) {
  center = c;
  to_z = tz;
  to_w = tw;
  radius = r;
  word = w;
  word_len = wl; //just a ball, no words, to start
}

int Ball::last_gen_index() const {
  return int(word[word_len-1]);
}

bool Ball::is_disjoint(const Ball& b) {
  return abs(center - b.center) > radius + b.radius;
}

bool Ball::is_disjoint(const cpx& ll, const cpx& ur) {
  cpx ul(ll.real(), ur.imag());
  cpx lr(ur.real(), ll.imag());
  cpx& bc = center;
  double& r = radius;
  if (ur.real() < bc.real() && ur.imag() < bc.imag()) { //upper right
    return !(abs(bc-ur)<r);
  } else if (ul.real() > bc.real() && ul.imag() < bc.imag()) {
    return !(abs(bc-ul)<r);
  } else if (ll.real() > bc.real() && ll.imag() > bc.imag()) {
    return !(abs(bc-ll)<r);
  } else if (lr.real() < bc.real() && lr.imag() > bc.imag()) {
    return !(abs(bc-lr)<r);
  } else if (ur.real() < bc.real()) {
    return !(ur.real() > bc.real() - r);
  } else if (ul.imag() < bc.imag()) {
    return !(ul.imag() > bc.imag() -r);
  } else if (ll.real() > bc.real()) {
    return !(ll.real() < bc.real() + r);
  } else if (lr.imag() > bc.imag()) {
    return !(lr.imag() < bc.imag() + r);
  }
  return false; //it's right over the square
}

bool Ball::is_contained(const cpx& ll, const cpx& ur) {
  return (ll.real() < center.real() - radius) && 
         (center.real() + radius < ur.real()) && 
         (ll.imag() < center.imag() - radius) && 
         (center.imag() + radius < ur.imag());
}


std::ostream& operator<<(std::ostream& os, const Ball& b) {
  return os << "Ball(" << b.center << "," << b.to_z << "," << b.to_w << "," << b.radius << "," << b.word << "," << b.word_len << ")";
}




bool Bitword::operator<(const Bitword& b) const {
  return w.to_ulong() < b.w.to_ulong();
}

Bitword Bitword::prefix(int n) const {
  std::string s = w.to_string();
  s = s.substr(64-len, n);
  std::bitset<64> b(s);
  return Bitword(b, n);
}

Bitword Bitword::suffix(int n) const {
  std::string s = w.to_string();
  s = s.substr(64-n, n);
  std::bitset<64> b(s);
  return Bitword(b, n);
}
  
std::string Bitword::str() const {
  return w.to_string().substr(64-len, len);
}

std::ostream& operator<<(std::ostream& os, const Bitword& b) {
  std::string w = b.w.to_string();
  w = w.substr(64-b.len, b.len);
  return os << w;
}





/****************  IFS functions ****************************/
ifs::ifs() {
}

ifs::ifs(cpx a, cpx b, int width, int mode) {
  initialize(a,b, width, mode);
  
}

//set everything and start graphics
void ifs::initialize(cpx a, cpx b, int width, int mode){
  // initialize z to a and w to b
  z=a;
  w=b;
  az = abs(z);
  aw = abs(w);
  sync=0;
  color_ifs=true;
  chunky_ifs=false;
  disconnection_depth=false;
  draw_contains_half = false;
  draw_trap_mode = false;
  find_trap_like_vectors = false;
  find_close_uv_words = false;
  step=0.01;	// size of adjustments to z and w
  seed=0.0;	// initial seed point for IFS
  center=0.0;	// in mandelbrot mode; center of screen, size of window, and mesh of accuracy
  wind=1.0;
  mesh=2;
  depth=10;	  // depth to iterate IFS or detect connectedness to
  trap_depth = depth;  //depth to search for traps 
  drawing_width = width;
  drawing_radius = drawing_width/2;  
  this->mode = mode;
}

//reset stuff when switching
void ifs::reinitialize(cpx a, cpx b) {
  z=a;
  w=b;
  az = abs(z);
  aw = abs(w);
}





//compute the image ball
Ball ifs::act_on_left(int index, const Ball& b) const {
  if (index == 0) {
    return Ball( z*b.center, z*b.to_z, z*b.to_w, az*b.radius, b.word, b.word_len+1 );
  } else {
    std::bitset<64> temp(0);
    temp[b.word_len] = 1;
    return Ball( (w*(b.center - 1.0)) + 1.0, w*b.to_z, w*b.to_w, aw*b.radius, b.word | temp, b.word_len+1 );
  }
}


Ball ifs::act_on_right(int index, const Ball& b) const {
  if (index == 0) {
    return Ball( b.center + b.to_z, z*b.to_z, z*b.to_w, az*b.radius, b.word<<1, b.word_len+1 );
  } else {
    return Ball( b.center + b.to_w, w*b.to_z, w*b.to_w, aw*b.radius, (b.word<<1)| std::bitset<64>(1), b.word_len+1 );
  }
}


//compute the image of a point under the action of a word
cpx ifs::apply_bitword(const Bitword& b, cpx x) const {
  cpx ans = x;
  for (int i=0; i<(int)b.len; ++i) {
    if (b.w[i] == 0) {
      ans = z*ans;
    } else {
      ans = w*(ans-1.0)+1.0;
    }
  }
  return ans;
}


double ifs::distance_from_balls(cpx p, const std::vector<Ball>& balls) {
  if (balls.size() == 0) return -1;
  double dist = abs(balls[0].center-p)-balls[0].radius;
  for (int i=1; i<(int)balls.size(); ++i) {
    double d = abs(balls[i].center-p) - balls[i].radius;
    if (d < dist) dist = d;
  }
  return dist;
}


double ifs::when_ray_hits_ball(cpx p, cpx v, const Ball& b) {
  cpx c = b.center;
  double r = b.radius;
  double t = cpx_dot(c-p, v) / cpx_dot(v,v);
  cpx x = p + t*v;
  double a1 = abs(c-x);
  if (a1 > r) return 1e12;
  double a2 = sqrt(r*r-a1*a1);
  return t - a2/abs(v);
}


double ifs::when_ray_hits_ball(cpx p, cpx v, const std::vector<Ball>& balls) {
  if (balls.size() == 0) return -1;
  double t = when_ray_hits_ball(p,v,balls[0]);
  for (int i=1; i<(int)balls.size(); ++i) {
    double tp = when_ray_hits_ball(p,v,balls[i]);
    if (tp < t) t = tp;
  }
  return t;
}



//take a list of balls and compute one more level of depth 
//(on the left)
void ifs::compute_next_ball_depth(std::vector<Ball>& balls, int current_depth) {
  std::vector<Ball> balls_temp;
  balls_temp.swap(balls);
  int L_cur = 1<<current_depth;
  int L_new = 1<<(current_depth+1);
  balls.resize(L_new);
  for (int j=0; j<L_cur; ++j) {
    //for each j, we want to append on the left both a 0 and a 1
    //(the highest bit position is the last function we applied)
    Ball parent_ball = balls_temp[j];
    balls[j] = act_on_left(0, parent_ball);
    balls[(j | (1<<current_depth))] = act_on_left(1, parent_ball);
  }
}

//take a seed ball and compute all the image balls
//create a list of image points
//each image point is indexed by the binary digits
//so 1011 means fgff, and it's a left action
//of course we need to know the word length to parse how many g's are in front
void ifs::compute_balls(std::vector<Ball>& balls, const Ball& ball_seed, int compute_depth) {
  std::vector<Ball> balls_temp;
  balls.resize(2);
  balls[0] = act_on_left(0, ball_seed);
  balls[1] = act_on_left(1, ball_seed);
  for (int i=1; i<compute_depth; ++i) {
    compute_next_ball_depth(balls, i);
  }
}


//these functions are the same as above, except they act on the right
void ifs::compute_next_ball_depth_right(std::vector<Ball>& balls, int current_depth) {
  std::vector<Ball> balls_temp;
  balls_temp.swap(balls);
  int L_cur = 1<<current_depth;
  int L_new = 1<<(current_depth+1);
  balls.resize(L_new);
  for (int j=0; j<L_cur; ++j) {
    //for each j, we want to append on the right both a 0 and a 1
    //(the lowest bit position is the last function we applied)
    Ball parent_ball = balls_temp[j];
    balls[j << 1] = act_on_right(0, parent_ball);
    balls[(j << 1) | 1] = act_on_right(1, parent_ball);
  }
}
void ifs::compute_balls_right(std::vector<Ball>& balls, const Ball& ball_seed, int compute_depth) {
  balls.resize(2);
  balls[0] = act_on_right(0, ball_seed);
  balls[1] = act_on_right(1, ball_seed);
  for (int i=1; i<compute_depth; ++i) {
    compute_next_ball_depth_right(balls, i);
  }
}







//given a ball, check if it is disjoint from the rectangle given
//this is lame as always with these functions
bool ifs::is_ball_disjoint(const Ball& b, const cpx& ll, const cpx& ur) {
  const cpx& bc = b.center;
  const double& r = b.radius;
  cpx ul(ll.real(), ur.imag());
  cpx lr(ur.real(), ll.imag());
  if (ur.real() < bc.real() && ur.imag() < bc.imag()) { //upper right
    return !(abs(bc-ur)<r);
  } else if (ul.real() > bc.real() && ul.imag() < bc.imag()) {
    return !(abs(bc-ul)<r);
  } else if (ll.real() > bc.real() && ll.imag() > bc.imag()) {
    return !(abs(bc-ll)<r);
  } else if (lr.real() < bc.real() && lr.imag() > bc.imag()) {
    return !(abs(bc-lr)<r);
  } else if (ur.real() < bc.real()) {
    return !(ur.real() > bc.real() - r);
  } else if (ul.imag() < bc.imag()) {
    return !(ul.imag() > bc.imag() -r);
  } else if (ll.real() > bc.real()) {
    return !(ll.real() < bc.real() + r);
  } else if (lr.imag() > bc.imag()) {
    return !(lr.imag() < bc.imag() + r);
  }
  return false; //it's right over the square
}


//return the lower left and upper right of a box (square) containing
//all the balls
void ifs::box_containing_balls(const std::vector<Ball>& balls, 
                               cpx& ll, 
                               cpx& ur) {
  ll = cpx(balls[0].center.real() - balls[0].radius, balls[0].center.imag() - balls[0].radius);
  ur = cpx(balls[0].center.real() + balls[0].radius, balls[0].center.imag() + balls[0].radius);
  for (int i=1; i<(int)balls.size(); ++i) {
    const cpx c = balls[i].center;
    const double r = balls[i].radius;
    if (c.real() + r > ur.real()) {
      ur = cpx(c.real()+r, ur.imag());
    }
    if (c.imag() + r > ur.imag()) {
      ur = cpx(ur.real(), c.imag() + r);
    }
    if (c.real() - r < ll.real()) {
      ll = cpx(c.real()-r, ll.imag());
    }
    if (c.imag() -r < ll.imag()) {
      ll = cpx(ll.real(), c.imag()-r);
    }
  }
  //std::cout << ll << " " << ur << "\n";
  double pwr = (ur.real() - ll.real())/2.0;
  double phr = (ur.imag() - ll.imag())/2.0;
  cpx center(ll.real() + pwr, ll.imag() + phr);
  double rad = (pwr > phr ? pwr : phr);
  ll = cpx(center.real() - rad, center.imag() - rad);
  ur = cpx(center.real() + rad, center.imag() + rad);  
}


void box_containing_points(const std::vector<cpx>& points, cpx& ll, cpx& ur) {
  ll = cpx(points[0].real(), points[0].imag());
  ur = cpx(points[0].real(), points[0].imag());
  for (int i=1; i<(int)points.size(); ++i) {
    const cpx c = points[i];
    if (c.real() > ur.real()) {
      ur = cpx(c.real(), ur.imag());
    }
    if (c.imag() > ur.imag()) {
      ur = cpx(ur.real(), c.imag());
    }
    if (c.real() < ll.real()) {
      ll = cpx(c.real(), ll.imag());
    }
    if (c.imag() < ll.imag()) {
      ll = cpx(ll.real(), c.imag());
    }
  }
  //std::cout << ll << " " << ur << "\n";
  double pwr = (ur.real() - ll.real())/2.0;
  double phr = (ur.imag() - ll.imag())/2.0;
  cpx center(ll.real() + pwr, ll.imag() + phr);
  double rad = (pwr > phr ? pwr : phr);
  ll = cpx(center.real() - rad, center.imag() - rad);
  ur = cpx(center.real() + rad, center.imag() + rad);  
}





//for all the balls, if they are outside the box, forget it
//if not, hit them on the right with z and w, and keep any 
//disks we touch the given box
void ifs::refine_balls_into_box(std::vector<Ball>& balls, 
                                const cpx& ll, 
                                const cpx& ur) {
  std::vector<Ball> bb;
  bb.swap(balls);
  balls.resize(0);
  int nb = (int)bb.size();
  for (int i=0; i<nb; ++i) {
    if (is_ball_disjoint(bb[i], ll, ur)) continue;
    Ball b1 = act_on_right(0,bb[i]);
    Ball b2 = act_on_right(1,bb[i]);
    if (!is_ball_disjoint(b1, ll, ur)) {
      balls.push_back(b1);
    }
    if (!is_ball_disjoint(b2, ll, ur)) {
      balls.push_back(b2);
    }
  }
}

//find two actions of the given length (words u and v) whose images of 
//p are as close as possible and which begin with different letters
void ifs::find_close_images_with_distinct_first_letters(const Ball& b, 
                                                        int length, 
                                                        Ball& zb, Ball& wb){
  zb = act_on_left(0,b);
  wb = act_on_left(1,b);
  //at every stage, replace b0 and b1 with their images which are closest
  for (int i=1; i<length; ++i) {
    Ball zbb[2] = { act_on_right(0, zb),
                    act_on_right(1, zb) };
    Ball wbb[2] = { act_on_right(0, wb),
                    act_on_right(1, wb) };
    //need to compare all pairs b0* - b1*
    double d[4] = { abs( zbb[0].center - wbb[0].center ),
                    abs( zbb[0].center - wbb[1].center ),
                    abs( zbb[1].center - wbb[0].center ),
                    abs( zbb[1].center - wbb[1].center ) };
    int min_i = 0;
    for (int j=1; j<4; ++j) {
      if (d[j] < d[min_i]) {
        min_i = j;
      }
    }
    zb = zbb[ (min_i>>1)&1 ];
    wb = wbb[ min_i&1 ];
  }
  
}


//find balls (i.e. words) such that the points p1 and p2 are close 
//when acted upon.  The ball is around so that we can act on the right
//if the ratio goes below the ratio goal, we just return it
void ifs::find_aligned_images_with_distinct_first_letters(const Ball& initial_ball, 
                                                          cpx p1, cpx p2, int search_depth,
                                                          Ball& zb, Ball& wb, 
                                                          double ratio_goal, double ratio_lower_limit) {
  //figure out how to get to p1 and p2 using the vector to_z
  //p1 = 0.5 + c1*to_z, so we can always get to the image of 
  //p1 by doing center + c1*to_z
  cpx c1 = (p1-initial_ball.center)/initial_ball.to_z;
  cpx c2 = (p2-initial_ball.center)/initial_ball.to_z;
  std::deque<std::pair<Ball, Ball> > stack(1);
  stack[0] = std::make_pair( act_on_left(0, initial_ball), 
                             act_on_left(1, initial_ball) );
  std::pair<Ball,Ball> best_pair = stack[0];
  double best_dist = abs((stack[0].first.center + c1*stack[0].first.to_z) -
                         (stack[0].second.center + c2*stack[0].second.to_z));
  best_dist /= stack[0].first.radius;
  while (stack.size() > 0) {
    std::pair<Ball,Ball> b = stack.back();
    stack.pop_back();
    
    //see whether the pair is good...
    double new_dist = abs((b.first.center + c1*b.first.to_z) -
                          (b.second.center + c2*b.second.to_z));
    new_dist /= b.first.radius;
    if (new_dist < best_dist && new_dist > ratio_lower_limit) {
      //std::cout << "Found better pair: " << b.first << " " << b.second << " p dist: " << new_dist << "\n";
      //std::cout << "{" << b.first.word_len << "," << new_dist << "}, ";
      if (new_dist < ratio_goal) {
        zb = b.first;
        wb = b.second;
        return;
      }
      best_dist = new_dist;
      best_pair = b;
    }
    
    //go down the stack a little
    if (b.first.word_len - initial_ball.word_len >= search_depth) continue;
    
    Ball zbb[2] = { act_on_right(0, b.first),
                    act_on_right(1, b.first) };
    Ball wbb[2] = { act_on_right(0, b.second),
                    act_on_right(1, b.second) };
    
    for (int i=0; i<4; ++i) {
      if (!zbb[i>>1].is_disjoint(wbb[i&1])) {
        stack.push_front(std::make_pair(zbb[i>>1], wbb[i&1]));
      }
    }

  }
  zb = best_pair.first;
  wb = best_pair.second;
}


//returns the u,v pair which are closest of the given length
//
//it does a breadth-first search, and keeps track of the closest
//pair it finds in each length.  If a pair under consideration
//move necessarily result in a farther pair, it can be discarded
void ifs::find_closest_uv_words(std::vector<std::pair<Bitword,Bitword> >& words, 
                                int uv_depth,
                                double last_step_tolerance,
                                int list_size_max) {
  words.resize(0);
  std::vector<double> distances(0);
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return;
  
  Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
  std::vector<std::pair<Ball, Ball> > pairs(1);
  std::vector<std::pair<Ball, Ball> > next_pairs;
  pairs[0] = std::make_pair(act_on_right(0,b), act_on_right(1,b));

  for (int i=2; i<=uv_depth; ++i) {
    //make all possible next pairs
    next_pairs.resize(4*pairs.size());
    for (int j=0; j<(int)pairs.size(); ++j) {
      Ball& bz = pairs[j].first;
      Ball& bw = pairs[j].second;
      Ball bzs[2] = {act_on_right(0, bz), act_on_right(1, bz)};
      Ball bws[2] = {act_on_right(0, bw), act_on_right(1, bw)};
      for (int k=0; k<4; ++k) {
        next_pairs[4*j + k] = std::make_pair( bzs[k>>1], bws[k&1] );
      }
    }
    //go through and find the closest
    double min_d = 1e10;
    for (int j=0; j<(int)next_pairs.size(); ++j) {
      double d = abs(next_pairs[j].first.center - next_pairs[j].second.center);
      if (d < min_d) min_d = d;
    }
    //make up the next pairs
    //we want all possibilities, so we allow duplicate lengths
    //(for safety, take a little range)
    double cutoff_dist = (i==uv_depth ? min_d+last_step_tolerance : min_d + pow(az, i)*2*min_r);
    //std::cout << "Cutoff for i=" << i << ": " << cutoff_dist << "\n";
    pairs.resize(0);
    for (int j=0; j<(int)next_pairs.size(); ++j) {
      double d = abs(next_pairs[j].first.center - next_pairs[j].second.center);
      if (d <= cutoff_dist && (list_size_max < 0 ||(int)pairs.size() < list_size_max)) {
        pairs.push_back(next_pairs[j]);
      }
    }  
  }
  //now we've got a bunch of pairs
  words.resize(pairs.size());
  for (int i=0; i<(int)pairs.size(); ++i) {
    words[i] = std::make_pair( Bitword(pairs[i].first.word, pairs[i].first.word_len),
                               Bitword(pairs[i].second.word, pairs[i].second.word_len) );
    //std::cout << abs(pairs[i].first.center - pairs[i].second.center) << " " << pairs[i].first.center << " - ";
  }
  //std::cout << "\n";
}


void ifs::find_closest_uv_words_along_path(const std::vector<cpx>& path, 
                                           bool closed_path, 
                                           int word_len) {
  std::vector<std::pair<Bitword,Bitword> > words;
  for (int i=0; i<(int)path.size(); ++i) {
    this->set_params(path[i], path[i]);
    find_closest_uv_words(words, word_len);
    for (int j=0; j<(int)words.size(); ++j) {
      std::cout << "(" << words[j].first << " " << words[j].second << ") ";
    }
    std::cout << "\n";
  }
}
  


//this tries to find 
//(1) a pair of balls of depth n such that |u1/2-v1/2| < epsilon
//(2) that every other pair of balls of depth n is (recursively)
//at least 12*epsilon far away
bool ifs::close_to_set_C(int n_depth, double epsilon) {
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return false;
  
  Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
  std::vector<std::pair<Ball, Ball> > pairs(1);
  std::vector<std::pair<Ball, Ball> > next_pairs;
  pairs[0] = std::make_pair(act_on_right(0,b), act_on_right(1,b));

  bool found_good_pair = false;

  int list_size_max = 1000;
  
  if (abs(z) > 1.0/sqrt(2.0)) return false;
  
  //std::cout << "Epsilon: " << epsilon << "\n";

  for (int i=2; i<=n_depth; ++i) {
  
    //check if any pair is close enough
    for (int j=0; j<(int)pairs.size(); ++j) {
      if ( abs(pairs[j].first.center - pairs[j].second.center) < 10*epsilon ) {
        found_good_pair = true;
        cpx c0 = pairs[j].first.center;
        cpx c1 = pairs[j].second.center;
        //std::cout << "Found the good pair: " << pairs[j].first << " " << pairs[j].second << "\n";
        //get rid of *all* pairs that involve either disk
        int k = 0;
        while (k < (int)pairs.size()) {
          if ( abs(pairs[k].first.center - c0) < 0.000001 ||
               abs(pairs[k].second.center - c1) < 0.000001 ) {
            pairs.erase(pairs.begin() + k);
          } else {
            ++k;
          }
        }
        break;
      }
    }
  
    //make all possible next pairs
    next_pairs.resize(4*pairs.size());
    for (int j=0; j<(int)pairs.size(); ++j) {
      Ball& bz = pairs[j].first;
      Ball& bw = pairs[j].second;
      Ball bzs[2] = {act_on_right(0, bz), act_on_right(1, bz)};
      Ball bws[2] = {act_on_right(0, bw), act_on_right(1, bw)};
      for (int k=0; k<4; ++k) {
        next_pairs[4*j + k] = std::make_pair( bzs[k>>1], bws[k&1] );
      }
    }
    pairs.resize(0);
    for (int j=0; j<(int)next_pairs.size(); ++j) {
      double d = abs(next_pairs[j].first.center - next_pairs[j].second.center);
      double r = next_pairs[j].first.radius;
      //make sure we don't carry around an absurd number
      if ( (d < 12*epsilon + 2*r) && ((int)pairs.size() < list_size_max)) {
        pairs.push_back(next_pairs[j]);
      }
      if (found_good_pair && d < 10*epsilon) return false;
    }  
  }
  
  //std::cout << "pairs left: " << pairs.size();
  
  //now, if found_good_pair and pairs is *empty*, then we're good
  return found_good_pair && (pairs.size() == 0);
}



void ifs::compute_uv_graph(std::vector<Point3d<int> >& uv_graph, 
                           std::vector<Ball>& balls, 
                           int uv_depth, 
                           int verbose) {
  //make all the balls
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return;
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  compute_balls_right(balls, initial_ball, uv_depth);
  
  //find the closest uv words
  std::vector<std::pair<Bitword,Bitword> > words;
  find_closest_uv_words(words, uv_depth, 0.000000001);
  if (verbose > 0) {
    std::cout << "Found the closest bitwords (depth " << uv_depth << "):\n";
    for (int i=0; i<(int)words.size(); ++i) {
      std::cout << words[i].first << " " << words[i].second << "\n";
    }
  }
  std::vector<Point2d<int> > closest_balls(0);
  std::vector<std::string> closest_ball_strs(0); //records the words starting with 0
  std::stringstream T;
  for (int i=0; i<(int)words.size(); ++i) {
    int bi1 = (int)words[i].first.w.to_ulong();
    int bi2 = (int)words[i].second.w.to_ulong();
    closest_balls.push_back(Point2d<int>(bi1,bi2));
    closest_ball_strs.push_back( words[i].first.w.to_string().substr(64-uv_depth, uv_depth) );
  }
  
  if (verbose>0) {
    std::cout << "Closest ball list:\n";
    for (int i=0; i<(int)closest_balls.size(); ++i) {
      std::cout << closest_balls[i] << " " << closest_ball_strs[i] << "\n";
    }
  }
  
  //go through and find all the edges
  //only draw an edge if it goes between i->j with i<j (to avoid double edges)
  //note that i (the ball index) *is* the ball word
  uv_graph.resize(0);
  for (int i=0; i<(int)balls.size(); ++i) {
    std::string bs = balls[i].word.to_string().substr(64-uv_depth, uv_depth);
    //std::cout << "Finding edges from ball " << i << " with name " << bs << "\n";
    for (int j=0; j<uv_depth; ++j) {
      int suf_len = uv_depth-j;
      //get the suffix; if it doesn't begin with a zero, continue
      std::string suf = bs.substr(j, suf_len);
      int suf_mask = 0;
      for (int k=0; k<suf_len; ++k) suf_mask |= (1 << k);
      //std::cout << "Looking for suffix " << suf << " with mask " << suf_mask << "\n";
      if (suf[0] != '0') continue; 
      for (int k=0; k<(int)closest_ball_strs.size(); ++k) {
        if (closest_ball_strs[k].substr(0, suf_len) == suf) {
          //get the first suf_len digits from the *other* ball
          int matching_closest_ball = closest_balls[k].y;
          int shifted_digits = matching_closest_ball >> j;
          int other_ball = (i ^ (i & suf_mask)) | shifted_digits;
          uv_graph.push_back(Point3d<int>(i,other_ball,suf_len));
          //std::cout << "Added the edge " << uv_graph.back() << " suffix: " << suf << " matching ball: " << matching_closest_ball <<"\n";
        }
      } 
    }
  }
}
  
  





//find the center of mass of the balls (I guess just average for now
cpx ifs::center_of_mass(const std::vector<Ball>& balls) {
  cpx a = 0.0;
  for (int i=0; i<(int)balls.size(); ++i) {
    a += balls[i].center;
  }
  return a;
}


//compute the minimal radius of a ball around 1/2 which contains 
//its image under both actions
bool ifs::minimal_enclosing_radius(double& r) {
  if (1.0 - az < 1e-4 || 1.0 - aw < 1e-4) return false;
  //initialize the chunky radius to contain the whole set
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
  double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  r = (z_restriction > w_restriction ? z_restriction : w_restriction);
  //std::cout << "z: " << z << " az: " << az << " w: " << w << " aw: " << aw << "\n"; 
  //std::cout << "Z restriction: " << z_restriction << "\n";
  //std::cout << "W restriction: " << w_restriction << "\n";
  //std::cout << "Computed minimal radius as " << r << "\n";
  return true;
}

cpx ifs::iterate(int index, cpx u){
  // apply generator to u
  if(index==0){
          return(u*z);
  } else {
          return(((u-1.0)*w)+1.0);
  }
}


//convert a complex to the point in the drawing
//1/2 is the center, and the box must contain [-1,2]x[-1.5,1.5]
Point2d<int> ifs::cpx_to_point(cpx w) {
  Point2d<int> p;
  double x = (w.real() + 1.0)/3.0;
  double y = (w.imag() + 1.5)/3.0;
  p.x = int(drawing_width*x);
  p.y = int(drawing_width*y);
  return p;
}

int ifs::cpx_to_radius(cpx w) {
  int r;
  r= int(drawing_radius*abs(w));
  return r;
}

//this produces a point in [-1,1]x[-1,1] depending 
//on where in the window it is
cpx ifs::point_to_cpx(const Point2d<int>& p) {
  cpx w;
  w = cpx( double(p.x-drawing_radius)/drawing_radius,
	         double(p.y-drawing_radius)/drawing_radius );
  return w;
}

Point2d<int> ifs::cpx_to_point_mandlebrot(cpx w) {
  cpx one_scaled = ((w-center)/wind); //this is a number in [-1,1]x[-1,1]
  return Point2d<int>( int( one_scaled.real()*drawing_radius + drawing_radius ),
                       int( one_scaled.imag()*drawing_radius + drawing_radius ));
}







void ifs::set_params(cpx Z, cpx W) {
  z = Z;
  az = abs(z);
  w = W;
  aw = abs(w);
}

void ifs::draw_ifs_to_array(std::vector<std::vector<Point3d<unsigned char> > >& bmp, 
                            const cpx& region_ll, const cpx& region_ur, int depth) {
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  
  //clear the array
  for (int i=0; i<(int)bmp.size(); ++i) {
    for (int j=0; j<(int)bmp[i].size(); ++j) {
      bmp[i][j] = Point3d<unsigned char>(255,255,255);
    }
  }
  
  double drawing_width = region_ur.real() - region_ll.real();
  int num_pixels = bmp.size();
  double pixel_width = drawing_width/double(num_pixels);
  
  std::deque<Ball> stack(1);
  stack[0] = initial_ball;
  while (stack.size() > 0) {
    Ball b = stack.back();
    stack.pop_back();
    if ( (depth > 0 && b.word_len >= depth) || 
         (depth < 0 && b.radius < pixel_width/2.0) ) {
      int x = int(num_pixels*((b.center.real() - region_ll.real()) / drawing_width));
      int y = int(num_pixels*((b.center.imag() - region_ll.imag()) / drawing_width));
      if (0 <= x && x < num_pixels && 0 <= y && y < num_pixels) {
        bmp[x][y] = (b.last_gen_index() == 0 ? Point3d<unsigned char>(0xFF, 0xAA, 0x00) :
                                               Point3d<unsigned char>(0x00, 0xAA, 0xFF));
      }
    } else {
      Ball bz = act_on_right(0, b);
      Ball bw = act_on_right(1, b);
      if (!is_ball_disjoint(bz, region_ll, region_ur)) stack.push_front(bz);
      if (!is_ball_disjoint(bw, region_ll, region_ur)) stack.push_front(bw);
    }
  }
    
}
























