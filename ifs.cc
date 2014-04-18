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


//first some ball functions
Ball::Ball() { 
  center = 0.5;
  to_z = to_w = 0;
  radius = 1.0;
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

std::ostream& operator<<(std::ostream& os, const Ball& b) {
  return os << "Ball(" << b.center << "," << b.to_z << "," << b.to_w << "," << b.radius << "," << b.word << "," << b.word_len << ")";
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
    balls[j] = act_on_right(0, parent_ball);
    balls[(j | (1<<current_depth))] = act_on_right(1, parent_ball);
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







//given a ball, check if it is disjoint from the square given
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
  if (fabs(az-1.0) < 1e-4 || fabs(aw-1.0) < 1e-4) return false;
  //initialize the chunky radius to contain the whole set
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
  double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  r = (z_restriction > w_restriction ? z_restriction : w_restriction);
  return true;
  //cout << "z: " << z << " az: " << az << " w: " << w << " aw: " << aw << "\n"; 
  //cout << "Z restriction: " << z_restriction << "\n";
  //cout << "W restriction: " << w_restriction << "\n";
  //cout << "Computed minimal radius as " << chunky_radius << "\n";
}

cpx ifs::iterate(int index, cpx u){
	// apply generator to u
	if(index==0){
		return(u*z);
	} else {
		return(((u-1.0)*w)+1.0);
	};
};


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







































