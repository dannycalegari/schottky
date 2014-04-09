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
  crad = 0.5;
  radius = 1.0;
  word = 0;
  word_len = 0;
}

Ball::Ball(cpx c, cpx cr, double r) {
  center = c;
  crad = cr;
  radius = r;
  word = 0;
  word_len = 0; //just a ball, no words, to start
}

Ball::Ball(cpx c, cpx cr, double r, int w, int wl) {
  center = c;
  crad = cr;
  radius = r;
  word = w;
  word_len = wl; //just a ball, no words, to start
}

int Ball::last_gen_index() const {
  return (word >> (word_len-1))&1;
}

std::ostream& operator<<(std::ostream& os, const Ball& b) {
  return os << "Ball(" << b.center << "," << b.crad << "," << b.radius << "," << b.word << "," << b.word_len << ")";
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
Ball ifs::act_on_left(int index, const Ball& b) {
  int word = b.word;
  int word_len = b.word_len;
  if (index == 0) {
    return Ball( z*b.center, z*b.crad, az*b.radius, word, word_len+1 );
  } else {
    return Ball( (w*(b.center - 1.0)) + 1.0, w*b.crad, aw*b.radius, word | (1 << word_len), word_len+1 );
  }
}


Ball ifs::act_on_right(int index, const Ball& b) {
  int word = b.word;
  int word_len = b.word_len;
  if (index == 0) {
    //the new center should be at the end of -crad (which points to 1)
    return Ball( b.center - b.crad, z*b.crad, az*b.radius, word<<1, word_len+1 );
  } else {
    return Ball( b.center + b.crad, w*b.crad, aw*b.radius, (word<<1)|1, word_len+1 );
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



//compute the minimal radius of a ball around 1/2 which contains 
//its image under both actions
double ifs::minimal_enclosing_radius() {
  //initialize the chunky radius to contain the whole set
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
  double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  return (z_restriction > w_restriction ? z_restriction : w_restriction);
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


cpx ifs::point_to_cpx(const Point2d<int>& p) {
  cpx w;
  w = cpx( double(p.x-drawing_radius)/drawing_radius,
	   double(p.y-drawing_radius)/drawing_radius );
  return w;
};



