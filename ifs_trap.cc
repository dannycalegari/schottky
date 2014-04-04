//Functions to check for a trap


//this records a ball
struct Ball {
  Ball() { 
    center = 0.5;
    radius = 1.0;
  }
  Ball(cpx c, double r) {
    center = c;
    radius = r;
  }
  cpx center;
  double radius;
};

//compute the image ball
//(left action)
Ball ifs::act_on_left(int index, const Ball& b) {
  if (index == 0) {
    return Ball( z*b.center, az*b.radius );
  } else {
    return Ball( (w*(b.center - 1.0)) + 1.0, aw*b.radius );
  }
}


//Trap computation grid functions

struct GridPixel {
  cpx center; //the center of this pixel
  
  //the data on balls *containing* the pixel
  //bool contained_in_z;
  int z_ball_containing;
  //double z_containing_distance; 
  //bool contained_in_w;
  int w_ball_containing;
  //double w_containing_distance;
  
  //the data on balls *touching* the pixel
  //bool touching_z;
  int z_ball_touching;
  double z_touching_distance;
  //bool touching_w;
  int w_ball_touching;
  double w_touching_distance;
};
  
  
struct TrapGrid {
  int pixel_width; //the number of pixels on a side(it's a square)
  cpx lower_left;
  cpx upper_right;
  double pixel_radius; //half the height of the pixels
  
  //the pixel grid grid[i][j] is the ith pixel in the x direction (horizontal), and 
  //the jth in the y direction.  grid[0][0] is the lower left pixel
  std::vector<std::vector<GridPixel> > grid;
  
  TrapGrid() {
    width = 0;
    grid.resize(0);
  }
  
  TrapGrid(const std::vector<Ball>& balls, int w, double rad_mul);
  void reset_grid(cpx ll, cpx ur);
  void fill_pixels(const std::vector<Ball>& balls, double rad_mul);
  void fill_pixels(const Ball& b, double rad_mul); 
  void pixel_indices(int& i, int&j, const cpx& u); //find the pixel with the given point
  cpx pixel_center(int i, int j);  //return the point which is the center of the pixel
  cpx pixel_center(const cpx& u);  //return the point which is the center of pixel containing the given point
};


//initialize the grid by finding the min/max of x,y coordinates
//and then filling in the balls
TrapGrid::TrapGrid(const std::vector<Ball>& balls, 
                   int w,
                   double rad_mul) {
  pixel_width = w;
  lower_left = cpx(100,100);
  upper_right = cpx(-100,-100);
  int nb = (int)balls.size();
  //compute the max/min
  for (int i=0; i<nb; ++i) {
    cpx rad_vec = rad_mul*cpx(balls[i].radius, balls[i].radius);
    cpx temp_ll = balls[i].center - rad_vec;
    cpx temp_ur = balls[i].center + rad_vec;
    if (temp_ll.real() < lower_left.real()) {
      lower_left = cpx(temp_ll.real(), lower_left.imag());
    }
    if (temp_ll.imag() < lower_left.imag()) {
      lower_left = cpx(lower_left.real(), temp_ll.imag());
    }
    if (temp_ur.real() > upper_right.real()) {
      upper_right = cpx(temp_ur.real(), upper_right.imag());
    }
    if (temp_ur.imag() > upper_right.imag()) {
      upper_right - cpx(upper_right.imag(), temp_ur.imag());
    }
  }
  //make it a square so we don't go insane?
  double putative_width = upper_right.real() - lower_left.real();
  double putative_height = upper_right.imag() - lower_left.imag();
  if (putative_width > putative_height) { //the width is larger -- use it
    double height_adjustment = (putative_width-putative_height)/2.0;
    upper_right = cpx(upper_right.real(), upper_right.imag() + height_adjustment);
    lower_left = cpx(lower_left.real(), lower_left.imag() - height_adjustment);
  } else { //the height is larger -- use it
    double width_adjustment = (putative_height-putative_width)/2.0;
    upper_right = cpx(upper_right.real() + width_adjustment, upper_right.imag());
    lower_left = cpx(lower_left.real() - width_adjustment, lower_left.imag());
  }
  
  //build the blank grid
  
  
  //initialize the pixels
  reset_grid(lower_left, upper_right);
  
  //fill in the pixels
  fill_pixels(balls, rad_mul);
  
}

//set the lower left and upper right, clear all the pixels, and 
//set their centers
void TrapGrid::reset_grid(cpx ll, cpx ur) {
}

//fill all the pixels with the given balls
void TrapGrid::fill_pixels(const std::vector<Ball>& balls, 
                           double rad_mul) {
  int nb = (int)balls.size();
  for (int i=0; i<nb; ++i) {
    fill_pixels(G, balls[i], rad_mul);
  }
}

//fill the pixels for a given ball
//this just naively checks everything in a square
void TrapGrid::fills_pixels(const Ball& ball, double rad_mul) {
  double r = ball.radius * rad_mul;
  double c = ball.center;
  
}

//find the pixel containing the given point
void TrapGrid::pixel_indices(int& i, int&j, const cpx& u) {
}

//return the point which is the center of the pixel
cpx TrapGrid::pixel_center(int i, int j) {
  return grid[i][j].center;
}

//return the point which is the center of pixel containing the given point
cpx TrapGrid::pixel_center(const cpx& u) {
  
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


bool ifs::find_trap() {
  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
  double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  double min_initial_radius = (z_restriction > w_restriction 
                                             ? z_restriction 
                                             : w_restriction);
  
  
  //starting depth will always be 8?
  int current_depth = 8;
  std::vector<Ball> balls;
  compute_balls(balls, Ball(0.5, min_initial_radius), current_depth);
  
  //now we need to effectively increase the radius of the 
  //*previous* step until the balls of the *current* step 
  //form a connected set for both f and g
  //we'll start by multiplying the radius by 1.1 (it's probably pretty small)
  double more_dramatic_action = (az < aw ? az : aw);
  double min_prev_rad = pow( more_dramatic_action, current_depth-1 ) * min_initial_radius;
  double prev_rad_mul = 1.1;
  
  //epsilon is the amount we *know* the disks contain an epsilon nbhd of L
  //the smallest absolute increase will be on the smallest ball
  double epsilon = min_prev_rad*(prev_rad_mul-1.0);
  
  //now fill out the grid using this prev_rad_mul
  //grid size is always 512?
  TrapGrid TG(*this, balls, 512, prev_rad_mul);

  
  
  return true;
}

void ifs::draw_trap() {
}















