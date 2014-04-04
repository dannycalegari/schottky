//Functions to check for a trap

#include "XGraphics.h"

//this records a ball
struct Ball {
  cpx center;
  double radius;
  int word;     //the f and g word we applied to get here (highest bit is most recent application)
  int word_len; //the number of significant bits
  Ball() { 
    center = 0.5;
    radius = 1.0;
    word = 0;
    word_len = 0;
  }
  Ball(cpx c, double r) {
    center = c;
    radius = r;
    word = 0;
    word_len = 0; //just a ball, no words, to start
  }
  Ball(cpx c, double r, int w, int wl) {
    center = c;
    radius = r;
    word = w;
    word_len = wl; //just a ball, no words, to start
  }
  int last_gen_index() const {
    return (word >> (word_len-1))&1;
  }
};

//compute the image ball
//(left action)
Ball ifs::act_on_left(int index, const Ball& b) {
  int word = b.word;
  int word_len = b.word_len;
  if (index == 0) {
    return Ball( z*b.center, az*b.radius, word, word_len+1 );
  } else {
    return Ball( (w*(b.center - 1.0)) + 1.0, aw*b.radius, word | (1 << word_len), word_len+1 );
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
  int num_pixels; //the number of pixels on a side(it's a square)
  cpx lower_left;  //the lower left corner of the lower left pixel
  cpx upper_right; //the upper right corner of the upper right pixel
  double pixel_diameter; //the height/width of a pixel
  double box_width; //the real width/height of the grid (upper_right.real() - lower_left.real())
  
  //the pixel grid grid[i][j] is the ith pixel in the x direction (horizontal), and 
  //the jth in the y direction.  grid[0][0] is the lower left pixel
  std::vector<std::vector<GridPixel> > grid;
  
  TrapGrid() {
    num_pixels = 0;
    grid.resize(0);
  }
  
  TrapGrid(const std::vector<Ball>& balls, int w, double rad_mul);
  void reset_grid(cpx ll, cpx ur);
  void fill_pixels(const std::vector<Ball>& balls, double rad_mul); 
  void pixel_indices(int& i, int&j, const cpx& u); //find the pixel with the given point
  cpx pixel_center(int i, int j);  //return the point which is the center of the pixel
  cpx pixel_center(const cpx& u);  //return the point which is the center of pixel containing the given point
};


//initialize the grid by finding the min/max of x,y coordinates
//and then filling in the balls
TrapGrid::TrapGrid(const std::vector<Ball>& balls, 
                   int w,
                   double rad_mul) {
  num_pixels = w;
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
  grid = std::vector<std::vector<GridPixel> >(num_pixels, std::vector<GridPixel>(num_pixels));
  
  //initialize the pixels
  reset_grid(lower_left, upper_right);
  
  //fill in the pixels
  fill_pixels(balls, rad_mul);
  
}

//set the lower left and upper right, clear all the pixels, and 
//set their centers
void TrapGrid::reset_grid(cpx ll, cpx ur) {
  lower_left = ll;
  upper_right = ur;
  box_width = ur.real() - ll.real();
  if (box_width < 1e-10) std::cout << "Grid precision is too low!\n";
  pixel_diameter = box_width/double(num_pixels);
  for (int i=0; i<num_pixels; ++i) {
    double rp = lower_left.real() + (double(i)+0.5)*pixel_diameter;
    for (int j=0; j<num_pixels; ++j) {
      grid[i][j].center = cpx(rp, lower_left.imag() + (double(j)+0.5)*pixel_diameter);
      grid[i][j].z_ball_containing = -1;
      grid[i][j].w_ball_containing = -1;
      grid[i][j].z_ball_touching = -1;
      grid[i][j].z_touching_distance = -1;
      grid[i][j].w_ball_touching = -1;
      grid[i][j].w_touching_distance = -1;
    }
  }
}

//fill the pixels for a given ball
//this just naively checks everything in a square
void TrapGrid::fill_pixels(const std::vector<Ball>& balls, double rad_mul) {
  int nb = (int)balls.size();
  for (int bi=0; bi<nb; ++bi) {
    double r = balls[bi].radius * rad_mul;
    const cpx& c = balls[bi].center;
    int ll_touching_i, ll_touching_j;
    int ur_touching_i, ur_touching_j;
    double pixel_radius = 0.5*pixel_diameter;
    pixel_indices(ll_touching_i, ll_touching_j, cpx(c.real()-r, c.imag()-r));
    pixel_indices(ur_touching_i, ur_touching_j, cpx(c.real()+r, c.imag()+r));
    for (int i=ll_touching_i; i<=ur_touching_i; ++i) {
      for (int j=ll_touching_j; j<=ur_touching_j; ++j) {
        //determine the quadrant that grid[i][j] is in relative to the center
        //first, the four most likely things:
        cpx p_c = grid[i][j].center;
        cpx p_ur(p_c.real() + pixel_radius, p_c.imag() + pixel_radius);
        cpx p_lr(p_c.real() + pixel_radius, p_c.imag() - pixel_radius);
        cpx p_ul(p_c.real() - pixel_radius, p_c.imag() + pixel_radius);
        cpx p_ll(p_c.real() - pixel_radius, p_c.imag() - pixel_radius);
        bool touching = false;
        bool containing = false;
        if (p_ur.real() < c.real() && p_ur.imag() < c.imag()) { 
          //pixel is to the lower left of the center
          //it suffices to check the upper right for touching and lower left for containing
          if (abs(c-p_ll) < r) containing = true;
          else if (abs(c-p_ur) < r) touching = true;
                
        } else if (p_lr.real() < c.real() && p_lr.imag() > c.imag()) {
          //pixel is to the upper left
          //it suffices to check lr for touching and upper left for containing
          if (abs(c-p_ul) < r) containing = true;
          else if (abs(c-p_lr) < r) touching = true;
          
        } else if (p_ll.real() > c.real() && p_ll.imag() > c.imag()) {
          //pixel is to the upper right
          //suffices to check ll for touching and ur for containing
          if (abs(c-p_ur) < r) containing = true;
          else if (abs(c-p_ll) < r) touching = true;
          
        } else if (p_ul.real() > c.real() && p_ul.imag() < c.imag()) {
          //pixel is to the lower right
          //it suffices to check ul for touching and lr for containing
          if (abs(c-p_lr) < r) containing = true;
          else if (abs(c-p_ul) < r) touching = true;
        
        } else {
          //those were the four easier cases
          if (p_ur.real() < c.real()) { 
            //the pixel is to the left
            //it suffices to check that p_ll.real() > c.real() - r for touching
            //and that it contains p_ll and p_ul
            if (abs(c-p_ll) < r && abs(c-p_ul)) containing = true;
            else if (p_ll.real() > c.real() - r) touching = true;
          
          } else if (p_ll.imag() > c.imag()) {
            //the pixel is above
            //it suffices to check that p_ll.imag() < c.imag() + r for touching
            //and that it contains the top ones for containing
            if (abs(c-p_ur)<r && abs(c-p_ul)<r) containing = true;
            else if (p_ll.imag() > c.imag() + r) touching = true;
            
          } else if (p_ll.real() > c.real()) {
            //the pixel is to the right
            //it suffices to check that p_ll.real() < c.real() + r for touching 
            //and that it contains the right ones for containing
            if (abs(c-p_ur)<r && abs(c-p_lr)<r) containing = true;
            else if (p_ll.real() > c.real() + r) touching = true;
            
          } else if (p_ur.imag() < c.imag()) {
            //the pixel is below
            //it suffices to check that p_ur.imag() > c.imag() - r for touching and 
            //that it contains the bottom ones for containing
            if (abs(c-p_ll)<r && abs(c-p_lr)<r) containing = true;
            else if (p_ur.imag() > c.imag() -r) touching = true;
          
          } else { 
            //the pixel is directly over the center of the circle
            //it *is* touching, and its containing if all the corners are contained 
            if (abs(c-p_ll)<r && abs(c-p_lr)<r && abs(c-p_ur)<r && abs(c-p_ul)<r) containing = true;
            else touching = true;
          }
        }
        if (containing) touching = true;
        if (balls[bi].last_gen_index() == 0) {                //it's in the z set
          if (containing) grid[i][j].z_ball_containing = bi;  //might as well overwrite it
          if (touching) {
            double d = abs(p_c-c);
            if (!grid[i][j].z_ball_touching || d < grid[i][j].z_touching_distance) {
              grid[i][j].z_ball_touching = bi;
              grid[i][j].z_touching_distance = d;
            }
          }
        } else {                                              //it's in the w set
          if (containing) grid[i][j].w_ball_containing = bi;  //might as well overwrite it
          if (touching) {
            double d = abs(p_c-c);
            if (!grid[i][j].w_ball_touching || d < grid[i][j].w_touching_distance) {
              grid[i][j].w_ball_touching = bi;
              grid[i][j].w_touching_distance = d;
            }
          }
        }
      }
    }//<-- end of loop over pixels
    
  }//<-- end of loop over balls
  
}


//find the pixel containing the given point
void TrapGrid::pixel_indices(int& i, int&j, const cpx& u) {
  double h_dist = u.real() - lower_left.real();
  double v_dist = u.imag() - lower_left.imag();
  if (h_dist < 0 || v_dist < 0 
                 || upper_right.real() - u.real() < 0 
                 || upper_right.imag() - u.imag() < 0) { //point is outside the grid
    i=-1;
    j=-1;
  }
  i = int(h_dist/pixel_diameter);
  j = int(v_dist/pixel_diameter);
}

//return the point which is the center of the pixel
cpx TrapGrid::pixel_center(int i, int j) {
  return grid[i][j].center;
}

//return the point which is the center of pixel containing the given point
cpx TrapGrid::pixel_center(const cpx& u) {
  int i,j;
  pixel_indices(i,j,u);
  return grid[i][j].center;
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
  
  
  //starting depth will always be 10?
  int current_depth = 10;
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
  TrapGrid TG(balls, 512, prev_rad_mul);
  
  TG.show();

  
  
  return true;
}

void ifs::draw_trap() {
}















