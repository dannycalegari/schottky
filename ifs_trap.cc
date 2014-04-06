//Functions to check for a trap


//Trap computation grid functions

struct GridPixel {
  cpx center; //the center of this pixel
  
  //the z and w ball status
  int z_ball_status; //0=not touching, 1=touching, 2=containing
  int closest_z_ball;
  double z_ball_distance;
  
  int w_ball_status;
  int closest_w_ball;
  double w_ball_distance;
  
  //which components this pixel is in
  int z_comp, z_cut_by_w_comp;
  int w_comp, w_cut_by_z_comp;
  int intersection_comp;
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
  
  //Connected component lists:
  //these are lists of pixels in the connected components of the pixels 
  //which are contained in z balls
  std::vector<std::vector<Point2d<int> > > z_components;
  //these are lists of pixels in the connected components of the pixels
  //which are contained in the w balls 
  std::vector<std::vector<Point2d<int> > > w_components;
  //these are the lists of pixels in the connected components of the pixels
  //contained in z balls and also not contained in w balls
  std::vector<std::vector<Point2d<int> > > z_cut_by_w_components;
  //these are the lists of pixels in the connected components of the pixels
  //contained in the w balls and also not contained in the z balls
  std::vector<std::vector<Point2d<int> > > w_cut_by_z_components;
  //these are connected components of intersection pixels 
  //(at least touched by z and w)
  std::vector<std::vector<Point2d<int> > > intersection_components;
  
  TrapGrid() {
    num_pixels = 0;
    grid.resize(0);
  }
  
  TrapGrid(const std::vector<Ball>& balls, int max_num_pixels, double rad_mul);
  void reset_grid(cpx ll, cpx ur);
  void fill_pixels(const std::vector<Ball>& balls, double rad_mul); 
  void pixel_indices(int& i, int&j, const cpx& u); //find the pixel with the given point
  cpx pixel_center(int i, int j);  //return the point which is the center of the pixel
  cpx pixel_center(const cpx& u);  //return the point which is the center of pixel containing the given point
  void compute_connected_components(); //compute the connected components
  void pursue_z_comp(int i, int j, int ind);
  void pursue_w_comp(int i, int j, int ind);
  void pursue_z_cut_by_w_comp(int i, int j, int ind);
  void pursue_w_cut_by_z_comp(int i, int j, int ind);
  void pursue_intersection_comp(int i, int j, int ind);
  void show();
};


//initialize the grid by finding the min/max of x,y coordinates
//and then filling in the balls
TrapGrid::TrapGrid(const std::vector<Ball>& balls, 
                   int max_num_pixels,
                   double rad_mul) {
  lower_left = cpx(100,100);
  upper_right = cpx(-100,-100);
  int nb = (int)balls.size();
  double min_radius = 1000;
  double max_radius = -1;
  //compute the max/min
  for (int i=0; i<nb; ++i) {
    double r = rad_mul*balls[i].radius;
    if (r < min_radius) min_radius = r;
    else if (r > max_radius) max_radius = r;
    cpx rad_vec = cpx(r, r);
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
      upper_right = cpx(upper_right.imag(), temp_ur.imag());
    }
  }
  //make it a square so we don't go insane?
  double putative_width = upper_right.real() - lower_left.real();
  double putative_height = upper_right.imag() - lower_left.imag();
  std::cout << "Putative ll, ur: " << lower_left << ", " << upper_right << "\n";
  std::cout << "Putative width and height: " << putative_width << ", " << putative_height << "\n";
  if (putative_width > putative_height) { //the width is larger -- use it
    double height_adjustment = (putative_width-putative_height)/2.0;
    std::cout << "Height adjustment: " << height_adjustment << "\n";
    upper_right = cpx(upper_right.real(), upper_right.imag() + height_adjustment);
    lower_left = cpx(lower_left.real(), lower_left.imag() - height_adjustment);
  } else { //the height is larger -- use it
    double width_adjustment = (putative_height-putative_width)/2.0;
    std::cout << "Width adjustment: " << width_adjustment << "\n";
    upper_right = cpx(upper_right.real() + width_adjustment, upper_right.imag());
    lower_left = cpx(lower_left.real() - width_adjustment, lower_left.imag());
  }
  std::cout << "Grid will run " << lower_left << " -- " << upper_right << "\n";
  
  //decide how many pixels to use
  //we want the smallest disk to contain 4ish pixels, so 
  //pixel_diameter should be at most like radius/1.5 ish
  if (max_radius / min_radius > 10) std::cout << "Badly balanced radii\n";
  num_pixels = int( (upper_right.real() - lower_left.real())/(min_radius/1.5) );
  
  std::cout << "num_pixels: " << num_pixels << "\n";
  
  //initialize the pixels
  reset_grid(lower_left, upper_right);
  
  //fill in the pixels
  fill_pixels(balls, rad_mul);
  
}

//set the lower left and upper right, decide how many pixels to have,
//clear all the pixels, set their centers, clear all the connected components
void TrapGrid::reset_grid(cpx ll, cpx ur) {
  lower_left = ll;
  upper_right = ur;
  box_width = ur.real() - ll.real();

  //build the blank grid
  grid = std::vector<std::vector<GridPixel> >(num_pixels, std::vector<GridPixel>(num_pixels));

  //reset the connected components
  z_components.resize(0); z_cut_by_w_components.resize(0);
  w_components.resize(0); w_cut_by_z_components.resize(0);
  intersection_components.resize(0);

  pixel_diameter = box_width/double(num_pixels);
  if (pixel_diameter < 1e-10) std::cout << "Grid precision [" << ll << "," << ur << "] is too low!\n";
  std::cout << "Box width: " << box_width << " and pixel diameter: " << pixel_diameter << "\n";
  for (int i=0; i<num_pixels; ++i) {
    double rp = lower_left.real() + (double(i)+0.5)*pixel_diameter;
    for (int j=0; j<num_pixels; ++j) {
      grid[i][j].center = cpx(rp, lower_left.imag() + (double(j)+0.5)*pixel_diameter);
      grid[i][j].z_ball_status = 0;
      grid[i][j].closest_z_ball = -1;
      grid[i][j].z_ball_distance = -1;
      grid[i][j].w_ball_status = 0;
      grid[i][j].closest_w_ball = -1;
      grid[i][j].w_ball_distance = -1;
    }
  }
}

//fill the pixels for a given ball
//this just naively checks everything in a square
void TrapGrid::fill_pixels(const std::vector<Ball>& balls, double rad_mul) {
  int nb = (int)balls.size();
  for (int bi=0; bi<nb; ++bi) {
    double r = balls[bi].radius * rad_mul;
    const cpx c = balls[bi].center;
    //std::cout << "Filling pixels for ball at " << c << ", " << r << "\n";
    int ll_touching_i, ll_touching_j;
    int ur_touching_i, ur_touching_j;
    double pixel_radius = 0.5*pixel_diameter;
    pixel_indices(ll_touching_i, ll_touching_j, cpx(c.real()-r, c.imag()-r));
    pixel_indices(ur_touching_i, ur_touching_j, cpx(c.real()+r, c.imag()+r));
    for (int i=ll_touching_i; i<=ur_touching_i; ++i) {
      if (i >= num_pixels) break;
      else if (i < 0) continue; 
      for (int j=ll_touching_j; j<=ur_touching_j; ++j) {
        if (j >= num_pixels) break;
        else if (j < 0) continue;
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
            //it suffices to check that p_lr.real() > c.real() - r for touching
            //and that it contains p_ll and p_ul for containing
            if (abs(c-p_ll) < r && abs(c-p_ul)<r) containing = true;
            else if (p_lr.real() > c.real() - r) touching = true;
          
          } else if (p_ll.imag() > c.imag()) {
            //the pixel is above
            //it suffices to check that p_ll.imag() < c.imag() + r for touching
            //and that it contains the top ones for containing
            if (abs(c-p_ur)<r && abs(c-p_ul)<r) containing = true;
            else if (p_ll.imag() < c.imag() + r) touching = true;
            
          } else if (p_ll.real() > c.real()) {
            //the pixel is to the right
            //it suffices to check that p_ll.real() < c.real() + r for touching 
            //and that it contains the right ones for containing
            if (abs(c-p_ur)<r && abs(c-p_lr)<r) containing = true;
            else if (p_ll.real() < c.real() + r) touching = true;
            
          } else if (p_ur.imag() < c.imag()) {
            //the pixel is below
            //it suffices to check that p_ur.imag() > c.imag() - r for touching and 
            //that it contains the bottom ones for containing
            if (abs(c-p_ll)<r && abs(c-p_lr)<r) containing = true;
            else if (p_ur.imag() > c.imag() - r) touching = true;
          
          } else { 
            //the pixel is directly over the center of the circle
            //it *is* touching, and its containing if all the corners are contained 
            if (abs(c-p_ll)<r && abs(c-p_lr)<r && abs(c-p_ur)<r && abs(c-p_ul)<r) containing = true;
            else touching = true;
          }
        }
        //std::cout << "Found pixel " << i << "," << j << " with center " << p_c << " is: ";
        //if (containing) std::cout << " containing\n";
        //else if (touching) std::cout << "touching\n";
        //else std::cout << "neither\n";
        if (containing) touching = true;
        if (balls[bi].last_gen_index() == 0) {                //it's in the z set
          if (containing) grid[i][j].z_ball_status = 2;  //might as well overwrite it
          if (touching) {
            if (grid[i][j].z_ball_status == 0) grid[i][j].z_ball_status = 1;
            double d = abs(p_c-c);
            if (grid[i][j].closest_z_ball<0 || d < grid[i][j].z_ball_distance) {
              grid[i][j].closest_z_ball = bi;
              grid[i][j].z_ball_distance = d;
            }
          }
        } else {                                              //it's in the w set
          if (containing) grid[i][j].w_ball_status = 2;  //might as well overwrite it
          if (touching) {
            if (grid[i][j].w_ball_status == 0) grid[i][j].w_ball_status = 1;
            double d = abs(p_c-c);
            if (grid[i][j].closest_w_ball<0 || d < grid[i][j].w_ball_distance) {
              grid[i][j].closest_w_ball = bi;
              grid[i][j].w_ball_distance = d;
            }
          }
        }
      }
    }//<-- end of loop over pixels
    
  }//<-- end of loop over balls
  
}


//fill out the connected components
//during this computation, the grid[i][j].*comp is -1 if nothing 
//has touched it, -2 if it's been stacked, and >=0 if it records 
//which component
void TrapGrid::compute_connected_components() { 
  //reset the connected components
  z_components.resize(0); z_cut_by_w_components.resize(0);
  w_components.resize(0); w_cut_by_z_components.resize(0);
  intersection_components.resize(0);
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      grid[i][j].z_comp = grid[i][j].w_comp = 
      grid[i][j].z_cut_by_w_comp = grid[i][j].w_cut_by_z_comp = 
      grid[i][j].intersection_comp = -1;
    }
  }
  //go through and pursue all the components we encounter
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      //does it form a new z component?
      if (grid[i][j].z_ball_status==2 && grid[i][j].z_comp == -1) {
        z_components.push_back(std::vector<Point2d<int> >(0));
        int ind = z_components.size()-1;
        pursue_z_comp(i,j, ind);
        if (grid[i][j].w_ball_status==0 && grid[i][j].z_cut_by_w_comp == -1) {
          z_cut_by_w_components.push_back(std::vector<Point2d<int> >(0));
          int ind = z_cut_by_w_components.size()-1;
          pursue_z_cut_by_w_comp(i,j,ind);
        }
      }
      if (grid[i][j].w_ball_status==2 && grid[i][j].w_comp == -1) {
        w_components.push_back(std::vector<Point2d<int> >(0));
        int ind = w_components.size()-1;
        pursue_w_comp(i,j, ind);
        if (grid[i][j].z_ball_status==0 && grid[i][j].w_cut_by_z_comp == -1) {
          w_cut_by_z_components.push_back(std::vector<Point2d<int> >(0));
          int ind = w_cut_by_z_components.size()-1;
          pursue_w_cut_by_z_comp(i,j,ind);
        }
      }
      if (grid[i][j].z_ball_status > 0 && grid[i][j].w_ball_status > 0) {
        intersection_components.push_back(std::vector<Point2d<int> >(0));
        int ind = intersection_components.size()-1;
        pursue_intersection_comp(i,j,ind);
      }
    }
  }
}

void TrapGrid::pursue_z_comp(int i, int j, int ind) {
}
void TrapGrid::pursue_w_comp(int i, int j, int ind) {
}
void TrapGrid::pursue_z_cut_by_w_comp(int i, int j, int ind) {
}
void TrapGrid::pursue_w_cut_by_z_comp(int i, int j, int ind) {
}
void TrapGrid::pursue_intersection_comp(int i, int j, int ind) {
}

//find the pixel containing the given point
//there is *no error checking!*
void TrapGrid::pixel_indices(int& i, int&j, const cpx& u) {
  double h_dist = u.real() - lower_left.real();
  double v_dist = u.imag() - lower_left.imag();
  //if (h_dist < 0 || v_dist < 0 
  //               || upper_right.real() - u.real() < 0 
  //               || upper_right.imag() - u.imag() < 0) { //point is outside the grid
  //  i=-1;
  //  j=-1;
  //}
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





bool ifs::find_trap() {
  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double min_initial_radius = minimal_enclosing_radius();
  
  
  //starting depth will always be 10?
  int current_depth = depth; //10;
  std::vector<Ball> balls;
  compute_balls(balls, Ball(0.5, min_initial_radius), current_depth);
  
  //std::cout << "Computed " << balls.size() << " balls\n";
  //for (int i=0; i<(int)balls.size(); ++i) {
  //  std::cout << i << ": " << balls[i] << "\n";
  //}

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
  
  //find the connected components
  TG.compute_connected_components();
  
  TG.show();

  
  
  return true;
}


void TrapGrid::show() {
  int pixel_group_width = -1;
  int num_drawing_pixels = -1;
  if (num_pixels < 512) {
    pixel_group_width = 512 / num_pixels;
    num_drawing_pixels = pixel_group_width*num_pixels;
  } else {
    pixel_group_width = 1;
    num_drawing_pixels = num_pixels;
  }
  XGraphics X2(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
  Point2d<int> p;
  int zt = X2.get_rgb_color(0,0,0.5);
  int zc = X2.get_rgb_color(0,0,1);
  int wt = X2.get_rgb_color(0,0.5,0);
  int wc = X2.get_rgb_color(0,1,0);
  int ztwt = X2.get_rgb_color(0.5,0,0);
  int ztwc = X2.get_rgb_color(0.5,0,0);
  int zcwt = X2.get_rgb_color(0.5,0,0);
  int zcwc = X2.get_rgb_color(1,0,0);
  int whi = X2.get_rgb_color(1,1,1);
  for (int i=0; i<num_pixels; ++i) {
    p.x = pixel_group_width*i;
    for (int j=0; j<num_pixels; ++j) {
      p.y = pixel_group_width*j;
      int col = 0;
      if (grid[i][j].z_ball_status==2) {
        if (grid[i][j].w_ball_status==2) {
          col = zcwc;
        } else if (grid[i][j].w_ball_status==1) {
          col = zcwt;
        } else {
          col = zc;
        }
      } else if (grid[i][j].z_ball_status == 1) {
        if (grid[i][j].w_ball_status==2) {
          col = ztwc;
        } else if (grid[i][j].w_ball_status==1) {
          col = ztwt;
        } else {
          col = zt;
        }
      } else {
        if (grid[i][j].w_ball_status==2) {
          col = wc;
        } else if (grid[i][j].w_ball_status==1) {
          col = wt;
        } else {
          col = whi;
        }
      }
      if (num_pixels<512) {
        X2.draw_box(p, pixel_group_width, col);
      } else {
        X2.draw_point(p, col);
      }
    }
  }
  (void)X2.wait_for_key();
}




void ifs::draw_trap() {

}




