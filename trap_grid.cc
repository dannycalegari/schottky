#include <deque>

#include "graphics.h"

#include "trap_grid.h"



/***************************************************************************
 * TrapGrid functions
 ****************************************************************************/

//blank TG
TrapGrid::TrapGrid() {
}


//set the lower left and upper right, decide how many pixels to have,
//clear all the pixels, set their centers, clear all the connected components
void TrapGrid::reset_grid(cpx ll, cpx ur, int np) {
  lower_left = ll;
  upper_right = ur;
  box_width = ur.real() - ll.real();
  num_pixels = np;

  //build the blank grid
  grid = std::vector<std::vector<GridPixel> >(num_pixels, std::vector<GridPixel>(num_pixels));

  //reset the connected components
  z_components.resize(0); z_cut_by_w_components.resize(0);
  w_components.resize(0); w_cut_by_z_components.resize(0);
  intersection_components.resize(0);

  pixel_diameter = box_width/double(num_pixels);
  if (pixel_diameter < 1e-10) std::cout << "Grid precision [" << ll << "," << ur << "] is too low!\n";
  //std::cout << "Box width: " << box_width << " and pixel diameter: " << pixel_diameter << "\n";
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
void TrapGrid::fill_pixels(const std::vector<Ball>& balls) {
  int nb = (int)balls.size();
  for (int bi=0; bi<nb; ++bi) {
    double r = balls[bi].radius;
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



//decide if the ball b is disjoint from z or w pixels
//unfortunately, this seems like it needs to be just as complicated 
//as the function above
bool TrapGrid::disjoint_from_z_or_w(const Ball& b, int z_or_w) {
  double r = b.radius;
  const cpx c = b.center;
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
      cpx p_c = grid[i][j].center;
      cpx p_ur(p_c.real() + pixel_radius, p_c.imag() + pixel_radius);
      cpx p_lr(p_c.real() + pixel_radius, p_c.imag() - pixel_radius);
      cpx p_ul(p_c.real() - pixel_radius, p_c.imag() + pixel_radius);
      cpx p_ll(p_c.real() - pixel_radius, p_c.imag() - pixel_radius);
      if (p_ur.real() < c.real() && p_ur.imag() < c.imag()) { 
        if (abs(c-p_ur) < r) {
          if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
            return false;
          }
        }
              
      } else if (p_lr.real() < c.real() && p_lr.imag() > c.imag()) {
        //pixel is to the upper left
        //it suffices to check lr for touching and upper left for containing
        if (abs(c-p_lr) < r) {
          if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
            return false;
          }
        }
        
      } else if (p_ll.real() > c.real() && p_ll.imag() > c.imag()) {
        //pixel is to the upper right
        //suffices to check ll for touching and ur for containing
        if (abs(c-p_ll) < r) {
          if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
            return false;
          }
        }
        
      } else if (p_ul.real() > c.real() && p_ul.imag() < c.imag()) {
        //pixel is to the lower right
        //it suffices to check ul for touching and lr for containing
        if (abs(c-p_ul) < r) {
          if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
            return false;
          }
        }
      
      } else {
        //those were the four easier cases
        if (p_ur.real() < c.real()) { 
          //the pixel is to the left
          //it suffices to check that p_lr.real() > c.real() - r for touching
          //and that it contains p_ll and p_ul for containing
          if (p_lr.real() > c.real() - r) {
            if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
              return false;
            }
          }
        
        } else if (p_ll.imag() > c.imag()) {
          //the pixel is above
          //it suffices to check that p_ll.imag() < c.imag() + r for touching
          //and that it contains the top ones for containing
          if (p_ll.imag() < c.imag() + r) {
            if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
                (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
              return false;
            }
          }
          
        } else if (p_ll.real() > c.real()) {
          //the pixel is to the right
          //it suffices to check that p_ll.real() < c.real() + r for touching 
          //and that it contains the right ones for containing
          if (p_ll.real() < c.real() + r) {
            if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
              return false;
            }
          }
          
        } else if (p_ur.imag() < c.imag()) {
          //the pixel is below
          //it suffices to check that p_ur.imag() > c.imag() - r for touching and 
          //that it contains the bottom ones for containing
          if (p_ur.imag() > c.imag() - r) {
            if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
                (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
              return false;
            }
          }
        
        } else { 
          //the pixel is directly over the center of the circle
          //it *is* touching, and its containing if all the corners are contained 
          if ((z_or_w == 0 && grid[i][j].z_ball_status > 0) ||
              (z_or_w == 1 && grid[i][j].w_ball_status > 0)) {
            return false;
          }
        }
      }
    }
  }//<-- end of loop over pixels
  return true;
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
      grid[i][j].i_comp = -1;
    }
  }
  //go through and pursue all the components we encounter
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      //does it form a new z component?
      if (grid[i][j].z_ball_status==2 && grid[i][j].z_comp == -1) {
        z_components.push_back(std::vector<Point2d<int> >(0));
        pursue_z_comp(i,j, z_components.size()-1);
      }
      if (grid[i][j].z_ball_status==2 && 
          grid[i][j].z_cut_by_w_comp == -1 &&
          grid[i][j].w_ball_status==0) {
        z_cut_by_w_components.push_back(std::vector<Point2d<int> >(0));
        pursue_z_cut_by_w_comp(i,j,z_cut_by_w_components.size()-1);
      }
      //how about a w component?
      if (grid[i][j].w_ball_status==2 && grid[i][j].w_comp == -1) {
        w_components.push_back(std::vector<Point2d<int> >(0));
        pursue_w_comp(i,j, w_components.size()-1);
      }
      if (grid[i][j].w_ball_status==2 && 
          grid[i][j].w_cut_by_z_comp == -1 &&
          grid[i][j].z_ball_status==0) {
        w_cut_by_z_components.push_back(std::vector<Point2d<int> >(0));
        pursue_w_cut_by_z_comp(i,j,w_cut_by_z_components.size()-1);
      }
      //or intersection component
      if (grid[i][j].z_ball_status >0 && 
          grid[i][j].w_ball_status >0 &&
          grid[i][j].i_comp == -1) {
        intersection_components.push_back(std::vector<Point2d<int> >(0));
        pursue_intersection_comp(i,j,intersection_components.size()-1);
      }
    }
  }
}

//These are all the component-finding functions.
//There are slight differences, so there are multiple functions
//I guess these will use breadth-first search?
void TrapGrid::pursue_z_comp(int i, int j, int ind) {
  std::deque<Point2d<int> > stack(1, Point2d<int>(i,j));
  grid[i][j].z_comp = -2;
  while (stack.size() > 0) {
    Point2d<int> p = stack.back();
    stack.pop_back();
    //mark p as done
    grid[p.x][p.y].z_comp = ind;
    z_components[ind].push_back(p);
    //go around to the others  
    for (int ii=-1; ii<=1; ++ii) {
      int iii = ii + p.x;
      if (iii < 0) continue;
      else if (iii >= num_pixels) break;
      for (int jj=-1; jj<=1; ++jj) {
        if (jj == 0 && ii == 0) continue;
        int jjj = jj + p.y;
        if (jjj < 0) continue;
        else if (jjj >= num_pixels) break;
        if (grid[iii][jjj].z_comp == -1 && grid[iii][jjj].z_ball_status == 2) {
          grid[iii][jjj].z_comp = -2;
          stack.push_front(Point2d<int>(iii,jjj));
        }
      }
    }
  }
}

void TrapGrid::pursue_w_comp(int i, int j, int ind) {
  std::deque<Point2d<int> > stack(1, Point2d<int>(i,j));
  grid[i][j].w_comp = -2;
  while (stack.size() > 0) {
    Point2d<int> p = stack.back();
    stack.pop_back();
    //mark p as done
    grid[p.x][p.y].w_comp = ind;
    w_components[ind].push_back(p);
    //go around to the others  
    for (int ii=-1; ii<=1; ++ii) {
      int iii = ii + p.x;
      if (iii < 0) continue;
      else if (iii >= num_pixels) break;
      for (int jj=-1; jj<=1; ++jj) {
        if (jj == 0 && ii == 0) continue;
        int jjj = jj + p.y;
        if (jjj < 0) continue;
        else if (jjj >= num_pixels) break;
        if (grid[iii][jjj].w_comp == -1 && grid[iii][jjj].w_ball_status == 2) {
          grid[iii][jjj].w_comp = -2;
          stack.push_front(Point2d<int>(iii,jjj));
        }
      }
    }
  }
}

//this finds components of z-contained pixels, where the components 
//aren't allowed to go through w-contained pixels.  
//because it is safer to assume that things *are* connected, 
//we *allow* these components to go through pixels which are 
//touched by w (including going diagonally when one side is contained 
//in w and one side is touched by w)
void TrapGrid::pursue_z_cut_by_w_comp(int i, int j, int ind) {
  std::deque<Point2d<int> > stack(1, Point2d<int>(i,j));
  grid[i][j].z_cut_by_w_comp = -2;
  int ii,jj;
  while (stack.size() > 0) {
    Point2d<int> p = stack.back();
    stack.pop_back();
    //mark p as done
    grid[p.x][p.y].z_cut_by_w_comp = ind;
    z_cut_by_w_components[ind].push_back(p);
    //figure out which sides run against the edge of the pixels
    bool right_max, top_max, left_min, bottom_min; 
    if (p.x == num_pixels - 1) {
      right_max = true;
      left_min = false;
    } else if (p.x == 0) {
      left_min = true;
      right_max = false;
    } else {
      left_min = right_max = false;
    }
    if (p.y == num_pixels - 1) {
      top_max = true;
      bottom_min = false;
    } else if (p.y == 0) {
      top_max = false;
      bottom_min = true;
    } else {
      top_max = bottom_min = false;
    }
    //right column
    if (!right_max) {
      ii = p.x+1;
      //bottom
      if (!bottom_min) {
        jj=p.y-1;
        if (grid[ii][jj].z_ball_status == 2
            && grid[ii][jj].z_cut_by_w_comp == -1
            && grid[ii][jj].w_ball_status < 2 
            && (grid[ii-1][jj].w_ball_status < 2 || grid[ii][jj+1].w_ball_status < 2)) {
          grid[ii][jj].z_cut_by_w_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }      
      //middle   
      jj = p.y;
      if (grid[ii][jj].z_ball_status == 2 
          && grid[ii][jj].w_ball_status < 2 
          && grid[ii][jj].z_cut_by_w_comp == -1) { 
        grid[ii][jj].z_cut_by_w_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
      //right top
      if (!top_max) {
        jj = p.y+1;
        if (grid[ii][jj].z_ball_status == 2
            && grid[ii][jj].z_cut_by_w_comp == -1
            && grid[ii][jj].w_ball_status < 2 
            && (grid[ii-1][jj].w_ball_status < 2 || grid[ii][jj-1].w_ball_status < 2)) {
          grid[ii][jj].z_cut_by_w_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
    }
    //top spot
    if (!top_max) {
      ii = p.x;
      jj = p.y + 1;
      if (grid[ii][jj].z_ball_status == 2
          && grid[ii][jj].z_cut_by_w_comp == -1
          && grid[ii][jj].w_ball_status < 2) {
        grid[ii][jj].z_cut_by_w_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
    }
    //left column
    if (!left_min) {
      ii = p.x-1;
      //top
      if (!top_max) {
        jj = p.y+1;
        if (grid[ii][jj].z_ball_status == 2
            && grid[ii][jj].z_cut_by_w_comp == -1
            && grid[ii][jj].w_ball_status < 2 
            && (grid[ii+1][jj].w_ball_status < 2 || grid[ii][jj-1].w_ball_status < 2)) {
          grid[ii][jj].z_cut_by_w_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
      //middle
      jj = p.y;
      if (grid[ii][jj].z_ball_status == 2 
          && grid[ii][jj].w_ball_status < 2 
          && grid[ii][jj].z_cut_by_w_comp == -1) { 
        grid[ii][jj].z_cut_by_w_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
      //bottom
      if (!bottom_min) {
        jj=p.y-1;
        if (grid[ii][jj].z_ball_status == 2
            && grid[ii][jj].z_cut_by_w_comp == -1
            && grid[ii][jj].w_ball_status < 2 
            && (grid[ii+1][jj].w_ball_status < 2 || grid[ii][jj+1].w_ball_status < 2)) {
          grid[ii][jj].z_cut_by_w_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
    }
    //bottom middle
    if (!bottom_min) {
      ii = p.x;
      jj = p.y - 1;
      if (grid[ii][jj].z_ball_status == 2
          && grid[ii][jj].z_cut_by_w_comp == -1
          && grid[ii][jj].w_ball_status < 2) {
        grid[ii][jj].z_cut_by_w_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
    }
     
  }//<-- end of while
}

void TrapGrid::pursue_w_cut_by_z_comp(int i, int j, int ind) {
  std::deque<Point2d<int> > stack(1, Point2d<int>(i,j));
  grid[i][j].w_cut_by_z_comp = -2;
  int ii,jj;
  while (stack.size() > 0) {
    Point2d<int> p = stack.back();
    stack.pop_back();
    //mark p as done
    grid[p.x][p.y].w_cut_by_z_comp = ind;
    w_cut_by_z_components[ind].push_back(p);
    //figure out which sides run against the edge of the pixels
    bool right_max, top_max, left_min, bottom_min; 
    if (p.x == num_pixels - 1) {
      right_max = true;
      left_min = false;
    } else if (p.x == 0) {
      left_min = true;
      right_max = false;
    } else {
      left_min = right_max = false;
    }
    if (p.y == num_pixels - 1) {
      top_max = true;
      bottom_min = false;
    } else if (p.y == 0) {
      top_max = false;
      bottom_min = true;
    } else {
      top_max = bottom_min = false;
    }
    //right column
    if (!right_max) {
      ii = p.x+1;
      //bottom
      if (!bottom_min) {
        jj=p.y-1;
        if (grid[ii][jj].w_ball_status == 2
            && grid[ii][jj].w_cut_by_z_comp == -1
            && grid[ii][jj].z_ball_status < 2 
            && (grid[ii-1][jj].z_ball_status < 2 || grid[ii][jj+1].z_ball_status < 2)) {
          grid[ii][jj].w_cut_by_z_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }      
      //middle   
      jj = p.y;
      if (grid[ii][jj].w_ball_status == 2 
          && grid[ii][jj].z_ball_status < 2 
          && grid[ii][jj].w_cut_by_z_comp == -1) { 
        grid[ii][jj].w_cut_by_z_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
      //right top
      if (!top_max) {
        jj = p.y+1;
        if (grid[ii][jj].w_ball_status == 2
            && grid[ii][jj].w_cut_by_z_comp == -1
            && grid[ii][jj].z_ball_status < 2 
            && (grid[ii-1][jj].z_ball_status < 2 || grid[ii][jj-1].z_ball_status < 2)) {
          grid[ii][jj].w_cut_by_z_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
    }
    //top spot
    if (!top_max) {
      ii = p.x;
      jj = p.y + 1;
      if (grid[ii][jj].w_ball_status == 2
          && grid[ii][jj].w_cut_by_z_comp == -1
          && grid[ii][jj].z_ball_status < 2) {
        grid[ii][jj].w_cut_by_z_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
    }
    //left column
    if (!left_min) {
      ii = p.x-1;
      //top
      if (!top_max) {
        jj = p.y+1;
        if (grid[ii][jj].w_ball_status == 2
            && grid[ii][jj].w_cut_by_z_comp == -1
            && grid[ii][jj].z_ball_status < 2 
            && (grid[ii+1][jj].z_ball_status < 2 || grid[ii][jj-1].z_ball_status < 2)) {
          grid[ii][jj].w_cut_by_z_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
      //middle
      jj = p.y;
      if (grid[ii][jj].w_ball_status == 2 
          && grid[ii][jj].z_ball_status < 2 
          && grid[ii][jj].w_cut_by_z_comp == -1) { 
        grid[ii][jj].w_cut_by_z_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
      //bottom
      if (!bottom_min) {
        jj=p.y-1;
        if (grid[ii][jj].w_ball_status == 2
            && grid[ii][jj].w_cut_by_z_comp == -1
            && grid[ii][jj].z_ball_status < 2 
            && (grid[ii+1][jj].z_ball_status < 2 || grid[ii][jj+1].z_ball_status < 2)) {
          grid[ii][jj].w_cut_by_z_comp = -2;
          stack.push_front(Point2d<int>(ii,jj));
        }
      }
    }
    //bottom middle
    if (!bottom_min) {
      ii = p.x;
      jj = p.y - 1;
      if (grid[ii][jj].w_ball_status == 2
          && grid[ii][jj].w_cut_by_z_comp == -1
          && grid[ii][jj].z_ball_status < 2) {
        grid[ii][jj].w_cut_by_z_comp = -2;
        stack.push_front(Point2d<int>(ii,jj));
      }
    }
     
  }//<-- end of while
}

//I think it's enough to not go diagonally
void TrapGrid::pursue_intersection_comp(int i, int j, int ind) {
  std::deque<Point2d<int> > stack(1, Point2d<int>(i,j));
  grid[i][j].i_comp = -2;
  while (stack.size() > 0) {
    Point2d<int> p = stack.back();
    stack.pop_back();
    //mark p as done
    grid[p.x][p.y].i_comp = ind;
    intersection_components[ind].push_back(p);
    //check the non-diagonals
    int ii = p.x + 1;
    int jj = p.y;
    if (ii < num_pixels 
        && grid[ii][jj].i_comp == -1
        && grid[ii][jj].w_ball_status >0 
        && grid[ii][jj].z_ball_status >0) {
      grid[ii][jj].i_comp = -2;
      stack.push_front(Point2d<int>(ii,jj));
    }
    ii = p.x-1;
    if (ii > 0
        && grid[ii][jj].i_comp == -1
        && grid[ii][jj].w_ball_status >0 
        && grid[ii][jj].z_ball_status >0) {
      grid[ii][jj].i_comp = -2;
      stack.push_front(Point2d<int>(ii,jj));
    }    
    ii = p.x;
    jj = p.y+1;
    if (jj < num_pixels
        && grid[ii][jj].i_comp == -1
        && grid[ii][jj].w_ball_status >0 
        && grid[ii][jj].z_ball_status >0) {
      grid[ii][jj].i_comp = -2;
      stack.push_front(Point2d<int>(ii,jj));
    }
    jj = p.y-1;
    if (jj > 0
        && grid[ii][jj].i_comp == -1
        && grid[ii][jj].w_ball_status >0 
        && grid[ii][jj].z_ball_status >0) {
      grid[ii][jj].i_comp = -2;
      stack.push_front(Point2d<int>(ii,jj));
    }    
  }
}


//for each intersection boundary, follow it around and get an ordered 
//list of what components appear (currently this only finds outside boundaries)
//I think it's enough to check non-diagonal pixels.
//This scans through pixels until it finds (i,j) which is *not* an intersection, 
//and (i,j+1) which *is*.  I think this must always exist for any boundary.
//Then it walks around
void TrapGrid::compute_intersection_boundaries() {
  intersection_boundaries.resize(intersection_components.size(), std::vector<Point3d<int> >(0));
  std::vector<bool> icomp_done(intersection_components.size(), false);
  for (int i=0; i<num_pixels; ++i) {
    for (int j=1; j<num_pixels; ++j) {
      if (grid[i][j-1].i_comp == -1 &&
          grid[i][j].i_comp >= 0 &&
          icomp_done[grid[i][j].i_comp] == false) {
        pursue_intersection_boundary(i,j,0,intersection_boundaries[grid[i][j].i_comp]);
        icomp_done[grid[i][j].i_comp] = true;
      }
    }
  }
  
}

//pursue a boundary component by following it around
//the indices around a pixel start 0=bottom and go ccw around
void TrapGrid::pursue_intersection_boundary(int i, int j, int ind, std::vector<Point3d<int> >& bd) {
  Point2d<int> start_pixel(i,j);
  int start_ind = ind;
  Point2d<int> current_pixel = start_pixel;
  int current_ind = start_ind;
  int i_offset[4] = {0,1,0,-1};
  int j_offset[4] = {-1,0,1,0};
  bd.resize(0);
  do {
    int ii = current_pixel.x + i_offset[current_ind];
    int jj = current_pixel.y + j_offset[current_ind];
    if (ii < 0 || ii >= num_pixels || jj < 0 || jj >= num_pixels) {
      current_ind = (current_ind+1)%4;
      continue;
    }
    int zcwc = grid[ii][jj].z_cut_by_w_comp;
    int wczc = grid[ii][jj].w_cut_by_z_comp;
    if (zcwc >= 0) {
      Point3d<int> this_pair(0, grid[ii][jj].z_comp, zcwc);
      if (bd.size() == 0 || 
         (bd.size() == 1 && bd[0] != this_pair) ||
         (bd.size() >1 && bd[bd.size()-1] != this_pair && bd[0] != this_pair)) {
        bd.push_back(this_pair);
      }
    } else if (wczc >= 0) {
      Point3d<int> this_pair(1, grid[ii][jj].w_comp, wczc);
      if (bd.size() == 0 || 
         (bd.size() == 1 && bd[0] != this_pair) ||
         (bd.size() >1 && bd[bd.size()-1] != this_pair && bd[0] != this_pair)) {
        bd.push_back(this_pair);
      }
    }
    if (grid[ii][jj].i_comp == -1) {
      current_ind = (current_ind+1)%4;
    } else {
      current_pixel = Point2d<int>(ii,jj);
      current_ind = (current_ind+3)%4; //i.e. -1 mod 4 to stay on the boundary
    }
    
  } while (current_pixel != start_pixel || current_ind != start_ind);
}

//this finds a list of all interleaved components
//the result is a 4-tuple of Point3d<int> triples as above such that 
//t[0], t[2] are the same z component but different z cut by w components
//t[1], t[3] are the same w component but different w cut by z components
//algorithm: for every intersection component and for every component we 
//encounter around the boundary, go around again and see if we see it again.
//if so, divide up the boundary into two pieces and iterate over all possibilities
//there has to be a better way.
//
//the good_components list the components that this function is allowed to use
bool TrapGrid::find_interleaved_components(std::vector<std::vector<Point3d<int> > >& interleaved_components,
                                   const std::vector<bool>& good_components_z,
                                   const std::vector<bool>& good_components_w) {
  interleaved_components.resize(0);
  std::vector<Point3d<int> > temp_ic(4);
  for (int b=0; b<(int)intersection_boundaries.size(); ++b) {
    std::vector<Point3d<int> >& ib = intersection_boundaries[b];
    int bL = (int)ib.size();
    for (int j0=0; j0<bL; ++j0) {
      if (ib[j0].x != 0 || !good_components_z[ib[j0].z]) continue;
      for (int j1=j0+1; j1<bL; ++j1) {
        if (ib[j1].x != 0 || !good_components_z[ib[j1].z]) continue;
        if (ib[j0].x == ib[j1].x && 
            ib[j0].y == ib[j1].y && 
            ib[j0].z != ib[j1].z) {
          for (int k0=j0+1; k0<j1; ++k0) {
            if (ib[k0].x != 1 || !good_components_w[ib[k0].z]) continue;
            for (int k1p=0; k1p<(j0-j1+bL); ++k1p) {
              int k1 = (j1+k1p)%bL;
              if (ib[k1].x != 1 || !good_components_w[ib[k1].z]) continue;
              if (ib[k0].x == ib[k1].x && 
                  ib[k0].y == ib[k1].y && 
                  ib[k0].z != ib[k1].z) {
                temp_ic[0] = ib[j0];
                temp_ic[1] = ib[k0];
                temp_ic[2] = ib[j1];
                temp_ic[3] = ib[k1];
                interleaved_components.push_back(temp_ic);
                return true;
              }
            }
          }
        }
      }
    }
  }
  return (interleaved_components.size() > 0);
}


//populate the z_distance and w_distance on the pixels
//here it is distance to a pixel which is just touched by w or z
//this sets the appropriate pixels to have distance 0 and then
//does a breadth-first search (it's important that it's breadth-first)
void TrapGrid::compute_distances() {
  std::deque<Point2d<int> > z_stack(0);
  std::deque<Point2d<int> > w_stack(0);
  //clear the distances
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      grid[i][j].z_distance = grid[i][j].w_distance = -1;
    }
  }
  //initialize the ones of length 0
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      if (grid[i][j].z_ball_status > 0) {
        grid[i][j].z_distance = 0;
        z_stack.push_front(Point2d<int>(i,j));
      }
      if (grid[i][j].w_ball_status > 0) {
        grid[i][j].w_distance = 0;
        w_stack.push_front(Point2d<int>(i,j));
      }
    }
  }
  //now do the stacks; if we encounter anything new, set its length
  //to be 1+ours and push it on the stack
  while (z_stack.size() > 0) {
    Point2d<int> p = z_stack.back();
    z_stack.pop_back();
    if (p.y > 0 && grid[p.x][p.y-1].z_distance == -1) {
      grid[p.x][p.y-1].z_distance = grid[p.x][p.y].z_distance + 1;
      z_stack.push_front(Point2d<int>(p.x, p.y-1));
    }
    if (p.x < num_pixels-1 && grid[p.x+1][p.y].z_distance == -1) {
      grid[p.x+1][p.y].z_distance = grid[p.x][p.y].z_distance + 1;
      z_stack.push_front(Point2d<int>(p.x+1, p.y));
    }
    if (p.y < num_pixels-1 && grid[p.x][p.y+1].z_distance == -1) {
      grid[p.x][p.y+1].z_distance = grid[p.x][p.y].z_distance + 1;
      z_stack.push_front(Point2d<int>(p.x, p.y+1));
    }
    if (p.x > 0 && grid[p.x-1][p.y].z_distance == -1) {
      grid[p.x-1][p.y].z_distance = grid[p.x][p.y].z_distance + 1;
      z_stack.push_front(Point2d<int>(p.x-1, p.y));
    }
  }
  while (w_stack.size() > 0) {
    Point2d<int> p = w_stack.back();
    w_stack.pop_back();
    if (p.y > 0 && grid[p.x][p.y-1].w_distance == -1) {
      grid[p.x][p.y-1].w_distance = grid[p.x][p.y].w_distance + 1;
      w_stack.push_front(Point2d<int>(p.x, p.y-1));
    }
    if (p.x < num_pixels-1 && grid[p.x+1][p.y].w_distance == -1) {
      grid[p.x+1][p.y].w_distance = grid[p.x][p.y].w_distance + 1;
      w_stack.push_front(Point2d<int>(p.x+1, p.y));
    }
    if (p.y < num_pixels-1 && grid[p.x][p.y+1].w_distance == -1) {
      grid[p.x][p.y+1].w_distance = grid[p.x][p.y].w_distance + 1;
      w_stack.push_front(Point2d<int>(p.x, p.y+1));
    }
    if (p.x > 0 && grid[p.x-1][p.y].w_distance == -1) {
      grid[p.x-1][p.y].w_distance = grid[p.x][p.y].w_distance + 1;
      w_stack.push_front(Point2d<int>(p.x-1, p.y));
    }
  }
        
    
}

//compute_distances must have already been called.  This finds 
//the point in a z_cut_by_w component (or vice versa) which 
//is farthest from the other component
Point2d<int> TrapGrid::farthest_from_other_component(int z_or_w, int cut_comp) {
  std::vector<Point2d<int> >& c = (z_or_w == 0 ? z_cut_by_w_components[cut_comp] 
                                               : w_cut_by_z_components[cut_comp]);
  int max_dist = 0;
  int max_ind = -1;
  for (int i=0; i<(int)c.size(); ++i) {
    int d = (z_or_w==0 ? grid[c[i].x][c[i].y].w_distance : grid[c[i].x][c[i].y].z_distance);
    if (max_ind == -1 || d > max_dist) {
      max_ind = i;
      max_dist = d;
    }
  }
  return c[max_ind];
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

//return the lower left and upper right pixels 
//of the smallest rectangle containing all the pixels in the list
void TrapGrid::compute_pixel_extents(const std::vector<Point2d<int> >& L, 
                                     Point2d<int>& ll, 
                                     Point2d<int>& ur) {
  ll = ur = L[0];
  for (int i=0; i<(int)L.size(); ++i) {
    if (L[i].x < ll.x) {
      ll.x = L[i].x;
    }
    if (L[i].y < ll.y) {
      ll.y = L[i].y;
    }
    if (L[i].x > ur.x) {
      ur.x = L[i].x;
    }
    if (L[i].y > ur.y) {
      ur.y = L[i].y;
    }
  }
}



void TrapGrid::show_connected_components() {
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
  std::vector<int> zc(z_components.size());
  std::vector<int> wc(w_components.size());
  std::vector<int> ic(intersection_components.size());
  
  for (int i=0; i<(int)z_components.size(); ++i) {
    zc[i] = X2.get_rgb_color(0,0,0.5 + 0.5*double(i)/(double(z_components.size())));
  }
  for (int i=0; i<(int)w_components.size(); ++i) {
    wc[i] = X2.get_rgb_color(0,0.5 + 0.5*double(i)/(double(w_components.size())),0);
  }
  for (int i=0; i<(int)intersection_components.size(); ++i) {
    ic[i] = X2.get_rgb_color(0.5 + 0.5*double(i)/(double(intersection_components.size())),0,0);
  }  
  int wh = X2.get_rgb_color(1,1,1);
  int col;
  for (int i=0; i<num_pixels; ++i) {
    p.x = pixel_group_width*i;
    for (int j=0; j<num_pixels; ++j) {
      p.y = pixel_group_width*j;
      col = wh;
      if (grid[i][j].z_comp != -1) {
        col = zc[grid[i][j].z_comp];
      }
      if (grid[i][j].w_comp != -1) {
        col = wc[grid[i][j].w_comp];
      }
      if (grid[i][j].i_comp != -1) {
        col = ic[grid[i][j].i_comp];
      }
      if (num_pixels<512) {
        X2.draw_box(p, pixel_group_width, col);
      } else {
        X2.draw_point(p, col);
      }
    }
  }
  
  XGraphics X3(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
  
  std::vector<int> zcwc(z_cut_by_w_components.size());
  std::vector<int> wczc(w_cut_by_z_components.size());
   for (int i=0; i<(int)z_cut_by_w_components.size(); ++i) {
    zcwc[i] = X3.get_rgb_color(0,0,0.5 + 0.5*double(i)/(double(z_cut_by_w_components.size())));
  }
  for (int i=0; i<(int)w_cut_by_z_components.size(); ++i) {
    wczc[i] = X3.get_rgb_color(0,0.5 + 0.5*double(i)/(double(w_cut_by_z_components.size())),0);
  }
  for (int i=0; i<num_pixels; ++i) {
    p.x = pixel_group_width*i;
    for (int j=0; j<num_pixels; ++j) {
      p.y = pixel_group_width*j;
      col = wh;
      if (grid[i][j].z_cut_by_w_comp != -1) {
        col = zcwc[grid[i][j].z_cut_by_w_comp];
      }
      if (grid[i][j].w_cut_by_z_comp != -1) {
        col = wczc[grid[i][j].w_cut_by_z_comp];
      }
      if (grid[i][j].i_comp != -1) {
        col = ic[grid[i][j].i_comp];
      }
      if (num_pixels<512) {
        X3.draw_box(p, pixel_group_width, col);
      } else {
        X3.draw_point(p, col);
      }
    }
  }
  (void)X3.wait_for_key();
}


void TrapGrid::show(std::vector<Point2d<int> >* marked_points,
                    std::vector<Ball>* b) {
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
  int bc = X2.get_rgb_color(0,0,0);
  if (b != NULL) {
    //now draw the balls
    double real_pixel_width = box_width/num_drawing_pixels;
    int bc = X2.get_rgb_color(0,0,0);
    for (int i=0; i<(int)b->size(); ++i) {
      int x = int( ((*b)[i].center.real()-lower_left.real())/real_pixel_width );
      int y = int( ((*b)[i].center.imag()-lower_left.imag())/real_pixel_width );
      double draw_radius = (*b)[i].radius / real_pixel_width;
      //std::cout << "Drawing disk of radius " << (*b)[i].radius << " at " << (*b)[i].center << "\n";
      //std::cout << "In the drawing, at " << Point2d<int>(x,y) << " radius: " << draw_radius << "\n";
      X2.draw_disk(Point2d<int>(x,y), draw_radius, bc);
    }
  }
  if (marked_points != NULL) {
    for (int i=0; i<(int)marked_points->size(); ++i) {
      if (pixel_group_width > 1) {
        Point2d<int> p((*marked_points)[i].x*pixel_group_width, 
                       (*marked_points)[i].y*pixel_group_width);
        X2.draw_box(p, pixel_group_width,  bc);
      } else {
        X2.draw_point((*marked_points)[i],  bc);
      }
    }
  }
  (void)X2.wait_for_key();
}

void TrapGrid::show_distance_functions() {
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
  int z_max_dist = 0;
  int w_max_dist = 0;
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      if (grid[i][j].z_distance > z_max_dist) z_max_dist = grid[i][j].z_distance;
      if (grid[i][j].w_distance > w_max_dist) w_max_dist = grid[i][j].w_distance;
    }
  }
  std::vector<int> z_cols(z_max_dist+1);
  for (int i=0; i<z_max_dist; ++i) {
    z_cols[i] = X2.get_rgb_color(0,0,double(i)/double(z_max_dist));
  }
  std::vector<int> w_cols(w_max_dist+1);
  for (int i=0; i<w_max_dist; ++i) {
    w_cols[i] = X2.get_rgb_color(0,double(i)/double(w_max_dist),0);
  }
  for (int i=0; i<num_pixels; ++i) {
    p.x = pixel_group_width*i;
    for (int j=0; j<num_pixels; ++j) {
      p.y = pixel_group_width*j;
      if (pixel_group_width > 1) {
        X2.draw_box(p, pixel_group_width, z_cols[grid[i][j].z_distance]);
      } else {
        X2.draw_point(p, z_cols[grid[i][j].z_distance]);
      }
    }
  }
  (void)X2.wait_for_key();
  X2.erase_field();
  for (int i=0; i<num_pixels; ++i) {
    p.x = pixel_group_width*i;
    for (int j=0; j<num_pixels; ++j) {
      p.y = pixel_group_width*j;
      if (pixel_group_width > 1) {
        X2.draw_box(p, pixel_group_width, w_cols[grid[i][j].w_distance]);
      } else {
        X2.draw_point(p, w_cols[grid[i][j].w_distance]);
      }
    }
  }
  (void)X2.wait_for_key();
}
















