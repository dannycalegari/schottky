//Functions to check for a trap


#include "trap_grid.h"

struct Trap {
};


bool ifs::find_trap(int verbose) {

  //first we need to rule out some silly stuff
  if (abs(z) > 0.999 || abs(w) > 0.999) return false;


  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double min_initial_radius = minimal_enclosing_radius();
  
  
  //starting depth will always be 10?
  int current_depth = depth; //10;
  std::vector<Ball> balls;
  compute_balls(balls, Ball(0.5, 0.5, min_initial_radius), current_depth);
  
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

  if (TG.grid_error) {
    std::cout << "grid error ";
    std::cout << "z=" << z << "; w=" << w << "min_initial_radius=" << min_initial_radius << "\n";
    return false;
  }

  //std::cout << "Trying to find trap...";
  
  //the grid is initialized
  while (true) {
  
    //show what it looks like
    if (verbose>0) TG.show(NULL, NULL);
        
    //find the connected components
    TG.compute_connected_components();
    
    if (verbose>0) {
      std::cout << "number of z components: " << TG.z_components.size() << "\n";
      std::cout << "number of w components: " << TG.w_components.size() << "\n";
      std::cout << "number of z cut by w components: " << TG.z_cut_by_w_components.size() << "\n";
      std::cout << "number of w cut by z components: " << TG.w_cut_by_z_components.size() << "\n";
      std::cout << "number of intersection components: " << TG.intersection_components.size() << "\n";
    }
        
    //find the boundaries of the intersection components
    TG.compute_intersection_boundaries();
    
    if (verbose>0) {
      for (int i=0; i<(int)TG.intersection_components.size(); ++i) {
        std::cout << "Boundary of intersection component " << i << ":\n";
        for (int j=0; j<(int)TG.intersection_boundaries[i].size(); ++j) {
          std::cout << TG.intersection_boundaries[i][j] << ",";
        }
        std::cout << "\n";
      }
    
      TG.show_connected_components();
    }
    
    std::vector<std::vector<Point3d<int> > > ic;
    if (!TG.find_interleaved_components(ic)) {
      if (verbose>0) std::cout << "Didn't find interleaved components\n";
      //std::cout << "no\n";
      return false;
    }
    if (verbose>0) {
      std::cout << "Found interleaved components: ";
      for (int i=0; i<(int)ic.size(); ++i) {
        for (int j=0; j<4; ++j) {
          std::cout << ic[i][j] << " ";
        }
        std::cout << "\n";
      }
    }
    
    TG.compute_distances();
    
    if (verbose>0) {
      std::cout << "Found distance functions\n";
      if (verbose>1) TG.show_distance_functions();
    }
    
    bool found_all_balls = true;
    std::vector<Ball> good_balls(4);
    for (int i=0; i<(int)ic.size(); ++i) {
      if (verbose>0) std::cout << "Trying the interleaved components number " << i << "\n";
    
      std::vector<Point2d<int> > farthest_points(4);
      for (int j=0; j<4; ++j) {
        farthest_points[j] = TG.farthest_from_other_component(ic[i][j].x, ic[i][j].z);
        if (verbose>0) std::cout << "Farthest point in " << ic[i][j] << " is " << farthest_points[j] << "\n";
      }
    
      if (verbose>0) {
        TG.show(&farthest_points, NULL);
      }

      found_all_balls = true;
      for (int j=0; j<4; ++j) {      
        Point2d<int>& p = farthest_points[j];
        int zw = ic[i][j].x;
        //if the pixel is touching a ball, it can't be good
        if ((zw == 0 ? TG.grid[p.x][p.y].w_distance == 0
                    : TG.grid[p.x][p.y].z_distance == 0)) {
          found_all_balls = false;
          break;
        }
        Ball b = balls[(zw==0 ? TG.grid[p.x][p.y].closest_z_ball
                              : TG.grid[p.x][p.y].closest_w_ball)];
        Ball b2 = b;
        if (verbose>0) std::cout << "Starting with ball " << b << "\n";
        for (int k=0; k<4; ++k) {
          b = act_on_right(0,b);
          b2 = act_on_right(1,b2);
        }
        if (verbose>0) std::cout << "Refined to balls:\n" << b << "\n" << b2 << "\n";
        if (TG.disjoint_from_z_or_w(b, 1-zw)) {
          good_balls[j] = b;
        } else if (TG.disjoint_from_z_or_w(b2, 1-zw)) {
          good_balls[j] = b2;
        } else {
          found_all_balls = false;
          break;
        }
      }
      if (found_all_balls) break; //we're done already
    }
    if (found_all_balls) {
      if (verbose>0) {
        std::cout << "Found four disjoint balls!\n";
        TG.show(NULL, &good_balls);
      }
      //std::cout << "yes!\n";
      return true;
    }
    
    //otherwise, we need to zoom in and stuff
    //to do this, find the largest connected component and zoom in
    //until it fills the middle third of the screen
    //std::cout << "no\n";
    return false;
    
    
  }//<- end of trap searching
  
  //std::cout << "yes\n";
  return  true;

}






void ifs::draw_trap() {

}




