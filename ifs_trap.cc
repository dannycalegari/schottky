//Functions to check for a trap


#include "trap_grid.h"

struct Trap {
};


bool ifs::find_trap() {
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
  
  TG.show(NULL, NULL);
  
  //find the connected components
  TG.compute_connected_components();
  
  std::cout << "number of z components: " << TG.z_components.size() << "\n";
  std::cout << "number of w components: " << TG.w_components.size() << "\n";
  std::cout << "number of z cut by w components: " << TG.z_cut_by_w_components.size() << "\n";
  std::cout << "number of w cut by z components: " << TG.w_cut_by_z_components.size() << "\n";
  std::cout << "number of intersection components: " << TG.intersection_components.size() << "\n";
  
  //find the boundaries of the intersection components
  TG.compute_intersection_boundaries();
  
  for (int i=0; i<(int)TG.intersection_components.size(); ++i) {
    std::cout << "Boundary of intersection component " << i << ":\n";
    for (int j=0; j<(int)TG.intersection_boundaries[i].size(); ++j) {
      std::cout << TG.intersection_boundaries[i][j] << ",";
    }
    std::cout << "\n";
  }
  
  TG.show_connected_components();
  
  std::vector<Point3d<int> > ic;
  if (!TG.find_interleaved_components(ic)) {
    std::cout << "Didn't find interleaved components\n";
    return false;
  }
  std::cout << "Found interleaved components: ";
  for (int i=0; i<4; ++i) {
    std::cout << ic[i] << " ";
  }
  std::cout << "\n";
  
  TG.compute_distances();
  
  std::cout << "Found distance functions\n";
  
  TG.show_distance_functions();
  
  std::vector<Point2d<int> > farthest_points(4);
  for (int i=0; i<4; ++i) {
    farthest_points[i] = TG.farthest_from_other_component(ic[i].x, ic[i].z);
    std::cout << "Farthest point in " << ic[i] << " is " << farthest_points[i] << "\n";
  }
  
  TG.show(&farthest_points, NULL);
  
  std::vector<Ball> good_balls(4);
  bool found_all_balls = true;
  for (int i=0; i<4; ++i) {
    Point2d<int>& p = farthest_points[i];
    if (ic[i].x == 0) {
      if (TG.grid[p.x][p.y].w_distance == 0) {
        found_all_balls = false;
        break;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_z_ball];
      std::cout << "Starting with ball " << b << "\n";
      int attempt = 0;
      while (true) {
        if (attempt > 5) {
          found_all_balls = false;
          break;
        }
        if (TG.disjoint_from_z_or_w(b, 1)) {
          good_balls[i] = b;
          break;
        }
        b = act_on_right(0, b);
        std::cout << "No good; trying next ball " << b << "\n";
        ++attempt;
      }
      if (!found_all_balls) break;
    
    } else { //p[0] == 1 (i.e. w component)
      if (TG.grid[p.x][p.y].z_distance == 0) {
        found_all_balls = false;
        break;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_w_ball];
      std::cout << "Starting with ball " << b << "\n";
      int attempt = 0;
      while (true) {
        if (attempt > 5) {
          found_all_balls = false;
          break;
        }
        if (TG.disjoint_from_z_or_w(b, 0)) {
          good_balls[i] = b;
          break;
        }
        b = act_on_right(1, b);
        std::cout << "No good; trying next ball " << b << "\n";
        ++attempt;
      }
      if (!found_all_balls) break;
    }
  }
    
  if (!found_all_balls) {
    std::cout << "Didn't find four disjoint balls\n";
    return false;
  }
  std::cout << "Found four disjoint balls!\n";
  TG.show(NULL, &good_balls);
  
  //at this point, the only thing we need to check is that 
  //there are two points epsilon apart
  
  
  
  return true;
}






void ifs::draw_trap() {

}




