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
  
  TG.show();
  
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
  
  std::vector<Point3d<int> > ic;
  if (TG.find_interleaved_components(ic)) {
    std::cout << "Found interleaved components: ";
    for (int i=0; i<4; ++i) {
      std::cout << ic[i] << " ";
    }
    std::cout << "\n";
  }
  
  TG.show_connected_components();

  
  
  return true;
}






void ifs::draw_trap() {

}




