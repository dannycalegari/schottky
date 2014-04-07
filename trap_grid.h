#ifndef __TRAP_GRID__
#define __TRAP_GRID__

#include <vector>

#include "ifs.h"

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
  //these are ordered lists of (original component, cut component) pairs
  //one for each boundary
  std::vector<std::vector<Point2d<int> > > intersection_boundaries;
  
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
  void compute_intersection_boundaries();
  void pursue_intersection_boundary(int i, int j, int ind, std::vector<Point2d<int> >& bd);
  void show();
  void show_connected_components();
};


#endif