#include "trap_grid.h"

//this does dot product as if they were vectors
double cpx_dot(cpx a, cpx b) {
  return a.real()*b.real() + a.imag()*b.imag();
}

//return a complex number b such that {b,a} is an oriented basis for R^2
cpx perp_to(cpx a) {
  return cpx(a.imag(), -a.real());
}


//returns the halfspace which is on the left when walking x1->x2
halfspace halfspace_on_left(cpx x1, cpx x2) {
  return halfspace(perp_to(x2-x1), x1);
}


//this function finds the extremal x and y points and throws out any point contained 
//in the diamond they make
//the return value is a list of the indices which *could* be in the convex hull
void heuristic_convex_hull(std::vector<int>& ch, const std::vector<cpx>& X) {
  //find the four extremal points
  int min_x_ind=0;
  int max_x_ind=0; 
  int min_y_ind=0; 
  int max_y_ind=0;
  double max_x = X[0].real();
  double min_x = X[0].real();
  double max_y = X[0].imag(); 
  double min_y = X[0].imag();
  for (int i=1; i<(int)X.size(); ++i) {
    double x = X[i].real();
    double y = X[i].imag();
    if (x < min_x) {
      min_x = x;
      min_x_ind = i;
    } else if (x > max_x) {
      max_x = x;
      max_x_ind = i;
    }
    if (y < min_y) {
      min_y = y;
      min_y_ind = i;
    } else if (y > max_y) {
      max_y = y;
      max_y_ind = i;
    }
  }
  if (min_x_ind == max_x_ind || min_x_ind == min_y_ind || min_x_ind == max_y_ind ||
      max_x_ind == min_y_ind || max_x_ind == max_y_ind || 
      min_y_ind == max_y_ind) { 
    //just return all of them -- we don't want to think about this case
    ch.resize(X.size());
    for (int i=0; i<(int)X.size(); ++i) {
      ch[i] = i;
    }
  }
  halfspace UR( perp_to(X[max_y_ind] - X[max_x_ind]), X[max_x_ind] );
  halfspace UL( perp_to(X[min_x_ind] - X[max_y_ind]), X[max_y_ind] );
  halfspace LL( perp_to(X[min_y_ind] - X[min_x_ind]), X[min_y_ind] );
  halfspace LR( perp_to(X[max_x_ind] - X[min_y_ind]), X[min_y_ind] );
  ch.resize(0);
  for (int i=0; i<(int)X.size(); ++i) {
    if (UR.contains(X[i]) && UL.contains(X[i]) &&
        LL.contains(X[i]) && LR.contains(X[i])) {
      ch.push_back(i);
    }
  }
}
  

//return the indices into the lists ch_1 and ch_2
//whose hyperplane to the left contains the other points
void extreme_indices(int& ind_1, int& ind_2, 
                    const std::vector<int>& ch_1,
                    const std::vector<int>& ch_2,
                    const std::vector<cpx>& X) {
  ind_1 = ind_2 = 0;
  while (true) {
    halfspace H = halfspace_on_left(X[ch_1[ind_1]], X[ch_2[ind_2]]);
    //find the point most violating
    int largest_1_violator = -1;
    double largest_1_violation = 0;
    for (int i=0; i<(int)ch_1.size(); ++i) {
      double v = H.val(X[ch_1[i]]);
      if (v > 0 && (largest_1_violator == -1 || v > largest_1_violation)) {
        largest_1_violator = i;
        largest_1_violation = v;
      }
    }
    int largest_2_violator = -1;
    double largest_2_violation = 0;
    for (int i=0; i<(int)ch_2.size(); ++i) {
      double v = H.val(X[ch_2[i]]);
      if (v > 0 && (largest_2_violator == -1 || v > largest_2_violation)) {
        largest_2_violator = i;
        largest_2_violation = v;
      }
    }
    if (largest_1_violator == -1 && largest_2_violator == -1) {
      return;
    }
    if (largest_1_violation > largest_2_violation) {
      ind_1 = largest_1_violator;
    } else {
      ind_2 = largest_2_violator;
    }
  }
}

  
bool cmp_pairs(const std::pair<double, int>& a, 
               const std::pair<double, int>& b) {
  return a.first < b.first;
}

bool cmp_pairs_reverse(const std::pair<double, int>& a, 
                       const std::pair<double, int>& b) {
  return a.first > b.first;
}


//this *assumes* that the indices in initial_indices are 
//sorted in increasing x order
void convex_hull_recurse(std::vector<int>& ch, 
                         const std::vector<int>& initial_indices,
                         const std::vector<cpx>& X) {
  if (initial_indices.size() < 2) {
    ch = initial_indices;
    return;
  
  } else if (initial_indices.size() == 2) {
    ch = initial_indices;
    return;
  }
  
  //if we're here, we should divide it, recurse, and recombine
  int cut_index = initial_indices.size()/2;
  int left_over = initial_indices.size()-cut_index;
  std::vector<int> indices_1(cut_index);
  std::vector<int> indices_2(left_over);
  for (int i=0; i<cut_index; ++i) {
    indices_1[i] = initial_indices[i];
  }
  for (int i=0; i<left_over; ++i) {
    indices_2[i] = initial_indices[cut_index+i];
  }
  
  //the sub convex hulls
  std::vector<int> ch_1;
  std::vector<int> ch_2;
  convex_hull_recurse(ch_1, indices_1, X);
  convex_hull_recurse(ch_2, indices_2, X);
  
  //recombine -- there are two halfspaces to find
  //they go top_2 - > top_1 and bottom_1 -> bottom_2
  int top_1;
  int top_2;
  int bottom_1;
  int bottom_2;
  extreme_indices(top_2, top_1, ch_2, ch_2, X);
  extreme_indices(bottom_1, bottom_2, ch_1, ch_2, X);
  
  //we get the convex hull by tracing out everything until the junctions, 
  //then skipping over them
  ch.resize(0);
  if (top_1 == bottom_1) {
    ch.push_back(ch_1[top_1]);
  } else {
    int i = top_1;
    while (true) {
      ch.push_back(ch_1[i]);
      if (i == bottom_1) break;
      i = (i == (int)ch_1.size()-1 ? 0 : i+1);
    }
  }
  
  if (top_2 == bottom_2) {
    ch.push_back(ch_2[top_2]);
  } else {
    int i = bottom_2;
    while (true) {
      ch.push_back(ch_2[i]);
      if (i == top_2) break;
      i = (i == (int)ch_2.size()-1 ? 0 : i+1);
    }
  }
  
}


void convex_hull(std::vector<int>& ch, 
                 const std::vector<cpx>& X) {

  //this sets ch to list the indices that might be in the convex hull
  heuristic_convex_hull(ch, X);

  //sort the points into increasing x order
  std::vector<std::pair<double, int> > x_vals(ch.size());
  for (int i=0; i<(int)ch.size(); ++i) {
    x_vals[i] = std::make_pair(X[ch[i]].real(), ch[i]);
  }
  std::sort(x_vals.begin(), x_vals.end(), cmp_pairs);
  
  //put them in the list
  std::vector<int> initial_indices(ch.size());
  for (int i=0; i<(int)ch.size(); ++i) {
    initial_indices[i] = x_vals[i].second;
  }
  
  convex_hull_recurse(ch, initial_indices, X);
  
}
    

//the halfspace H[i] is the halfspace between balls i and i+1 in the convex hull
//boundary_points[2*i] and boundary_points[2*i+1] live on the ball ch[i]
void ball_convex_hull(std::vector<int>& ch,
                      std::vector<cpx>& boundary_points,
                      std::vector<halfspace>& H,
                      const std::vector<Ball>& balls) {
  std::vector<cpx> centers(balls.size());
  for (int i=0; i<(int)balls.size(); ++i) {
    centers[i] = balls[i].center;
  }
  convex_hull(ch, centers);
  
  //make the halfspaces
  if (ch.size() < 2) {
    H.resize(0);
    return;
  }
  H.resize(ch.size());
  boundary_points.resize(2*ch.size());
  for (int i=0; i<(int)ch.size(); ++i) {
    int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
    cpx x1 = balls[ch[i]].center;
    double r = balls[ch[i]].radius;
    cpx x2 = balls[ch[ip1]].center;
    cpx v = perp_to(x2-x1); //vector pointing out of x1
    v = r*(v / abs(v));  //now it exactly touches the point on the disk
    boundary_points[2*i+1] = x1+v;
    boundary_points[2*ip1] = x2+v; //this only works if the balls are the same size
    H[i] = halfspace_on_left(x1+v, x2+v);
  } 
}



//finds the num_TL_balls-many longest convex hulll sections, and 
//then finds the largest trap ball it can stick in there
void ifs::trap_like_balls_from_balls(std::vector<Ball>& TLB, 
                                     int num_TL_balls, 
                                     const std::vector<Ball>& balls,
                                     int verbose) {
  //get the convex hull of the balls
  std::vector<int> ch;
  std::vector<cpx> boundary_points;
  std::vector<halfspace> H;
  ball_convex_hull(ch, boundary_points, H, balls);
  
  //find the largest gaps
  std::vector<std::pair<double, int> > ch_gap_pairs(ch.size());
  for (int i=0; i<(int)ch.size(); ++i) {
    int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
    ch_gap_pairs[i] = std::make_pair(abs(boundary_points[2*ip1]-boundary_points[2*i+1]), i);
  }
  std::sort(ch_gap_pairs.begin(), ch_gap_pairs.end(), cmp_pairs_reverse);
  
  std::vector<int> ch_gaps(0);
  int current_gap_ind = 0;
  TLB.resize(0);
  while (true) {
    if ((int)TLB.size() >= num_TL_balls || current_gap_ind > (int)ch.size()) {
      break;
    }
    int i = ch_gap_pairs[current_gap_ind].second;
    cpx x1 = boundary_points[2*i+1];
    cpx x2 = boundary_points[2*(i==(int)ch.size()-1 ? 0 : i+1)];
    cpx v = -perp_to(x2-x1); //p should point towards the other disks
    v = v/abs(v);
    cpx p1 = 0.25*x1 + 0.75*x2;
    cpx p2 = 0.5*x1 + 0.5*x2;
    cpx p3 = 0.75*x1 + 0.25*x2;
    double t1 = when_ray_hits_ball(p1, v, balls);
    double t2 = when_ray_hits_ball(p2, v, balls);
    double t3 = when_ray_hits_ball(p3, v, balls);
    t1 /= 2.0;
    t2 /= 2.0;
    t3 /= 2.0;
    double d1 = distance_from_balls(p1 + t1*v, balls);
    double d2 = distance_from_balls(p2 + t2*v, balls);
    double d3 = distance_from_balls(p3 + t3*v, balls);
    cpx best_center;
    double best_radius;
    if (d1 > d2 && d1 > d3) {
      best_center = p1+t1*v;
      best_radius = d1;
    } else if (d2 > d3) {
      best_center = p2 + t2*v;
      best_radius = d2;
    } else {
      best_center = p3 + t3*v;
      best_radius = d3;
    }
    TLB.push_back(Ball(best_center, best_radius));
  }
  if (verbose > 0) {
    cpx ll, ur;
    box_containing_balls(balls, ll, ur);  
    int num_drawing_pixels = 512;
    XGraphics X2(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
    //draw the convex hull
    for (int i=0; i<(int)ch.size(); ++i) {
      int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
    }
    (void)X2.wait_for_key();
  }
}
    
  
  
    

bool ifs::trap_like_balls(std::vector<Ball>& TLB, int verbose) {
  
  int n_depth = depth;
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return false;
  
  if (!circ_connected(min_r)) {
    if (verbose>0) {
      std::cout << "Not even connected\n";
    }
    return false;
  }
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  compute_balls_right(balls, initial_ball, n_depth);
  
  trap_like_balls_from_balls(TLB, 10, balls, verbose);
 
  return true;
  
}
  













  
  
  
  
