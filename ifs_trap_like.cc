#include <algorithm>

#include "ifs.h"

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



void show_stuff(const std::vector<cpx>* points,
                const std::vector<cpx>* red_points,
                const std::vector<cpx>* blue_points,
                const std::vector<cpx>* connected_points,
                const std::vector<cpx>* connected_points_2) {
  cpx ll,ur;
  std::vector<cpx> all_points;
  if (points != NULL) {
    all_points.insert(all_points.end(), points->begin(), points->end());
  }
  if (red_points != NULL) {
    all_points.insert(all_points.end(), red_points->begin(), red_points->end());
  }
  if (blue_points != NULL) {
    all_points.insert(all_points.end(), blue_points->begin(), blue_points->end());
  }
  box_containing_points(all_points, ll, ur);
  double drawing_width = ur.real() - ll.real();
  int num_drawing_pixels = 512;
  double pixel_diameter = drawing_width / num_drawing_pixels;
  XGraphics X2(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
  int bcol = X2.get_rgb_color(0,0,0);
  if (points != NULL) {
    for (int i=0; i<(int)points->size(); ++i) {
      Point2d<int> p( ((*points)[i].real() - ll.real())/pixel_diameter, 
                      ((*points)[i].imag() - ll.imag())/pixel_diameter );
      X2.draw_dot( p, bcol);
    }
  }
  if (red_points != NULL) {
    int rcol = X2.get_rgb_color(1,0.2,0);
    for (int i=0; i<(int)red_points->size(); ++i) {
      Point2d<int> p( ((*red_points)[i].real() - ll.real())/pixel_diameter, 
                      ((*red_points)[i].imag() - ll.imag())/pixel_diameter );
      X2.draw_dot( p, rcol);
    }
  } 
  if (blue_points != NULL) {
    int blcol = X2.get_rgb_color(0,0.8,1);
    for (int i=0; i<(int)blue_points->size(); ++i) {
      Point2d<int> p( ((*blue_points)[i].real() - ll.real())/pixel_diameter, 
                      ((*blue_points)[i].imag() - ll.imag())/pixel_diameter );
      X2.draw_dot( p, blcol);
    }
  } 
  if (connected_points != NULL) {
    for (int i=0; i<(int)connected_points->size(); ++i) {
      int ip1 = (i == (int)connected_points->size()-1 ? 0 : i+1);
      Point2d<int> p1( ((*connected_points)[i].real() - ll.real())/pixel_diameter,
                       ((*connected_points)[i].imag() - ll.imag())/pixel_diameter );
      Point2d<int> p2( ((*connected_points)[ip1].real() - ll.real())/pixel_diameter,
                       ((*connected_points)[ip1].imag() - ll.imag())/pixel_diameter );
      X2.draw_line(p1, p2, bcol);
    }
  }  
  if (connected_points_2 != NULL) {
    for (int i=0; i<(int)connected_points_2->size(); ++i) {
      int ip1 = (i == (int)connected_points_2->size()-1 ? 0 : i+1);
      Point2d<int> p1( ((*connected_points_2)[i].real() - ll.real())/pixel_diameter,
                       ((*connected_points_2)[i].imag() - ll.imag())/pixel_diameter );
      Point2d<int> p2( ((*connected_points_2)[ip1].real() - ll.real())/pixel_diameter,
                       ((*connected_points_2)[ip1].imag() - ll.imag())/pixel_diameter );
      X2.draw_line(p1, p2, bcol);
    }
  }
  X2.wait_for_key();
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
    if (UR.strictly_contains(X[i]) && UL.strictly_contains(X[i]) &&
        LL.strictly_contains(X[i]) && LR.strictly_contains(X[i])) {
      continue;
    } else {
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
  
  //std::vector<cpx> rp(indices_1.size());
  //for (int i=0; i<(int)rp.size(); ++i) {
  //  rp[i] = X[indices_1[i]];
  //}
  //std::vector<cpx> bp(indices_2.size());
  //for (int i=0; i<(int)bp.size(); ++i) {
  //  bp[i] = X[indices_2[i]];
  //}
  //show_stuff(NULL, &rp, &bp, NULL, NULL);
  
  //the sub convex hulls
  std::vector<int> ch_1;
  std::vector<int> ch_2;
  convex_hull_recurse(ch_1, indices_1, X);
  convex_hull_recurse(ch_2, indices_2, X);
  
  //std::vector<cpx> h_1(ch_1.size());
  //for (int i=0; i<(int)ch_1.size(); ++i) {
  //  h_1[i] = X[ch_1[i]];
  //}  
  //std::vector<cpx> h_2(ch_2.size());
  //for (int i=0; i<(int)ch_2.size(); ++i) {
  //  h_2[i] = X[ch_2[i]];
  //}
  //show_stuff(NULL, &rp, &bp, &h_1, &h_2);
  
  
  //recombine -- there are two halfspaces to find
  //they go top_2 - > top_1 and bottom_1 -> bottom_2
  int top_1;
  int top_2;
  int bottom_1;
  int bottom_2;
  extreme_indices(top_2, top_1, ch_2, ch_1, X);
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
  
  //std::vector<cpx> h(ch.size());
  //for (int i=0; i<(int)ch.size(); ++i) {
  //  h[i] = X[ch[i]];
  //}
  //show_stuff(NULL, &rp, &bp, &h, NULL);
  
}


void convex_hull(std::vector<int>& ch, 
                 const std::vector<cpx>& X) {

  //std::cout << "Got input points: \n";
  //show_stuff(&X, NULL, NULL, NULL, NULL);
  
  //this sets ch to list the indices that might be in the convex hull
  heuristic_convex_hull(ch, X);
  
  //std::cout << "Did the heuristic: \n";
  //std::vector<cpx> temp_X(ch.size());
  //for (int i=0; i<(int)ch.size(); ++i) {
  //  temp_X[i] = X[ch[i]];
  //}
  //show_stuff(&temp_X, NULL, NULL, NULL, NULL);

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
                                     int num_ball_trials,
                                     const std::vector<Ball>& balls,
                                     double radius_subtraction,
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
  std::vector<Ball> TLB_untranslated(0);
  while (true) {
    if ((int)TLB_untranslated.size() >= num_TL_balls || current_gap_ind >= (int)ch.size()) {
      break;
    }
    int i = ch_gap_pairs[current_gap_ind].second;
    cpx x1 = boundary_points[2*i+1];
    cpx x2 = boundary_points[2*(i==(int)ch.size()-1 ? 0 : i+1)];
    cpx v = -perp_to(x2-x1); //p should point towards the other disks
    v = v/abs(v);
    cpx best_center = 0;
    double best_radius = -1;
    double step = 1.0/double(num_ball_trials+1);
    for (double j=step; j<1.0; j+=step) {
      cpx p = j*x1 + (1-j)*x2;
      double t = when_ray_hits_ball(p,v,balls);
      t/=2.0;
      double d = distance_from_balls(p+t*v,balls);
      if (d > best_radius) {
        best_radius = d;
        best_center = p+t*v;
      }
    }
    //subtract off what we need to
    double alpha = best_radius - radius_subtraction;
    if (verbose>0) {
      std::cout << "Found ball of radius " << best_radius << 
                   " which gives a TLB of radius " << alpha << "\n";
    }
    //std::cout << "Found best ball " << best_center << " " << best_radius << "\n";
    if (alpha > 0) {
      TLB_untranslated.push_back(Ball(best_center, alpha));
      TLB.push_back(Ball(best_center-x1, alpha));
      TLB.push_back(Ball(best_center-x2, alpha));
    }
    ++current_gap_ind;
  }
  if (verbose > 0) {
    cpx ll, ur;
    box_containing_balls(balls, ll, ur);  
    double drawing_width = ur.real() - ll.real();
    int num_drawing_pixels = 512;
    double pixel_diameter = drawing_width / double(num_drawing_pixels);
    XGraphics X2(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
    int bcol = X2.get_rgb_color(0,0,0);
    int blue_color = X2.get_rgb_color(0, 0.8, 1);
    int red_color = X2.get_rgb_color(1, 0.2, 0);
    //draw the convex hull
    for (int i=0; i<(int)ch.size(); ++i) {
      int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
      cpx p1 = boundary_points[2*i+1];
      cpx p2 = boundary_points[2*ip1];
      Point2d<int> p1i( (p1.real() - ll.real())/pixel_diameter, 
                        (p1.imag() - ll.imag())/pixel_diameter );
      Point2d<int> p2i( (p2.real() - ll.real())/pixel_diameter, 
                        (p2.imag() - ll.imag())/pixel_diameter );
      X2.draw_line(p1i, p2i, bcol);
    }
    //draw all the disks
    for (int i=0; i<(int)balls.size(); ++i) {
      Point2d<int> p( (balls[i].center.real() - ll.real())/pixel_diameter, 
                      (balls[i].center.imag() - ll.imag())/pixel_diameter );
      double r = balls[i].radius/pixel_diameter;
      X2.draw_disk(p, r, blue_color);
    }
    //draw the points on the convex hull
    for (int i=0; i<(int)boundary_points.size(); ++i) {
      Point2d<int> p( (boundary_points[i].real() - ll.real())/pixel_diameter, 
                      (boundary_points[i].imag() - ll.imag())/pixel_diameter );
      X2.draw_dot(p, bcol);
    }
    //draw the trap like balls
    for (int i=0; i<(int)TLB_untranslated.size(); ++i) {
      Point2d<int> p( (TLB_untranslated[i].center.real() - ll.real())/pixel_diameter, 
                      (TLB_untranslated[i].center.imag() - ll.imag())/pixel_diameter );
      double r = TLB_untranslated[i].radius/pixel_diameter;
      X2.draw_disk(p, r, red_color);
      X2.draw_dot(p, bcol);
    }
    (void)X2.wait_for_key();
    if (TLB.size() > 0) {
      box_containing_balls(TLB, ll, ur);
      drawing_width = ur.real() - ll.real();
      num_drawing_pixels = 512;
      pixel_diameter = drawing_width / double(num_drawing_pixels);
      XGraphics X3(num_drawing_pixels, num_drawing_pixels, 1, Point2d<float>(0,0));
      for (int i=0; i<(int)TLB.size(); ++i) {
        Point2d<int> p( (TLB[i].center.real() - ll.real())/pixel_diameter, 
                        (TLB[i].center.imag() - ll.imag())/pixel_diameter );
        double r = TLB[i].radius/pixel_diameter;
        X3.draw_disk(p, r, red_color);
        X3.draw_dot(p, bcol);
      }
      X3.wait_for_key();
    }
  }
}
    
  
  
    

bool ifs::trap_like_balls(std::vector<Ball>& TLB, 
                          double initial_radius, 
                          double radius_subtraction,
                          int n_depth,
                          int verbose) {
  
  int old_depth = depth;
  depth = n_depth;
  if (!circ_connected(initial_radius)) {
    if (verbose>0) {
      std::cout << "Not even connected\n";
    }
    depth = old_depth;
    return false;
  }
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,initial_radius);
  std::vector<Ball> balls(0);
  compute_balls_right(balls, initial_ball, n_depth);
  
  if (verbose>0) {
    std::cout << "Computed all balls\n";
  }
  
  trap_like_balls_from_balls(TLB, 5, 3, balls, radius_subtraction, verbose);
 
  if (verbose>0) {
    std::cout << "Found trap like balls\n";
  }
 
  depth = old_depth;
  return true;
  
}




void ifs::get_TLB_constants(cpx ll, cpx ur, double& K, double& C, double& A, double& Z) {
  cpx z0 = 0.5*(ur+ll);
  double dr = ur.real() - z0.real();
  double di = ur.imag() - z0.imag();
  K=C=A=0;
  Z = 5;
  for (int i=0; i<4; ++i) {
    cpx z( z0.real() + ((i&1) == 1 ? -1 : 1)*dr, 
           z0.imag() + (((i>>1)&1) == 1 ? -1 : 1)*di );
    double a = abs(z)/abs(z0);
    if (a > A) A = a;
    double c = 1.0/((abs(z)-1)*(abs(z)-1));
    if (c > C) C = c;
    double k = abs(z-1.0)/(2.0*(1.0-abs(z)));
    if (k > K) K = k;
    double ZZ = abs(z);
    if (ZZ < Z) Z = ZZ;
  }  
}

//for the current z value, produce a bunch of promising uv words, 
//plus some trap like balls for this value
bool ifs::TLB_for_region(std::vector<Ball>& TLB, 
                        cpx ll, cpx ur, int n_depth, double* TLB_C, double* TLB_Z, int verbose) {
  
  cpx backup_z = z;
  cpx backup_w = w;
  z = (ll+ur)/2.0;
  az = abs(z);
  w = z; aw = az;
  
  //see the paper for a description of the constants
  double K,C,A,Z;
  get_TLB_constants(ll,ur,K,C,A,Z);
  if (TLB_C != NULL) *TLB_C = C;
  if (TLB_Z != NULL) *TLB_Z = Z;
  double dr = 0.5*(ur.real() - ll.real());
  double di = 0.5*(ur.real() - ll.real());
  double d = (di > dr ? di : dr);
  cpx z0 = 0.5*(ll+ur);
  double Rz0 = pow(A, (double)n_depth)*K + 4*K + 3.0*(1.0/pow(abs(z0), n_depth))*C*sqrt(2)*d;
  double radius_subtraction = pow(abs(z0), n_depth)*Rz0 + 2*C*sqrt(2)*d;
  
  if (Rz0 > 100000) return false;
  
  if (verbose>0) {
    std::cout << "Finding TLB in box " << ll << " " << ur << " to depth " << n_depth << "\n";
    std::cout << "A=" << A << " K=" << K << " C=" << C << "\n";
    std::cout << "Rz0: " << pow(A, (double)n_depth)*K << " + " 
              << 4*K << " + " << 3.0*(1.0/pow(abs(z0), n_depth))*C*sqrt(2)*d << " = " << Rz0 << "\n";
  }
  
  double min_r;
  if (!minimal_enclosing_radius(min_r)) {
    if (verbose>0) {
      std::cout << "No minimal radius\n";
    }
    return false;
  }
  if (verbose>0) std::cout << "Minimal radius: " << min_r << "\n";
  int old_depth = depth;
  depth = n_depth;
  if (!circ_connected(min_r)) {
    if (verbose>0) {
      std::cout << "Not connected at\n";
    }
    depth = old_depth;
    return false;
  }
  depth = old_depth;
  
  if (!trap_like_balls(TLB, Rz0, radius_subtraction, n_depth, verbose)) {
    z = backup_z; az = abs(z);
    w = backup_w; aw = abs(w);
    return false;
  }
  
  //now we need to find uv words 
  //z^m(u(1/2)-v(1/2)) needs to land within Cz*box_diag_rad to be feasible for this box
  //std::cout << "Finding close uv words to depth " << uv_depth << "\n";
  //find_close_uv_words(words, TLB, Cz*box_diag_rad, 10000, uv_depth);
  //words.resize(0);
  if (verbose>0) {
    std::cout << "Box radius: " << d << "\n";
    std::cout << "Rz0: " << Rz0 << "\n";
    std::cout << "Found the TLB: \n";
    for (int i=0; i<(int)TLB.size(); ++i) {
      std::cout << TLB[i] << "\n";
    }
    //std::cout << "Found the uv words: \n";
    //for (int i=0; i<(int)words.size(); ++i) {
    //  std::cout << i << ": " << words[i].first << "\n" << words[i].second << "\n";
    //}
  }
  z = backup_z; az = abs(z);
  w = backup_w; aw = abs(w);
  return true;
}



int ifs::check_TLB_and_uv_words(const std::vector<Ball>& TLB, 
                                const std::vector<std::pair<Bitword,Bitword> >& words) {
  for (int i=0; i<(int)words.size(); ++i) {
    cpx u12 = apply_bitword(words[i].first, 0.5);
    cpx v12 = apply_bitword(words[i].second, 0.5);
    cpx zm = pow(z, -words[i].first.len);
    cpx x = zm*(u12-v12);
    for (int j=0; j<(int)TLB.size(); ++j) {
      if (abs(TLB[j].center - x) < TLB[j].radius) {
        return words[i].first.len;
      }
    }
  }
  return -1;
}


int ifs::check_TLB(const std::vector<Ball>& TLB, 
                   double* TLB_C, double* TLB_Z, 
                   double& trap_radius, 
                   std::vector<std::pair<Bitword,Bitword> >* trap_w, 
                   int uv_depth) {
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return -1;
  Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
  std::deque<std::pair<Ball, Ball> > stack(1);
  stack[0] = std::make_pair(act_on_right(0,b), act_on_right(1,b));
  if (stack[0].first.is_disjoint(stack[0].second)) return -1;
  
  int current_looking_depth = stack[0].first.word_len;
  bool found_one = false;
  double best_radius = 0;
  if (trap_w != NULL) trap_w->resize(0);
  while (stack.size() > 0) {
    Ball bz = stack.back().first;
    Ball bw = stack.back().second;
    stack.pop_back();
    
    if (current_looking_depth < bz.word_len) {
      if (found_one) {
        trap_radius = best_radius;
        return current_looking_depth;
      } else {
        current_looking_depth = bz.word_len;
      }
    }
    
    //we are assuming they are not disjoint if they got pushed on
    //so check the displacement vector
    cpx d = bz.center - bw.center;
    d *= pow(z, -bz.word_len);
    for (int i=0; i<(int)TLB.size(); ++i) {
      double dist = abs(TLB[i].center - d) - TLB[i].radius;
      if (dist < -0.001) {
        if (TLB_C == NULL) return bz.word_len;
        double ep = (pow(*TLB_Z, bz.word_len)/(2.0*(*TLB_C)))*(TLB[i].radius - abs(TLB[i].center - d));
        if (ep > best_radius) {
          best_radius = ep;
        }
        if (trap_w != NULL) {
          if (ep > best_radius) {
            trap_w->insert(trap_w->begin(), std::make_pair(Bitword(bz.word, bz.word_len),
                                                           Bitword(bw.word, bw.word_len)));
          } else {
            trap_w->push_back(std::make_pair(Bitword(bz.word, bz.word_len),
                                             Bitword(bw.word, bw.word_len)));
          }
        }
        found_one = true;
      }
    }
    
    //if the word length is too big, we can't push the children
    if (bz.word_len >= uv_depth) continue;
    
    //now push on any children
    //if they are not disjoint, put them on the stack
    Ball bzs[2] = {act_on_right(0, bz), act_on_right(1, bz)};
    Ball bws[2] = {act_on_right(0, bw), act_on_right(1, bw)};
    for (int i=0; i<4; ++i) {
      if ( !bzs[i>>1].is_disjoint(bws[i&1]) ) { 
        stack.push_front(std::make_pair(bzs[i>>1], bws[i&1]));
      }
    }
  }
  return -1;
  
}


int ifs::check_limit_TLB(const std::vector<Ball>& TLB, 
                         double* TLB_C, double* TLB_Z, 
                         double& trap_radius, 
                         std::vector<std::pair<Bitword,Bitword> >* trap_w,
                         int n_limit) {
  //first, find the trap vectors
  std::vector<std::pair<Bitword,Bitword> > trap_pairs(0);
  int diff = check_TLB(TLB, TLB_C, TLB_Z, trap_radius, &trap_pairs, n_limit);
  if (diff < 0) return -1;
  
  int n;
  int trap_word_len;
  int old_n = 0;
  if (trap_w != NULL) trap_w->resize(0);
  bool found_one = false;
  
  cpx old_z = z;
  cpx old_w = w;
  cpx omega = cpx(0.3718586800741364, 0.5194111537479428);
  cpx ep = old_z - omega;
  set_params(omega, omega);
  
  int verbose = 0;
  
  if (verbose>0) {
    std::cout << "Trap balls: \n";
    for (int i=0; i<(int)TLB.size(); ++i) {
      std::cout << TLB[i] << "\n";
    }
    std::cout << "Got " << trap_pairs.size() << " trap pairs\n";
  }
  
  double best_trap_radius = 0;
  
  //for each one, see whether it is a limit vector:
  for (int i=0; i<(int)trap_pairs.size(); ++i) {
    std::pair<Bitword,Bitword> p = trap_pairs[i];
    trap_word_len = p.first.len;
    
    if (verbose>0) std::cout << "Checking trap word pair: " << p.first << " " << p.second << "\n";
  
    //first, pick off the prefixes from the words
    if (p.first.prefix(8).str() != std::string("01000111") ||
        p.second.prefix(8).str() != std::string("10111000")) {
      continue;
    }
    
    //get the suffixes -- these are u and v
    n = trap_word_len-8;
    Bitword u = p.first.suffix(n);
    Bitword v = p.second.suffix(n);
    if (n > old_n && found_one ) break;
    
    if (verbose>0) std::cout << "Got the suffixes: " << u << "  " << v << "\n";
    
    //compute p_u(omega)-p_v(omega)
    cpx puzpvz = apply_bitword(u,0.5)-apply_bitword(v,0.5);
    
    //put it together into the sum
    cpx right_summand = pow(omega, -n)*( puzpvz + 1.0 );
    cpx left_summand = 2.0*ep*pow(omega, -n-8)*(1.0-2.0*omega+5.0*pow(omega,4)-8.0*pow(omega,7));
    cpx putative_vector = left_summand + right_summand;
    
    if (verbose>0) {
      std::cout << "n: " << n << "\n";
      std::cout << "u(1/2): " << apply_bitword(u,0.5) << "\nv(1/2): " << apply_bitword(u,0.5) << "\n";
      std::cout << "puzpvz: " << puzpvz << "\n";
      std::cout << "epsilon: " << ep << "\n";
      std::cout << "w^(-n-8): " << pow(omega, -n-8) << "\n";
      std::cout << "w^(-n): " << pow(omega, -n) << "\n";
      std::cout << "Left summand: " << left_summand << "\n";
      std::cout << "Right summand: " << right_summand << "\n";
      std::cout << "Putative vector: " << putative_vector << "\n";
    }
    
    //check the vector
    for (int j=0; j<(int)TLB.size(); ++j) {
      double dist = TLB[j].radius - abs(putative_vector - TLB[j].center);
      if (dist < 0) continue;
      if (dist > best_trap_radius) best_trap_radius = dist / abs(left_summand/ep);
      found_one = true;
      if (verbose>0) std::cout << "Vector " << putative_vector << " is inside " << TLB[j] << "\n";
      if (trap_w != NULL) {
        trap_w->push_back(p);
      } else {
        //goto BREAKALL;
      }
    }
    old_n = n;
  }
  
  BREAKALL:
  
  set_params(old_z, old_w);
  trap_radius = best_trap_radius;
  if (trap_radius < 1e-12) return -1;
  
  return (found_one ? n+8 : -1);
}


int ifs::check_limit_TLB_recursive(const std::vector<Ball>& TLB, 
                                   double* TLB_C, double* TLB_Z, 
                                   double& trap_radius, 
                                   std::vector<std::pair<Bitword,Bitword> >* trap_w, 
                                   int uv_depth) {
  
  cpx old_z = z;
  cpx old_w = w;
  cpx omega = cpx(0.3718586800741364, 0.5194111537479428);
  cpx ep = old_z - omega;
  set_params(omega, omega);
  int verbose = 0;
  
  
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return -1;
  bool found_one = false;
  double best_radius;
  int n=10;
  if (trap_w != NULL) trap_w->resize(0);
  for (n = 10; n < uv_depth; ++n) {
    Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
    std::deque<std::pair<Ball, Ball> > stack(1);
    stack[0] = std::make_pair(act_on_right(0,b), act_on_right(1,b));
    if (stack[0].first.is_disjoint(stack[0].second)) return -1;
  
    std::vector<Ball> target_list(TLB.size());
    cpx left_summand = 2.0*ep*pow(omega, -n-8)*(1.0-2.0*omega+5.0*pow(omega,4)-8.0*pow(omega,7));
    for (int i=0; i<(int)TLB.size(); ++i) {
      target_list[i] = Ball(pow(omega, n)*(TLB[i].center - left_summand) - 1.0,
                            abs(pow(omega, n))*TLB[i].radius);
    }    

    while (stack.size() > 0) {
      Ball bz = stack.back().first;
      if (bz.word_len == n) break;
      Ball bw = stack.back().second;
      stack.pop_back();
      
      //so check the displacement vector
      cpx d = bz.center - bw.center;
      bool potentially_good_pair = false;
      for (int i=0; i<(int)target_list.size(); ++i) {
        if (abs(d-target_list[i].center) < 2*bz.radius + 2*target_list[i].radius) {
          potentially_good_pair = true;
          break;
        }
      }
      //if there might be the possibility that it works, push on the children
      if (potentially_good_pair) {
        Ball bzs[2] = {act_on_right(0, bz), act_on_right(1, bz)};
        Ball bws[2] = {act_on_right(0, bw), act_on_right(1, bw)};
        for (int i=0; i<4; ++i) {
          stack.push_front(std::make_pair(bzs[i>>1], bws[i&1]));
        }
      }
    }
    
    //the remaining balls should be checked
    best_radius = 0;
    for (int i=0; i<(int)stack.size(); ++i) {
      Ball bz = stack[i].first;  Ball bw = stack[i].second;
      cpx puzpvz = bz.center - bw.center;
      cpx right_summand = pow(omega, -n)*( puzpvz + 1.0 );
      cpx putative_vector = left_summand + right_summand;
      for (int j=0; j<(int)TLB.size(); ++j) {
        double dist = TLB[j].radius - abs(putative_vector - TLB[j].center);
        if (dist < 0) continue;
        double this_radius = dist / abs(left_summand/ep);
        std::cout << "Found limit trap with radius " << this_radius << " compared to best radius " << best_radius << "\n";
        if (this_radius > best_radius) {
          best_radius = this_radius;
          if (trap_w != NULL) {
            trap_w->insert(trap_w->begin(), std::make_pair(Bitword(bz.word, bz.word_len),
                                                          Bitword(bw.word, bw.word_len)) );
          }
        } else if (trap_w != NULL) {
          trap_w->push_back(std::make_pair(Bitword(bz.word, bz.word_len),
                                           Bitword(bw.word, bw.word_len)) );
        }
        if (best_radius > 1e-13) {
          std::cout << "Our best radius is good enough\n";
          found_one = true;
        } else {
          std::cout << "Found one, but it's not good enough\n";
        }
      }
    }
    if (found_one) break;
  }
  
  set_params(old_z, old_w);
  
  if (found_one) {
    std::cout << "Returning radius " << best_radius << "\n";
    trap_radius = best_radius;
    return n+8;
  } else {
    std::cout << "Found no trap\n";
    return -1;
  }
  
}























bool ifs::find_TLB_along_loop(const std::vector<cpx>& loop, 
                              bool draw_it, 
                              int verbose) {
  
  //trap parameters
  //int uv_depth = depth;


  int nL = loop.size();
  if (verbose>0) {
    std::cout << "Finding traps along the loop:\n";
    for (int i=0; i<nL; ++i) {
      std::cout << i << ": " << loop[i] << "\n";
    }
  }
  if (verbose>0) std::cout << "Finding traps along the vertices:\n";
  //this is a list of the balls along each segment of the path
  std::vector<std::vector<std::pair<cpx,double> > > trap_list(loop.size());
  
  //get graphics stuff
  //int rcol = X.get_rgb_color(1,0,0);
  double pixel_width = (2*wind)/double(drawing_width);
  double difficulty;
  
  //find the TLB and stuff
  std::vector<Ball> TLB;
  cpx center = z;
  double TLB_Z, TLB_C;
  TLB_for_region(TLB, center-cpx(wind,wind), center+cpx(wind,wind),
                  15, &TLB_C, &TLB_Z, 0);
  
  if (verbose>0) {
    std::cout << "Got the TLB constants: " << TLB_C << " " << TLB_Z << "\n";
  }
  
  //get the traps at the vertices
  for (int i=0; i<nL; ++i) {
    trap_list[i].resize(1);
    z = loop[i]; az = abs(z);
    w = z; aw = az;
    double epsilon;
    if ( (difficulty = check_TLB(TLB,&TLB_C, &TLB_Z,epsilon,NULL,depth)) < 0 ) {
      if (verbose>0) std::cout << "Failed to find a trap at vertex " << i << "\n";
      return false;
    }
    trap_list[i][0] = std::make_pair(z, epsilon);
    if (draw_it) {
      Point2d<int> p = cpx_to_point_mandlebrot(z);
      double r = trap_list[i][0].second / pixel_width;
      if (r<2) r = 2;
      double gamount = difficulty/100.0;
      int c = X.get_rgb_color(0, gamount, 1.0);
      X.draw_disk(p,r,c);
      X.flush();
    }
    if (verbose>0) std::cout << i << ": " << trap_list[i][0].first << ", " << trap_list[i][0].second << "\n";
  }
  
  //for each interval, go along it, placing the center of the 
  //next trap at exactly the edge of the previous one
  if (verbose>0) std::cout << "Pixel width: " << pixel_width << "\n";
  for (int i=0; i<nL; ++i) {
    //vector of length 1 pointing along the path
    cpx d = trap_list[(i+1)%nL][0].first - trap_list[i][0].first;
    d = d / abs(d);
    while (true) {
      std::pair<cpx,double> last_trap = trap_list[i].back();
      //detect if we are done
      double d_to_end = abs(trap_list[(i+1)%nL][0].first - last_trap.first);
      if (d_to_end < last_trap.second ||
          d_to_end < trap_list[(i+1)%nL][0].second) {
        if (verbose>1) std::cout << "Done this edge\n";
        break;
      }
      //otherwise, go to the end of our current ball
      cpx new_z = last_trap.first + last_trap.second*d;
      //set it
      z = new_z; az = abs(z);
      w = z; aw = az;
      //run it
      trap_list[i].resize(trap_list[i].size()+1);
      trap_list[i].back().first = z;
      if ( (difficulty = check_TLB(TLB, &TLB_C, &TLB_Z, trap_list[i].back().second, NULL, depth)) < 0 ) {
        if (verbose>0) std::cout << "Failed to find trap at " << z << "\n";
        return false;
      }
      //display it
      if (draw_it) {
        Point2d<int> p = cpx_to_point_mandlebrot(z);
        double r = trap_list[i].back().second / pixel_width;
        if (r<2) r = 2;
        double gamount = difficulty/100.0;
        int c = X.get_rgb_color(0, gamount, 1.0);
        X.draw_disk(p,r,c);
        X.flush();
      }
      if (verbose>1) {
        std::cout << "Found new trap " << trap_list[i].back().first << ", " << trap_list[i].back().second << "\n";
      }
    }
  }
  if (verbose>0) {
    std::cout << "Loop certified:\n";
    for (int i=0; i<nL; ++i) {
      std::cout << i << ": " << loop[i] << "\n";
    }
  }
  return true;
    
  
}









  
  
