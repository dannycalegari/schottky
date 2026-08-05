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



/* a debugging window: plot several point sets in different colors and wait for
 * a keypress.  Nothing calls it, but it is handy to drop into a computation. */
#ifndef IFS_NO_GRAPHICS
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
#endif /* IFS_NO_GRAPHICS */




  


  
bool cmp_pairs(const std::pair<double, int>& a, 
               const std::pair<double, int>& b) {
  return a.first < b.first;
}

bool cmp_pairs_reverse(const std::pair<double, int>& a, 
                       const std::pair<double, int>& b) {
  return a.first > b.first;
}


//signed area of triangle O,A,B (>0 iff O->A->B is a left turn / CCW)
static double hull_cross(const cpx& O, const cpx& A, const cpx& B) {
  return (A.real()-O.real())*(B.imag()-O.imag())
       - (A.imag()-O.imag())*(B.real()-O.real());
}

//comparator for sorting point indices by (real, then imag)
struct hull_index_less {
  const std::vector<cpx>* X;
  hull_index_less(const std::vector<cpx>& x) : X(&x) {}
  bool operator()(int a, int b) const {
    if ((*X)[a].real() != (*X)[b].real()) return (*X)[a].real() < (*X)[b].real();
    return (*X)[a].imag() < (*X)[b].imag();
  }
};

//Robust convex hull via Andrew's monotone chain.
//Returns indices of X on the hull in COUNTERCLOCKWISE order, which is REQUIRED
//by the single caller ball_convex_hull (it uses perp_to(x2-x1) as the OUTWARD
//normal, which only points outward for a CCW traversal).  Collinear interior
//points are dropped (cross <= 0 is popped); exact-duplicate coordinates are
//collapsed.  This replaced an earlier heuristic + divide-and-conquer hull pair
//whose merge step (an unguarded while(true)) failed to terminate on the
//near-collinear covers of the dimension-2 attractor that occur on the circle
//|s|=1/sqrt(2).
void convex_hull(std::vector<int>& ch,
                 const std::vector<cpx>& X) {
  int N = (int)X.size();
  ch.resize(0);
  if (N == 0) return;

  //sort indices by (real, then imag)
  std::vector<int> idx(N);
  for (int i=0; i<N; ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), hull_index_less(X));

  //collapse exact-duplicate coordinates (keep first)
  std::vector<int> pts;
  pts.reserve(N);
  for (int k=0; k<N; ++k) {
    int i = idx[k];
    if (!pts.empty() && X[pts.back()].real()==X[i].real()
                     && X[pts.back()].imag()==X[i].imag()) continue;
    pts.push_back(i);
  }
  int n = (int)pts.size();
  if (n < 3) { ch = pts; return; }  //0,1,2 distinct points: hull is all of them

  //monotone chain: lower hull (L->R) then upper hull (R->L) gives CCW
  std::vector<int> h(2*n);
  int m = 0;
  for (int k=0; k<n; ++k) {
    while (m>=2 && hull_cross(X[h[m-2]], X[h[m-1]], X[pts[k]]) <= 0.0) --m;
    h[m++] = pts[k];
  }
  int lower = m+1;
  for (int k=n-2; k>=0; --k) {
    while (m>=lower && hull_cross(X[h[m-2]], X[h[m-1]], X[pts[k]]) <= 0.0) --m;
    h[m++] = pts[k];
  }
  h.resize(m-1);  //drop the last point (== first)
  ch = h;
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
#ifndef IFS_NO_GRAPHICS
  //show the hull, the balls and the trap-like balls we just found
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
#endif /* IFS_NO_GRAPHICS */
}


  
    

//keep only the maximal balls: drop b when some kept ball contains it.  Same greedy
//pass as prune_tlb.py -- sort by radius descending, keep a ball unless an already
//kept one swallows it.  O(n^2) in the kept count, which is fine at ~2000 balls.
static void prune_to_maximal(std::vector<Ball>& B) {
  std::vector<std::pair<double,int> > order(B.size());
  for (int i=0; i<(int)B.size(); ++i) order[i] = std::make_pair(B[i].radius, i);
  std::sort(order.begin(), order.end());
  std::reverse(order.begin(), order.end());
  std::vector<Ball> kept;
  for (int k=0; k<(int)order.size(); ++k) {
    const Ball& b = B[order[k].second];
    bool swallowed = false;
    for (int j=0; j<(int)kept.size(); ++j) {
      //b lies inside kept[j] exactly when |c_b - c_j| + r_b <= r_j
      if (abs(b.center - kept[j].center) + b.radius <= kept[j].radius) { swallowed = true; break; }
    }
    if (!swallowed) kept.push_back(b);
  }
  B.swap(kept);
}

bool ifs::trap_like_balls_many(std::vector<Ball>& TLB,
                               cpx ll, cpx ur, int n_depth,
                               int ngaps, int ntrials, int nradial,
                               bool prune,
                               int verbose) {
  TLB.resize(0);
  //the region-valid radius and subtraction, exactly as TLB_for_region computes them:
  //every ball has to stay trap-like over the whole parameter box, so the subtraction
  //absorbs both the depth-n truncation and the box width.
  double KK, CC, AA, ZZ;
  get_TLB_constants(ll, ur, KK, CC, AA, ZZ);
  cpx z0 = 0.5*(ll+ur);
  double Rz0 = std::pow(AA,(double)n_depth)*KK + 4*KK
             + 3.0*(1.0/std::pow(abs(z0),n_depth))*CC*sqrt(2.0)*abs(0.5*(ur-ll).real());
  double radius_subtraction = std::pow(abs(z0),n_depth)*Rz0
                            + 2*CC*sqrt(2.0)*abs(0.5*(ur-ll).real());
  double min_r;
  if (!minimal_enclosing_radius(min_r)) {
    if (verbose>0) std::cout << "trap_like_balls_many: no minimal enclosing radius\n";
    return false;
  }
  Ball initial_ball(0.5, (z-1.0)/2.0, (1.0-w)/2.0, Rz0);
  std::vector<Ball> balls;
  compute_balls_right(balls, initial_ball, n_depth);

  std::vector<int> ch; std::vector<cpx> bp; std::vector<halfspace> H;
  ball_convex_hull(ch, bp, H, balls);
  if (ch.size() < 2) {
    if (verbose>0) std::cout << "trap_like_balls_many: degenerate hull\n";
    return false;
  }
  //visit the widest gaps first: those are where a trap-like ball has room to be large
  std::vector<std::pair<double,int> > gp(ch.size());
  for (int i=0; i<(int)ch.size(); ++i) {
    int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
    gp[i] = std::make_pair(abs(bp[2*ip1]-bp[2*i+1]), i);
  }
  std::sort(gp.begin(), gp.end());
  std::reverse(gp.begin(), gp.end());

  for (int gi=0; gi<(int)gp.size() && gi<ngaps; ++gi) {
    int i = gp[gi].second;
    int ip1 = (i==(int)ch.size()-1 ? 0 : i+1);
    cpx x1 = bp[2*i+1], x2 = bp[2*ip1];
    cpx v = -perp_to(x2-x1);
    v = v/abs(v);                          //unit inward normal to the gap
    for (int j=1; j<=ntrials; ++j) {
      double lam = (double)j/(double)(ntrials+1);
      cpx p = lam*x1 + (1.0-lam)*x2;
      double thit = when_ray_hits_ball(p, v, balls);
      for (int k=1; k<=nradial; ++k) {
        cpx q = p + (thit*(double)k/(double)(nradial+1))*v;
        double alpha = distance_from_balls(q, balls) - radius_subtraction;
        if (alpha > 0) {
          //the two hull contact points give the two translates, as in CKW Def 8.2.3
          TLB.push_back(Ball(q-x1, alpha));
          TLB.push_back(Ball(q-x2, alpha));
        }
      }
    }
  }
  if (verbose>0)
    std::cout << "trap_like_balls_many: " << TLB.size() << " balls from "
              << balls.size() << " at depth " << n_depth << ", subtraction "
              << radius_subtraction << "\n";
  if (prune) {
    int before = (int)TLB.size();
    prune_to_maximal(TLB);
    if (verbose>0)
      std::cout << "  pruned " << before << " -> " << TLB.size() << " maximal\n";
  }
  return TLB.size() > 0;
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
        //Test before updating best_radius, not after.  The comparison used to sit below
        //the assignment, so it was always false and the best (u,v) pair was never moved
        //to the front of trap_w -- callers that take trap_w->front() as "the" certifying
        //pair were getting whichever pair happened to be found first instead.
        bool is_best = (ep > best_radius);
        if (is_best) best_radius = ep;
        if (trap_w != NULL) {
          if (is_best) {
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


//Guided level-synchronous beam variant of check_TLB.  At each depth level we
//keep only the beam_width most-promising pairs (smallest distance of the
//renormalized displacement to the trap-like balls) and expand those.  A pair
//is accepted as a trap by the SAME strict-containment test (dist < -0.001) and
//its parameter-disk radius is the SAME eps = (Z^m/2C)(radius - |center - d|)
//as in check_TLB, so a reported trap is exactly as trustworthy as check_TLB's.
//Pruning only affects COMPLETENESS: if the beam discards the branch that led to
//the only trap we return -1 (uncertified), never a false positive.
int ifs::check_TLB_bestfirst(const std::vector<Ball>& TLB,
                             double* TLB_C, double* TLB_Z,
                             double& trap_radius,
                             std::vector<std::pair<Bitword,Bitword> >* trap_w,
                             int uv_depth, int beam_width) {
  if (trap_w != NULL) trap_w->resize(0);
  if (TLB.size() == 0) return -1;
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return -1;

  struct BFNode { Ball bz; Ball bw; double key; };

  Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
  Ball f0 = act_on_right(0, b);   //u starts with f
  Ball g0 = act_on_right(1, b);   //v starts with g
  if (f0.is_disjoint(g0)) return -1;

  std::vector<BFNode> frontier;
  {
    cpx d = f0.center - g0.center;  d *= pow(z, -f0.word_len);
    double key = 1e300;
    for (int i=0; i<(int)TLB.size(); ++i) {
      double dd = abs(TLB[i].center - d) - TLB[i].radius;
      if (dd < key) key = dd;
    }
    BFNode nd; nd.bz = f0; nd.bw = g0; nd.key = key;
    frontier.push_back(nd);
  }

  for (int depth = f0.word_len; depth <= uv_depth && !frontier.empty(); ++depth) {
    //1) is there a trap at this level? (all frontier nodes have word_len==depth)
    bool found = false;
    double best_eps = 0.0;
    for (int n=0; n<(int)frontier.size(); ++n) {
      if (frontier[n].key >= -0.001) continue;   //displacement outside every TLB ball
      cpx d = frontier[n].bz.center - frontier[n].bw.center;
      d *= pow(z, -frontier[n].bz.word_len);
      for (int i=0; i<(int)TLB.size(); ++i) {
        double dist = abs(TLB[i].center - d) - TLB[i].radius;
        if (dist < -0.001) {
          if (TLB_C == NULL) { trap_radius = 0.0; return frontier[n].bz.word_len; }
          double ep = (pow(*TLB_Z, frontier[n].bz.word_len)/(2.0*(*TLB_C)))
                      * (TLB[i].radius - abs(TLB[i].center - d));
          if (ep > best_eps) {
            best_eps = ep;
            if (trap_w != NULL) {
              trap_w->insert(trap_w->begin(),
                std::make_pair(Bitword(frontier[n].bz.word, frontier[n].bz.word_len),
                               Bitword(frontier[n].bw.word, frontier[n].bw.word_len)));
            }
          }
          found = true;
        }
      }
    }
    if (found) { trap_radius = best_eps; return depth; }

    //2) no trap here: expand to the next level, keeping the best beam_width pairs
    if (depth >= uv_depth) break;
    std::vector<BFNode> children;
    children.reserve(frontier.size()*4);
    for (int n=0; n<(int)frontier.size(); ++n) {
      Ball bzs[2] = { act_on_right(0, frontier[n].bz), act_on_right(1, frontier[n].bz) };
      Ball bws[2] = { act_on_right(0, frontier[n].bw), act_on_right(1, frontier[n].bw) };
      for (int k=0; k<4; ++k) {
        Ball cz = bzs[k>>1];
        Ball cw = bws[k&1];
        if (cz.is_disjoint(cw)) continue;
        cpx d = cz.center - cw.center;  d *= pow(z, -cz.word_len);
        double key = 1e300;
        for (int i=0; i<(int)TLB.size(); ++i) {
          double dd = abs(TLB[i].center - d) - TLB[i].radius;
          if (dd < key) key = dd;
        }
        BFNode nd; nd.bz = cz; nd.bw = cw; nd.key = key;
        children.push_back(nd);
      }
    }
    if (beam_width > 0 && (int)children.size() > beam_width) {
      std::nth_element(children.begin(), children.begin()+beam_width, children.end(),
                       [](const BFNode& a, const BFNode& b){ return a.key < b.key; });
      children.resize(beam_width);
    }
    frontier.swap(children);
  }
  return -1;
}


//MIXED-LENGTH scaffold (see header).  u,v advance independently up to |len diff|<=kmax;
//displacement renormalized by s^{-min}.  SCAFFOLD: tests same-length trap-like balls with
//no protrusion check, so kmax>=1 is exploratory only.  kmax=0 should match check_TLB.
int ifs::check_TLB_mixed(const std::vector<Ball>& TLB,
                         double* TLB_C, double* TLB_Z,
                         double& trap_radius,
                         std::vector<std::pair<Bitword,Bitword> >* trap_w,
                         int uv_depth, int kmax) {
  if (trap_w != NULL) trap_w->resize(0);
  if (TLB.size()==0) return -1;
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return -1;
  Ball b(0.5,(z-1.0)/2.0,(1.0-w)/2.0,1.01*min_r);
  std::deque<std::pair<Ball,Ball> > stack;
  stack.push_back(std::make_pair(act_on_right(0,b), act_on_right(1,b)));
  bool found=false; double best_eps=0.0; int best_len=-1;
  long budget=0; const long BUD=20000000;
  while (!stack.empty()) {
    if (++budget>BUD) break;
    Ball bz=stack.back().first, bw=stack.back().second; stack.pop_back();
    int a = (bz.word_len<bw.word_len ? bz.word_len : bw.word_len);   // min length
    cpx d = (bz.center - bw.center) * pow(z, -a);
    for (int i=0;i<(int)TLB.size();++i) {
      double dist = abs(TLB[i].center - d) - TLB[i].radius;
      if (dist < -0.001) {
        if (TLB_C==NULL) { trap_radius=0.0; return a; }
        double ep = (pow(*TLB_Z, a)/(2.0*(*TLB_C)))*(TLB[i].radius - abs(TLB[i].center - d));
        if (ep > best_eps) {
          best_eps=ep; best_len=a;
          if (trap_w!=NULL) trap_w->insert(trap_w->begin(),
            std::make_pair(Bitword(bz.word,bz.word_len), Bitword(bw.word,bw.word_len)));
        }
        found=true;
      }
    }
    //expand: extend-both, extend-u, extend-v (respect kmax and prune disjoint)
    Ball zc[2]={act_on_right(0,bz), act_on_right(1,bz)};
    Ball wc[2]={act_on_right(0,bw), act_on_right(1,bw)};
    std::pair<Ball,Ball> kids[8]; int nk=0;
    bool cz = bz.word_len<uv_depth, cw = bw.word_len<uv_depth;
    if (cz && cw) for(int gz=0;gz<2;++gz) for(int gw=0;gw<2;++gw) kids[nk++]=std::make_pair(zc[gz],wc[gw]);
    if (cz)       for(int gz=0;gz<2;++gz) kids[nk++]=std::make_pair(zc[gz],bw);
    if (cw)       for(int gw=0;gw<2;++gw) kids[nk++]=std::make_pair(bz,wc[gw]);
    for (int c=0;c<nk;++c) {
      int ld = kids[c].first.word_len - kids[c].second.word_len; if (ld<0) ld=-ld;
      if (ld>kmax) continue;
      if (kids[c].first.is_disjoint(kids[c].second)) continue;
      stack.push_back(kids[c]);
    }
  }
  if (found) { trap_radius=best_eps; return best_len; }
  return -1;
}


/* HEXAHOLE-SPECIFIC.  check_limit_TLB and check_limit_TLB_recursive implement the
   limit traps of CKW section 9 at the ONE renormalization point omega used there:
   both hard-code omega = 0.3718586800741364 + 0.5194111537479428i below and then
   call set_params(omega,omega), overwriting whatever parameter the caller had.  So
   they answer a question about omega no matter what they are asked, which is why
   nothing in the interactive program calls them any more.

   For limit traps at an ARBITRARY renormalization point -- which is what covering a
   fundamental domain of E_sigma needs -- use funddom.c, or its C API funddom_core.h
   (fd_core_from_lm / fd_solver_init / fd_level).  That path takes the renormalization
   data (sigma, a, b, Delta, P'(sigma)) as input instead of assuming it, and it is what
   produces the coverage figures.  These two functions are kept as the reference
   implementation of the section 9 argument in its original form. */
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
  
  //BREAKALL:
  
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
        //std::cout << "Found limit trap with radius " << this_radius << " compared to best radius " << best_radius << "\n";
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
          //std::cout << "Our best radius is good enough\n";
          found_one = true;
        } else {
          //std::cout << "Found one, but it's not good enough\n";
        }
      }
    }
    if (found_one) break;
  }
  
  set_params(old_z, old_w);
  
  if (found_one) {
    //std::cout << "Returning radius " << best_radius << "\n";
    trap_radius = best_radius;
    return n+8;
  } else {
    //std::cout << "Found no trap\n";
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
#ifndef IFS_NO_GRAPHICS
    if (draw_it) {
      Point2d<int> p = cpx_to_point_mandelbrot(z);
      double r = trap_list[i][0].second / pixel_width;
      if (r<2) r = 2;
      double gamount = difficulty/100.0;
      int c = X.get_rgb_color(0, gamount, 1.0);
      X.draw_disk(p,r,c);
      X.flush();
    }
#endif
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
#ifndef IFS_NO_GRAPHICS
      if (draw_it) {
        Point2d<int> p = cpx_to_point_mandelbrot(z);
        double r = trap_list[i].back().second / pixel_width;
        if (r<2) r = 2;
        double gamount = difficulty/100.0;
        int c = X.get_rgb_color(0, gamount, 1.0);
        X.draw_disk(p,r,c);
        X.flush();
      }
#endif
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









  
  
