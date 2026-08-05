//Functions to check for a trap

#include <set>
#include <deque>
#include <string>
#include <algorithm>

#include "trap_grid.h"

struct Trap {
};





bool ifs::find_trap_given_balls_old(const std::vector<Ball>& initial_balls, 
                                int max_refinements,
                                int max_pixels,
                                bool far_trap_points,
                                double* minimum_trap_distance,
                                int verbose) {
  
  
  int current_refinements = 0;
  cpx ll, ur;
  
  //copy the balls
  std::vector<Ball> balls = initial_balls;  
  
  //these record the trap
  std::vector<Ball> trap_balls(4);
  
  
  while (true) {
  
    //if there are no balls, give up
    if (balls.size() == 0) return false;
    
    //find the ball extents
    box_containing_balls(balls, ll, ur);  
    
    if (verbose>0) std::cout << "New box: " << ll << " " << ur << "\n";
    
    //find the average radius
    double av_radius = 0;
    for (int i=0; i<(int)balls.size(); ++i) {
      av_radius += balls[i].radius;
    }
    av_radius /= (double)balls.size();
    if (verbose>0) std::cout << "Average radius: " << av_radius << "\n";
    
    //initialize the trap
    TrapGrid TG;
    
    //the trap will cover all the balls, and it will have as many pixels 
    //as it needs so that the average ball contains about 4 pixels
    double desired_pixel_diameter = av_radius/2.5;
    int np = int( ((ur.real() - ll.real())/desired_pixel_diameter) + 1);
    if (verbose > 0) std::cout << "Desired pixel radius: " << desired_pixel_diameter << "\nDesired pixels: " << np << "\n";
    if (np > max_pixels) np = max_pixels;
    TG.reset_grid(ll, ur, np);
    
    //fill the trap grid
    TG.fill_pixels(balls); 
  
    //show what it looks like
    if (verbose>0) {
      TG.show(NULL, NULL, NULL, NULL, NULL);
      if (verbose>1) TG.show(NULL, NULL, &balls, NULL, NULL);
    }
        
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
    
      if (verbose>1) TG.show_connected_components();
    }
    
    //find the distance functions -- might as well do this now
    TG.compute_distances();
    if (verbose>0) {
      std::cout << "Found distance functions\n";
      if (verbose>1) TG.show_distance_functions();
    }
    
    //now figure out which components are good -- into which ones can we fit 
    //a ball?
    std::vector<int> good_components[2]; //this records -1 if not good or distance to the other if good
    std::vector<Ball> good_balls[2];
    good_components[0].resize(TG.z_cut_by_w_components.size());
    good_balls[0].resize(TG.z_cut_by_w_components.size());
    good_components[1].resize(TG.w_cut_by_z_components.size());
    good_balls[1].resize(TG.w_cut_by_z_components.size());
    for (int i=0; i<(int)TG.z_cut_by_w_components.size(); ++i) {
      Point2d<int> p = TG.farthest_from_other_component(0, i);
      if (TG.grid[p.x][p.y].w_distance == 0) {
        good_components[0][i] = -1;
        continue;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_z_ball];
      Ball b2 = b;
      for (int k=0; k<5; ++k) {
        b = act_on_right(0,b);
        b2 = act_on_right(1,b2);
      }
      if (TG.disjoint_from_z_or_w(b, 1)) {
        good_balls[0][i] = b;
        good_components[0][i] = TG.grid[p.x][p.y].w_distance;
      } else if (TG.disjoint_from_z_or_w(b2, 1)) {
        good_balls[0][i] = b2;
        good_components[0][i] = TG.grid[p.x][p.y].w_distance;
      } else {
        good_components[0][i] = -1;
      }
    }
    for (int i=0; i<(int)TG.w_cut_by_z_components.size(); ++i) {
      Point2d<int> p = TG.farthest_from_other_component(1, i);
      if (TG.grid[p.x][p.y].z_distance == 0) {
        good_components[1][i] = -1;
        continue;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_w_ball];
      Ball b2 = b;
      for (int k=0; k<5; ++k) {
        b = act_on_right(0,b);
        b2 = act_on_right(1,b2);
      }
      if (TG.disjoint_from_z_or_w(b, 0)) {
        good_balls[1][i] = b;
        good_components[1][i] = TG.grid[p.x][p.y].z_distance;
      } else if (TG.disjoint_from_z_or_w(b2, 0)) {
        good_balls[1][i] = b2;
        good_components[1][i] = TG.grid[p.x][p.y].z_distance;
      } else {
        good_components[1][i] = -1;
      }
    } 
    
    if (verbose>0) {
      std::cout << "Found good components:\n";
      for (int i=0; i<(int)good_components[0].size(); ++i) {
        if (good_components[0][i]>=0) {
          std::cout << "(0," << i << ") ";
        }
      }
      std::cout << "\n";
      for (int i=0; i<(int)good_components[1].size(); ++i) {
        if (good_components[1][i]>=0) {
          std::cout << "(1," << i << ") ";
        }
      }
      std::cout << "\n";
    }
    //now let's find interleaved components
    //if we find interleaved good components, we're done
    std::vector<std::vector<Point3d<int> > > ic;
    if (far_trap_points) {
      if (!TG.find_interleaved_components(ic, good_components[0], good_components[1],far_trap_points)) {
        if (verbose>0) std::cout << "Didn't find far interleaved components\n";
        if (!TG.find_interleaved_components(ic, good_components[0], good_components[1],false)) {
          if (verbose>0) std::cout << "Didn't find any interleaved components\n";
          goto ZOOM_OR_REFINE;
        }
      } else {
        if (verbose>0) std::cout << "Found far interleaved components\n";
      }
    } else if (!TG.find_interleaved_components(ic, good_components[0], good_components[1],false)) {
      if (verbose>0) std::cout << "Didn't find any interleaved components\n";
      goto ZOOM_OR_REFINE; //sorry
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
    
    if (ic.size() > 0) {
      trap_balls[0] = good_balls[ic[0][0].x][ic[0][0].z];
      trap_balls[1] = good_balls[ic[0][1].x][ic[0][1].z];
      trap_balls[2] = good_balls[ic[0][2].x][ic[0][2].z];
      trap_balls[3] = good_balls[ic[0][3].x][ic[0][3].z];
      if (verbose>0) {
        std::cout << "Found four disjoint balls!\n";
        TG.show(NULL, NULL, &trap_balls, NULL, NULL);
      }
      //find the minimum trap distance
      //this will take some time, but I guess it's worth it?
      //std::cout << "Computing minimum distance\n";
      *minimum_trap_distance = 100;
      for (int i=0; i<4; ++i) {
        cpx ci = trap_balls[i].center;
        double ri = trap_balls[i].radius;
        for (int j=0; j<(int)balls.size(); ++j) {
          //if the last gen index is *equal* to the component index, 
          //then we don't want to compare -- we want to compare to the other comp
          if (balls[j].last_gen_index() == ic[0][i].x) continue;
          double bd = abs(ci - balls[j].center)-ri;
          if (bd < *minimum_trap_distance) {
            //std::cout << "Trap ball at " << ci << " new closest ball at " << balls[j].center << "\n";
            *minimum_trap_distance = bd;
            //std::cout << "New minimum: " << *minimum_trap_distance << "\n";
          }
        }
      }
          
      return true;
    }
    
    //otherwise, we need to zoom in and stuff
    //to do this, find the largest connected component and zoom in
    //until it fills the middle third of the screen
  ZOOM_OR_REFINE:
    
    //maybe we have to give up
    if (current_refinements >= max_refinements) {
      if (verbose>0) std::cout << "Exceeded depth limit\n";
      break;
    } else if (TG.intersection_components.size() == 0) {
      if (verbose>0) std::cout << "No intersection to refine!\n";
      break;
    }
    //otherwise, we're going for it
    ++current_refinements;
  
    int biggest_component=-1;
    int biggest_component_size = 0;
    for (int i=0; i<(int)TG.intersection_components.size(); ++i) {
      if (biggest_component == -1 || (int)TG.intersection_components[i].size() > biggest_component_size) {
        biggest_component_size = TG.intersection_components[i].size();
        biggest_component = i;
      }
    }
    Point2d<int> comp_ll, comp_ur;
    TG.compute_pixel_extents(TG.intersection_components[biggest_component], comp_ll, comp_ur);
    cpx comp_ll_cpx(TG.lower_left.real() + comp_ll.x*TG.pixel_diameter,
                    TG.lower_left.imag() + comp_ll.y*TG.pixel_diameter);
    cpx comp_ur_cpx(TG.lower_left.real() + comp_ur.x*TG.pixel_diameter,
                    TG.lower_left.imag() + comp_ur.y*TG.pixel_diameter);
    cpx comp_center_cpx = (comp_ll_cpx + comp_ur_cpx)/2.0;
    if (verbose>0) {
      std::cout << "Computed intersection component extents " << comp_ll_cpx << " " << comp_ur_cpx << "\n";
    }
    double putative_width = comp_ur_cpx.real() - comp_ll_cpx.real();
    double putative_height = comp_ur_cpx.imag() - comp_ll_cpx.imag();
    double true_radius = (putative_width > putative_height ? putative_width/2.0 
                                                           : putative_height/2.0);
    //the box should be 6/4 of the size of the component
    true_radius *= 1.5;
    cpx new_ll(comp_center_cpx.real()-true_radius, comp_center_cpx.imag()-true_radius);
    cpx new_ur(comp_center_cpx.real()+true_radius, comp_center_cpx.imag()+true_radius);
    
    if (verbose>0) {
      std::cout << "New box: " << new_ll << "-" << new_ur << "\n";
      TG.show(NULL,NULL,NULL,&new_ll, &new_ur);
    }
    
    //might as well refine the depth -- compute the balls which fit into
    //this box
    refine_balls_into_box(balls, new_ll, new_ur);
    
    //when we loop around, these are the new balls
    
  }//<- end of trap searching
  
  return false;
  
}















bool ifs::find_trap_given_balls(const std::vector<Ball>& balls,
                                int max_pixels,
                                double* min_trap_distance,
                                int verbose) {
  //if there are no balls, give up
  if (balls.size() == 0) return false;
  
  //find the ball extents
  cpx ll, ur;
  box_containing_balls(balls, ll, ur);  
  if (verbose>0) std::cout << "Grid: " << ll << " " << ur << "\n";
  
  //find the average radius
  double av_radius = 0;
  for (int i=0; i<(int)balls.size(); ++i) {
    av_radius += balls[i].radius;
  }
  av_radius /= (double)balls.size();
  if (verbose>0) std::cout << "Average radius: " << av_radius << "\n";
  
  //initialize the trap
  TrapGrid TG;
  
  //the trap will cover all the balls, and it will have as many pixels 
  //as it needs so that each ball contains a reasonable number of pixels
  double desired_pixel_diameter = av_radius/2.5;
  int np = int( ((ur.real() - ll.real())/desired_pixel_diameter) + 1);
  if (verbose > 0) std::cout << "Desired pixel diameter: " << desired_pixel_diameter << "\nDesired pixels: " << np << "\n";
  if (np > max_pixels) np = max_pixels;
  TG.reset_grid(ll, ur, np);
  
  //fill the trap grid
  TG.fill_pixels(balls); 
  
  //find the distance functions
  TG.compute_distances();
  
  //compute the boundary
  //the boundary is a list of 3 tuples (i,j), (distance from w if z, -distance from z if w, and 0 if intersection)
  std::vector<Point3d<int> > boundary(0);
  TG.compute_boundary(boundary);
  
  //compute the good pixels (anything on the boundary is good, 
  //and anything completely contained inside only one of z or w is good 
  TG.compute_good_pixels(boundary);
  
  
  if (verbose>0) {
    std::cout << "Boundary: ";
    for (int i=0; i<(int)boundary.size(); ++i) {
      std::cout << boundary[i] << " ";
    }
    std::cout << "\n";
    TG.show(NULL, &boundary, NULL, NULL, NULL);
  }
  
  //prune the boundary and find trap balls
  std::vector<Ball> trap_balls(4);
  if (!TG.prune_boundary(*this, balls, boundary, trap_balls)) {
    if (verbose>0) {
      std::cout << "No interleaved components\n";
    }
    return false;
  }
  if (verbose>0) {
    std::cout << "Found interleaved components and trap balls:\n";
    for (int i=0; i<4; ++i) {
      std::cout << boundary[i] << " " << trap_balls[i] << "\n";
    }
    TG.show(NULL,NULL,&trap_balls,NULL,NULL);
  }
  
  //find how far the balls are from the other component
  //(but only if we care)
  if (min_trap_distance == NULL) return true;
  *min_trap_distance = 10000;
  for (int i=0; i<4; ++i) {
    cpx ci = trap_balls[i].center;
    double ri = trap_balls[i].radius;
    int other_gen = 1 - trap_balls[i].last_gen_index();
    for (int j=0; j<(int)balls.size(); ++j) {
      if (balls[j].last_gen_index() != other_gen) continue;
      double d = abs(ci - balls[j].center) - ri;
      if (d < *min_trap_distance) {
        //std::cout << "New minimum distance of " << d << " between balls " << trap_balls[i] << " and " << balls[j] << "\n";
        *min_trap_distance = d;
      }
    }
  }
  if (verbose>0) std::cout << "Found minimum trap distance of " << *min_trap_distance << "\n";
  return true;
} 






bool ifs::find_trap(int max_uv_depth, int max_n_depth, int max_pixels, double Cz, double* epsilon, double* difficulty, int verbose) {

  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double min_initial_radius;
  if (!minimal_enclosing_radius(min_initial_radius)) {
    //std::cout << "initial radius is infinite\n";
    return false;
  }
  //make sure it's connected at the minimal possible radius
  if (verbose>0) std::cout << "Checking connectedness with minimal initial radius of " << min_initial_radius << "\n";
  int old_depth = depth;
  depth = max_uv_depth+max_n_depth+2;  //this ensures there are points within pr-r = 2r-r = r
  if (!circ_connected(min_initial_radius)) {
    if (verbose>0) {
      std::cout << "Not even connected\n";
    }
    depth = old_depth;
    return false;
  }
  depth = old_depth;
  
  double ratio_goal = 0.05;
  double ratio_lower_limit = 0.01;
  double current_ratio = 4;
  double d = 6.0;
  double p = 2.0;
  int current_n_depth = 10; //I guess always start here?
  int current_uv_depth = -1;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_initial_radius*d);
  Ball zb, wb;
  
  while (true) {
  
    //find actions u and v which start with z and w such that 
    //find_aligned_images_with_distinct_first_letters(initial_ball, z_cm, w_cm, uv_depth, zb, wb);
    if (abs(z) > (1.0/sqrt(2))) {
      if (verbose>0) {
        std::cout << "|z| is larger than 1/sqrt(2), so we're done\n";
      }
      if (epsilon != NULL) *epsilon = 0.005;
      if (difficulty != NULL) *difficulty = 1.0;
      return true;
    }
    if (current_ratio > ratio_goal) {
      if (verbose>0) std::cout << "current ratio: " << current_ratio << " trying to get below " << ratio_goal << "\n";
      std::cout.flush();
      find_aligned_images_with_distinct_first_letters(initial_ball, 0, 0, max_uv_depth, zb, wb, ratio_goal, ratio_lower_limit);
    }
    current_ratio = abs(zb.center-wb.center)/zb.radius;
    current_uv_depth = zb.word_len;
    
    if (verbose>0) {
      std::cout << "Found initial balls\n";
      std::cout << "z-ball: " << zb << "\nw-ball: " << wb << "\n";
      std::cout << "Ratio: " << current_ratio << "\n";
    }
    
    //now act on the *right* a bunch
    std::vector<Ball> ZB;
    compute_balls_right(ZB, zb, current_n_depth);
    std::vector<Ball> WB;
    compute_balls_right(WB, wb, current_n_depth);
    
    //now concatenate the lists and use that (we swap to avoid unnecessary work)
    std::vector<Ball> balls(0);
    balls.swap(ZB);
    balls.insert(balls.end(), WB.begin(), WB.end());
    
    //this returns the minimum distance between the trap
    //points and the other *centers*
    double min_trap_dist = -1;
    
    bool got_trap = find_trap_given_balls(balls, max_pixels, NULL, (verbose>0 ? verbose-1 : verbose));
    if (verbose>0) {
      if (got_trap) {
        std::cout << "Got trap at " << z << " with min trap dist " << min_trap_dist << "\n";
      } else {
        std::cout << "Failed to find trap\n";
      }
    }
    
    if (got_trap) {
      
      double ball_radius = balls[0].radius;
      double r = pow(az, current_n_depth + current_uv_depth)*min_initial_radius;
      double e = (1.0/Cz)*0.25*r*(d-p);
      if (verbose>0) {
        std::cout << "ball radius: " << ball_radius << " r = " << r << " d = " << d << " p= " << p << " epsilon = " << e << "\n";
      }
      if (epsilon != NULL) *epsilon = e;
      if (difficulty != NULL) *difficulty = ((-log10(e))-3.0)*15;
      return true;
    }
    //if we didn't get it, think about what to do
    if (current_n_depth < max_n_depth) {
      current_n_depth += (current_n_depth < max_n_depth-1 ? 3 : 1);
      if (verbose>0) std::cout << "Increasing depth to " << current_n_depth << "\n";
    } else if (current_ratio < ratio_goal && current_ratio > 2*ratio_lower_limit) {
      ratio_goal = (ratio_goal + 2*ratio_lower_limit)/3.0;
      if (verbose>0) std::cout << "Decreasing ratio to " << ratio_goal << "\n";
    } else {
      return false;
    }
    
  }
  
  return false;

}




// Topological interleaving test on a bag of z(=u) and w(=v) balls.
//
// Reuses schottky's TrapGrid pipeline (fill_pixels / connected components /
// intersection boundaries) but decides via find_interleaved_components with
// ALL cut-components marked good and far_trap_points=false.  This is the pure
// arc-alternation criterion: the boundary of uL union vL has a z-component and
// a w-component that interleave (4 points in cyclic order z,w,z,w) -- CKW
// Def 7.1.3 linking.  It does NOT require fitting a robust ball into each
// protrusion (that is what the stricter find_trap_given_balls / prune_boundary
// do, and what makes them marginal on the circle |s|=1/sqrt2).  A raster arc
// alternation of >= 4 is, per the standing convention here, sufficient for the
// linking certificate; robustness/soundness of the resulting parameter ball is
// handled separately (as for the same-length covering).
//
// Renormalized-frame note (see find_trap_mixed): the grid box adapts to the
// actual copies, which sit at their true relative sizes |s|^{|u|-|v|}, so the
// interleaving is detected in the O(1) renormalized frame.
//
// Optionally returns the number of z_cut_by_w and w_cut_by_z components (a
// coarse measure of how much the copies cut each other).
bool ifs::trap_interleaves_topological(const std::vector<Ball>& balls,
                                       int max_pixels,
                                       int* n_zcut, int* n_wcut) {
  if (n_zcut) *n_zcut = 0;
  if (n_wcut) *n_wcut = 0;
  if (balls.size() == 0) return false;
  cpx ll, ur;
  box_containing_balls(balls, ll, ur);
  double av_radius = 0;
  for (int i=0; i<(int)balls.size(); ++i) av_radius += balls[i].radius;
  av_radius /= (double)balls.size();
  double desired_pixel_diameter = av_radius/2.5;
  int np = int(((ur.real() - ll.real())/desired_pixel_diameter) + 1);
  if (np > max_pixels) np = max_pixels;
  if (np < 4) np = 4;
  TrapGrid TG;
  TG.reset_grid(ll, ur, np);
  TG.fill_pixels(balls);
  TG.compute_connected_components();
  TG.compute_intersection_boundaries();
  TG.compute_distances();
  if (n_zcut) *n_zcut = (int)TG.z_cut_by_w_components.size();
  if (n_wcut) *n_wcut = (int)TG.w_cut_by_z_components.size();
  //mark every cut-component "good" (nonnegative) so the interleaving test is
  //purely topological (no robust-ball requirement)
  std::vector<int> good_z(TG.z_cut_by_w_components.size(), 1);
  std::vector<int> good_w(TG.w_cut_by_z_components.size(), 1);
  if (good_z.empty() || good_w.empty()) return false;  //copies don't cut each other
  std::vector<std::vector<Point3d<int> > > ic;
  return TG.find_interleaved_components(ic, good_z, good_w, false);
}


// ---------------------------------------------------------------------------
// CKW round-disk POINT certificate (Def 7.1.3).  Given explicit words u,v at
// the current parameter, build the round enclosing disk D=disk(1/2,r_D) with
// f D, g D subset D (r_D = |s-1|/(2(1-|s|))); then uD=disk(u(1/2),|s|^|u| r_D)
// and vD=disk(v(1/2),|s|^|v| r_D) are round disks containing uLambda, vLambda.
// Detect the topological interleaving, extract the 4 deepest poke-out points
// (one per interleaved cut-component), and test the CKW disjointness
// p^± in uLambda \ vD, q^± in vLambda \ uD by round-disk membership.  We use
// the RIGOROUS margin  (dist to the opposite disk center) - r_opposite -
// ball_radius : if positive, a genuine limit-set point inside that ball is
// provably outside the opposite disk.  Returns true iff all four rigorous
// margins are positive (=> a certificate a rigorous verifier can check).
bool ifs::ckw_point_certificate(const std::string& us, const std::string& vs,
                                int n_depth, int max_pixels, double dmult,
                                double* min_margin_out, int verbose) {
  if (min_margin_out) *min_margin_out = -1.0;
  double mir; if (!minimal_enclosing_radius(mir)) return false;
  Ball init(0.5, (z-1.0)/2.0, (1.0-w)/2.0, mir*dmult);
  Ball zb = act_on_left(us[0]-'0', init);
  for (size_t i=1;i<us.size();++i) zb = act_on_right(us[i]-'0', zb);
  Ball wb = act_on_left(vs[0]-'0', init);
  for (size_t i=1;i<vs.size();++i) wb = act_on_right(vs[i]-'0', wb);
  int lz = zb.word_len, lw = wb.word_len;
  int Lfinal = std::max(lz,lw) + n_depth;
  std::vector<Ball> ZB, WB, balls;
  compute_balls_right(ZB, zb, Lfinal-lz);
  compute_balls_right(WB, wb, Lfinal-lw);
  balls.swap(ZB); balls.insert(balls.end(), WB.begin(), WB.end());

  cpx ll, ur; box_containing_balls(balls, ll, ur);
  double avr=0; for (int i=0;i<(int)balls.size();++i) avr += balls[i].radius;
  avr /= (double)balls.size();
  int np = int(((ur.real()-ll.real())/(avr/2.5)) + 1);
  if (np > max_pixels) np = max_pixels;  if (np < 4) np = 4;
  TrapGrid TG; TG.reset_grid(ll, ur, np); TG.fill_pixels(balls);
  TG.compute_connected_components(); TG.compute_intersection_boundaries(); TG.compute_distances();
  std::vector<int> gz(TG.z_cut_by_w_components.size(),1), gw(TG.w_cut_by_z_components.size(),1);
  if (gz.empty() || gw.empty()) return false;
  std::vector<std::vector<Point3d<int> > > ic;
  if (!TG.find_interleaved_components(ic, gz, gw, false)) return false;

  //round disk D and its word-images (centers are the seed-ball centers = w(1/2))
  double az_ = abs(z);
  double rD  = 0.5*abs(z-1.0)/(1.0-az_);
  double ruD = pow(az_, lz)*rD, rvD = pow(az_, lw)*rD;
  cpx cuD = zb.center, cvD = wb.center;

  double minm = 1e18; bool allpos = true;
  for (int k=0;k<4;++k) {
    int zw = ic[0][k].x;      //0 = u-copy point, 1 = v-copy point
    int cc = ic[0][k].z;      //cut-component index
    Point2d<int> px = TG.farthest_from_other_component(zw, cc);
    int bi = (zw==0 ? TG.grid[px.x][px.y].closest_z_ball
                    : TG.grid[px.x][px.y].closest_w_ball);
    if (bi < 0) { allpos=false; continue; }
    cpx p = balls[bi].center; double br = balls[bi].radius;
    double margin = (zw==0 ? abs(p-cvD)-rvD : abs(p-cuD)-ruD) - br; //rigorous
    if (margin < minm) minm = margin;
    if (margin <= 0) allpos = false;
    if (verbose>0) std::cout << "  pt" << k << " copy=" << (zw==0?"u":"v")
                             << " margin(rig)=" << margin << " (ballr=" << br << ")\n";
  }
  if (min_margin_out) *min_margin_out = minm;
  return allpos;
}


// ---------------------------------------------------------------------------
// Mixed-length geometric trap finder.
//
// Soundness is identical to the same-length find_trap: it reuses
// find_trap_given_balls (TrapGrid + find_interleaved_components), which
// certifies a trap by exhibiting 4 interleaved trap balls on the boundary of
// uL union vL (CKW Def 7.1.3, which uses min(|u|,|v|) and is valid for
// |u| != |v|).  The ONLY difference from find_trap is the seed search: the
// two sub-copies are advanced INDEPENDENTLY (|len(u)-len(v)| <= kmax) instead
// of in lockstep, so mixed-length near-coincident pairs are found.
//
// Renormalized-frame note: we grid uL and vL at their TRUE positions and
// sizes, so their diameters are automatically in the ratio |s|^{|u|-|v|}.
// find_trap_given_balls sizes the grid via box_containing_balls and sets the
// pixel diameter to av_radius/2.5, so the grid zooms to fit both copies --
// i.e. interleaving is detected in the O(1) renormalized frame with the copies
// scaled s^{|u|-|v|} relative to one another.  (This is why we do NOT hit the
// earlier fixed-box failure where tiny copies vanished sub-pixel.)  We extend
// both copies to a common final word length so every leaf ball has the same
// radius (uniform grid resolution for both copies).
//
// u has outermost letter 0 (=f), v outermost letter 1 (=g); fill_pixels uses
// last_gen_index() (the outermost letter) to classify z(=u) vs w(=v) balls,
// so mixed lengths are transparent.  Returns true and fills u_out,v_out
// (outermost-first strings) and min_trap_dist if a trap is found.
bool ifs::find_trap_mixed(int max_uv_depth, int n_depth, int max_pixels,
                          int kmax, double ratio_goal, int min_uv,
                          double dmult,
                          std::string* u_out, std::string* v_out,
                          double* min_trap_dist_out, int verbose,
                          std::vector<std::pair<std::string,std::string> >* cand_out) {
  double min_initial_radius;
  if (!minimal_enclosing_radius(min_initial_radius)) return false;

  //connectedness at the working depth (mirrors find_trap)
  int old_depth = depth;
  depth = max_uv_depth + n_depth + 2;
  bool conn = circ_connected(min_initial_radius);
  depth = old_depth;
  if (!conn) { if (verbose>0) std::cout << "not even connected\n"; return false; }

  if (abs(z) > 1.0/sqrt(2) + 1e-12) {
    if (u_out) *u_out = "-";  if (v_out) *v_out = "-";
    return true;  //outside disk: trivially interior (not our case on the circle)
  }

  Ball initial_ball(0.5, (z-1.0)/2.0, (1.0-w)/2.0, min_initial_radius*dmult);
  Ball z0 = act_on_left(0, initial_ball);   //outermost letter 0 => u
  Ball w0 = act_on_left(1, initial_ball);   //outermost letter 1 => v

  //---- independent-advance DFS: collect near-coincident (zb,wb) candidates ----
  typedef std::pair<Ball,Ball> BP;
  std::vector<BP> cand;
  std::deque<BP> st;
  st.push_back(std::make_pair(z0, w0));
  std::set<std::pair<std::pair<unsigned long long,int>,
                     std::pair<unsigned long long,int> > > seen;
  long budget = 8000000;
  const int cand_cap = 800;
  long dbg_nodes=0; int dbg_maxlen=0; double dbg_minratio=1e18; long dbg_pushed=0;
  while (!st.empty() && budget-- > 0 && (int)cand.size() < cand_cap) {
    BP b = st.back(); st.pop_back();
    Ball& zb = b.first; Ball& wb = b.second;
    if (zb.is_disjoint(wb)) continue;
    ++dbg_nodes; if (zb.word_len>dbg_maxlen) dbg_maxlen=zb.word_len;

    std::pair<std::pair<unsigned long long,int>,std::pair<unsigned long long,int> > key(
        std::make_pair(zb.word.to_ullong(), zb.word_len),
        std::make_pair(wb.word.to_ullong(), wb.word_len));
    if (!seen.insert(key).second) continue;

    int lz = zb.word_len, lw = wb.word_len;
    //alignment ratio measured against the TRUE copy radius (mir*az^len), i.e.
    //ball radius / dmult, so ratio_goal is independent of the ball fattening
    double rmin = std::min(zb.radius, wb.radius)/dmult;
    double ratio = abs(zb.center - wb.center)/rmin;
    if (lz >= min_uv && lw >= min_uv && abs(lz-lw) <= kmax && ratio < dbg_minratio) dbg_minratio=ratio;
    if (lz >= min_uv && lw >= min_uv && abs(lz-lw) <= kmax && ratio < ratio_goal) {
      cand.push_back(b);
    }

    //advance BOTH in lockstep (keeps |lendiff| the same; this is the ONLY move
    //for kmax=0, reproducing the same-length lockstep search)
    if (lz < max_uv_depth && lw < max_uv_depth) {
      Ball za[2] = { act_on_right(0, zb), act_on_right(1, zb) };
      Ball wa[2] = { act_on_right(0, wb), act_on_right(1, wb) };
      for (int a=0;a<2;++a) for (int c=0;c<2;++c)
        if (!za[a].is_disjoint(wa[c])) st.push_back(std::make_pair(za[a], wa[c]));
    }
    //for kmax>=1, also advance u (z) alone (diff += 1) and v (w) alone (diff -= 1)
    if (kmax >= 1 && lz < max_uv_depth && (lz+1) - lw <= kmax) {
      Ball c0 = act_on_right(0, zb), c1 = act_on_right(1, zb);
      if (!c0.is_disjoint(wb)) st.push_back(std::make_pair(c0, wb));
      if (!c1.is_disjoint(wb)) st.push_back(std::make_pair(c1, wb));
    }
    if (kmax >= 1 && lw < max_uv_depth && (lw+1) - lz <= kmax) {
      Ball c0 = act_on_right(0, wb), c1 = act_on_right(1, wb);
      if (!c0.is_disjoint(zb)) st.push_back(std::make_pair(zb, c0));
      if (!c1.is_disjoint(zb)) st.push_back(std::make_pair(zb, c1));
    }
  }

  if (verbose>0) std::cout << "collected " << cand.size() << " candidate mixed pairs"
                           << "  [nodes=" << dbg_nodes << " maxlen=" << dbg_maxlen
                           << " pushed=" << dbg_pushed << " minratio(filtered)="
                           << dbg_minratio << " budget_left=" << budget << "]\n";
  if (cand.empty()) return false;

  //sort candidates by alignment ratio (closest-aligned first) -- no lambda for
  //compiler-portability
  std::vector<std::pair<double,int> > order(cand.size());
  for (int i=0; i<(int)cand.size(); ++i) {
    double r = abs(cand[i].first.center - cand[i].second.center) /
               (std::min(cand[i].first.radius, cand[i].second.radius)/dmult);
    order[i] = std::make_pair(r, i);
  }
  std::sort(order.begin(), order.end());

  //---- test candidates with the geometric interleaving check ----
  const int test_cap = 150;
  int ntest = std::min((int)order.size(), test_cap);
  // If cand_out is given, collect up to cand_max passing candidates (in
  // alignment order) for the caller to rank by certificate margin; otherwise
  // return the first passing pair (legacy behavior).
  const int cand_max = cand_out ? 8 : 1;
  int collected = 0;
  for (int oi=0; oi<ntest; ++oi) {
    const BP& c = cand[order[oi].second];
    Ball zb = c.first, wb = c.second;
    int lz = zb.word_len, lw = wb.word_len;
    int Lfinal = std::max(lz, lw) + n_depth;
    std::vector<Ball> ZB, WB, balls;
    compute_balls_right(ZB, zb, Lfinal - lz);
    compute_balls_right(WB, wb, Lfinal - lw);
    balls.swap(ZB);
    balls.insert(balls.end(), WB.begin(), WB.end());
    int nzc=0, nwc=0;
    bool got = trap_interleaves_topological(balls, max_pixels, &nzc, &nwc);
    if (got) {
      std::string us = Bitword(zb.word, lz).str();
      std::string vs = Bitword(wb.word, lw).str();
      if (collected==0) {
        if (u_out) *u_out = us;
        if (v_out) *v_out = vs;
        if (min_trap_dist_out) *min_trap_dist_out = (double)std::min(nzc,nwc);
        if (verbose>0)
          std::cout << "TRAP u=" << us << " v=" << vs << " (|u|=" << lz
                    << " |v|=" << lw << ") ratio=" << order[oi].first
                    << " zcut=" << nzc << " wcut=" << nwc << "\n";
      }
      if (cand_out) cand_out->push_back(std::make_pair(us,vs));
      ++collected;
      if (collected >= cand_max) return true;
    }
  }
  return collected>0;
}




bool ifs::find_traps_along_loop(const std::vector<cpx>& loop,
                                bool draw_it, 
                                int verbose) {
  
  //trap parameters
  int max_pixels = 1024;
  int max_uv_depth = 25;
  int max_n_depth = 20;
  
  
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
  
  //int tv = (verbose > 0 ? verbose -1 : 0);
  double Cz = 0;
  for (int i=0; i<(int)loop.size(); ++i) {
    if (1.0/(1.0-abs(loop[i])) > Cz) {
      Cz = 1.0/(1.0-abs(loop[i]));
    }
  }
  if (verbose>0) {
    std::cout << "Found Cz = " << Cz << "\n";
  }
  
  //get graphics stuff
  //int rcol = X.get_rgb_color(1,0,0);
  double pixel_width = (2*wind)/double(drawing_width);
  double difficulty;
  
  //get the traps at the vertices
  for (int i=0; i<nL; ++i) {
    trap_list[i].resize(1);
    z = loop[i]; az = abs(z);
    w = z; aw = az;
    double epsilon;
    if (!find_trap(max_uv_depth, max_n_depth, max_pixels, Cz, &epsilon, &difficulty, verbose)) {
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
      int c = X.get_rgb_color(1.0, gamount, 0.0);
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
        if (verbose>0) std::cout << "Done this edge\n";
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
      if (!find_trap(max_uv_depth, max_n_depth, max_pixels, Cz, &trap_list[i].back().second, &difficulty, verbose)) {
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
        int c = X.get_rgb_color(1.0, gamount, 0.0);
        X.draw_disk(p,r,c);
        X.flush();
      }
#endif
      if (verbose>0) {
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


/* Empty, and deliberately left so.  This belongs to ifs::user_interface() -- the
   pre-2014 text-driven interface in ifs_draw.cc and ifs_interface.cc -- which nothing
   calls any more: the only reference to it is a commented-out line in schottky.cc.  No
   version of the program has ever drawn a trap through this function, and filling it in
   would add code that cannot be reached or tested.

   The live implementation is IFSGui::draw_limit_trap() in ifs_gui.cc, reached from the
   "Trap" checkbox on the limit-set pane.  If the text interface is ever revived, copy
   that, replacing the Widget/Pixmap calls with this class's own X member. */
void ifs::draw_trap() {

}


























