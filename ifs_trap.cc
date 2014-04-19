//Functions to check for a trap


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






bool ifs::find_trap(int max_uv_depth, int max_n_depth, int max_pixels, double Cz, double* epsilon, int verbose) {

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
  depth = max_uv_depth+max_n_depth+2;
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
  double p = 2.0;
  int current_n_depth = 10; //I guess always start here?
  int current_uv_depth = -1;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_initial_radius*p);
  Ball zb, wb;
  
  while (true) {
  
    //find actions u and v which start with z and w such that 
    //find_aligned_images_with_distinct_first_letters(initial_ball, z_cm, w_cm, uv_depth, zb, wb);
    if (abs(z) > (1.0/sqrt(2))) {
      if (verbose>0) {
        std::cout << "|z| is larger than 1/sqrt(2), so we're done\n";
      }
      if (epsilon != NULL) *epsilon = 0.005;
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
    double min_trap_dist;
    
    bool got_trap = find_trap_given_balls(balls, max_pixels, &min_trap_dist, (verbose>0 ? verbose-1 : verbose));
    if (verbose>0) {
      if (got_trap) {
        std::cout << "Got trap at " << z << " with min trap dist " << min_trap_dist << "\n";
      } else {
        std::cout << "Failed to find trap\n";
      }
    }
    
    if (got_trap) {
      double ball_radius = balls[0].radius;
      double dr = min_trap_dist; //this is the largest radius we could have at this point, which we define as dr
      double r = pow(az, current_n_depth + current_uv_depth)*min_initial_radius;
      double d = dr/r;
      double e = (1.0/Cz)*0.25*r*(d-p);
      if (verbose>0) {
        std::cout << "ball radius: " << ball_radius << " r = " << r << " d = " << d << " p= " << p << " epsilon = " << e << "\n";
      }
      if (epsilon != NULL) *epsilon = e;
      return true;
    }
    //if we didn't get it, think about what to do
    if (current_n_depth < max_n_depth) {
      current_n_depth += (current_n_depth < max_n_depth-1 ? 3 : 1);
      if (verbose>0) std::cout << "Increasing depth to " << current_n_depth << "\n";
    } else {
      return false;
    }
    
  }
  
  return false;

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
  
  
  //get the traps at the vertices
  for (int i=0; i<nL; ++i) {
    trap_list[i].resize(1);
    z = loop[i]; az = abs(z);
    w = z; aw = az;
    double epsilon;
    if (!find_trap(max_uv_depth, max_n_depth, max_pixels, 3.42, &epsilon, verbose)) {
      if (verbose>0) std::cout << "Failed to find a trap at vertex " << i << "\n";
      return false;
    }
    trap_list[i][0] = std::make_pair(z, epsilon);
    if (verbose>0) std::cout << i << ": " << trap_list[i][0].first << ", " << trap_list[i][0].second << "\n";
  }
  
  //for each interval, go along it, placing the center of the 
  //next trap at exactly the edge of the previous one
  int rcol = X.get_rgb_color(1,0,0);
  double pixel_width = (2*wind)/double(drawing_width);
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
      if (!find_trap(max_uv_depth, max_n_depth, max_pixels, 3.42, &trap_list[i].back().second, verbose)) {
        if (verbose>0) std::cout << "Failed to find trap at " << z << "\n";
        return false;
      }
      //display it
      if (draw_it) {
        Point2d<int> p = cpx_to_point_mandlebrot(z);
        double r = trap_list[i].back().second / pixel_width;
        if (r<2) r = 2;
        X.draw_disk(p,r,rcol);
        X.flush();
      }
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


void ifs::draw_trap() {

}


























