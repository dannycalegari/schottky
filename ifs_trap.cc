//Functions to check for a trap


#include "trap_grid.h"

struct Trap {
};





bool ifs::find_trap_given_balls(const std::vector<Ball>& initial_balls, 
                                int max_refinements,
                                int max_pixels,
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
      TG.show(NULL, NULL, NULL, NULL);
      if (verbose>1) TG.show(NULL, &balls, NULL, NULL);
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
    std::vector<bool> good_components[2];
    std::vector<Ball> good_balls[2];
    good_components[0].resize(TG.z_cut_by_w_components.size());
    good_balls[0].resize(TG.z_cut_by_w_components.size());
    good_components[1].resize(TG.w_cut_by_z_components.size());
    good_balls[1].resize(TG.w_cut_by_z_components.size());
    for (int i=0; i<(int)TG.z_cut_by_w_components.size(); ++i) {
      Point2d<int> p = TG.farthest_from_other_component(0, i);
      if (TG.grid[p.x][p.y].w_distance == 0) {
        good_components[0][i] = false;
        continue;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_z_ball];
      Ball b2 = b;
      for (int k=0; k<3; ++k) {
        b = act_on_right(0,b);
        b2 = act_on_right(1,b2);
      }
      if (TG.disjoint_from_z_or_w(b, 1)) {
        good_balls[0][i] = b;
        good_components[0][i] = true;
      } else if (TG.disjoint_from_z_or_w(b2, 1)) {
        good_balls[0][i] = b2;
        good_components[0][i] = true;
      } else {
        good_components[0][i] = false;
      }
    }
    for (int i=0; i<(int)TG.w_cut_by_z_components.size(); ++i) {
      Point2d<int> p = TG.farthest_from_other_component(1, i);
      if (TG.grid[p.x][p.y].z_distance == 0) {
        good_components[1][i] = false;
        continue;
      }
      Ball b = balls[TG.grid[p.x][p.y].closest_w_ball];
      Ball b2 = b;
      for (int k=0; k<3; ++k) {
        b = act_on_right(0,b);
        b2 = act_on_right(1,b2);
      }
      if (TG.disjoint_from_z_or_w(b, 0)) {
        good_balls[1][i] = b;
        good_components[1][i] = true;
      } else if (TG.disjoint_from_z_or_w(b2, 0)) {
        good_balls[1][i] = b2;
        good_components[1][i] = true;
      } else {
        good_components[1][i] = false;
      }
    } 
    
    if (verbose>0) {
      std::cout << "Found good components:\n";
      for (int i=0; i<(int)good_components[0].size(); ++i) {
        if (good_components[0][i]) {
          std::cout << "(0," << i << ") ";
        }
      }
      std::cout << "\n";
      for (int i=0; i<(int)good_components[1].size(); ++i) {
        if (good_components[1][i]) {
          std::cout << "(1," << i << ") ";
        }
      }
      std::cout << "\n";
    }
    //now let's find interleaved components
    //if we find interleaved good components, we're done
    std::vector<std::vector<Point3d<int> > > ic;
    if (!TG.find_interleaved_components(ic, good_components[0], good_components[1])) {
      if (verbose>0) std::cout << "Didn't find interleaved components\n";
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
      if (verbose>0) {
        trap_balls[0] = good_balls[ic[0][0].x][ic[0][0].z];
        trap_balls[1] = good_balls[ic[0][1].x][ic[0][1].z];
        trap_balls[2] = good_balls[ic[0][2].x][ic[0][2].z];
        trap_balls[3] = good_balls[ic[0][3].x][ic[0][3].z];
        std::cout << "Found four disjoint balls!\n";
        TG.show(NULL, &trap_balls, NULL, NULL);
      }
      //std::cout << "yes!\n";
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
      TG.show(NULL,NULL,&new_ll, &new_ur);
    }
    
    //might as well refine the depth -- compute the balls which fit into
    //this box
    refine_balls_into_box(balls, new_ll, new_ur);
    
    //when we loop around, these are the new balls
    
  }//<- end of trap searching
  
  return false;
  
}








bool ifs::find_trap(int max_uv_depth, int max_n_depth, double* epsilon, int verbose) {

  int uv_depth = 2*depth;
  int n_depth = depth;

  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double min_initial_radius;
  if (!minimal_enclosing_radius(min_initial_radius)) {
    //std::cout << "initial radius is infinite\n";
    return false;
  }
  //make sure it's connected at the minimal radius
  if (verbose>0) std::cout << "Checking connectedness with minimal initial radius of " << min_initial_radius << "\n";
  if (!circ_connected(min_initial_radius)) {
    if (verbose>0) {
      std::cout << "Not even connected\n";
    }
    return false;
  }
  
  while (true) {
  
    //find actions u and v which start with z and w such that 
    //u(1/2) and v(1/2) are well-aligned
    Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_initial_radius*2);
    Ball zb, wb;
    
    //find_aligned_images_with_distinct_first_letters(initial_ball, z_cm, w_cm, uv_depth, zb, wb);
    if (abs(z) > (1.0/sqrt(2))) {
      if (verbose>0) {
        std::cout << "|z| is larger than 1/sqrt(2), so we're done\n";
      }
      return true;
    }
    find_aligned_images_with_distinct_first_letters(initial_ball, 0, 0, uv_depth, zb, wb, 0.05);
    //find_close_images_with_distinct_first_letters(initial_ball, uv_depth, zb, wb);
    
    
    if (verbose>0) {
      std::cout << "Found initial balls\n";
      std::cout << "z-ball: " << zb << "\nw-ball: " << wb << "\n";
    }
    
    //now act on the *right* a bunch
    std::vector<Ball> ZB;
    compute_balls_right(ZB, zb, n_depth);
    std::vector<Ball> WB;
    compute_balls_right(WB, wb, n_depth);
    
    //now concatenate the lists and use that (we swap to avoid unnecessary work)
    std::vector<Ball> balls(0);
    balls.swap(ZB);
    balls.insert(balls.end(), WB.begin(), WB.end());
    
    //std::cout << "Computed " << balls.size() << " balls " << balls[0] << "\n";
    
    //trap finding parameters:
    int max_refinements = 0;
    int max_pixels = 800;
    
    bool got_trap = find_trap_given_balls(balls, max_refinements, max_pixels, verbose);
    
    if (got_trap) {
      if (epsilon != NULL) *epsilon = pow(az, uv_depth + n_depth);
      return true;
    }
    
  }
  
  return false;

}




bool ifs::find_traps_along_loop(const std::vector<cpx>& loop, 
                                bool draw_it, 
                                int verbose) {
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
  
  int tv = (verbose > 0 ? verbose -1 : 0);
  
  //get the traps at the vertices
  for (int i=0; i<nL; ++i) {
    trap_list[i].resize(1);
    z = loop[i]; az = abs(z);
    w = z; aw = az;
    double epsilon;
    if (!find_trap(&epsilon, tv)) {
      if (verbose>0) std::cout << "Failed to find a trap at vertex " << i << "\n";
      return false;
    }
    trap_list[i][0] = std::make_pair(z, epsilon);
    if (verbose>0) std::cout << i << ": " << trap_list[i][0].first << ", " << trap_list[i][0].second << "\n";
  }
  
  //for each interval, go along it, placing the center of the 
  //next trap at exactly the edge of the previous one
  int rcol = X.get_rgb_color(1,0,0);
  double pixel_radius = double(drawing_radius)/wind;
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
      if (!find_trap(&trap_list[i].back().second, tv)) {
        if (verbose>0) std::cout << "Failed to find trap at " << z << "\n";
        return false;
      }
      //display it
      if (draw_it) {
        Point2d<int> p = cpx_to_point_mandlebrot(z);
        double r = trap_list[i].back().second / pixel_radius;
        X.draw_disk(p,r,rcol);
      }
      if (verbose>0) {
        std::cout << "Found new trap " << trap_list[i].back().first << ", " << trap_list[i].back().second << "\n";
      }
    }
  }
  if (verbose>0) std::cout << "Loop certified!\n";
  return true;
    
  
}


void ifs::draw_trap() {

}


























