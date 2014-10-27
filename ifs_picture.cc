

void ifs::draw_ifs_to_array(std::vector<std::vector<Point3d<unsigned char> > >& bmp, 
                            const cpx& region_ll, const cpx& region_ur, 
                            int depth) {
  double min_r;
  if (!minimal_enclosing_radius(min_r)) return;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  
  //clear the array
  for (int i=0; i<(int)bmp.size(); ++i) {
    for (int j=0; j<(int)bmp[i].size(); ++j) {
      bmp[i][j] = Point3d<unsigned char>(255,255,255);
    }
  }
  
  double drawing_width = region_ur.real() - region_ll.real();
  int num_pixels = bmp.size();
  double pixel_width = drawing_width/double(num_pixels);
  
  std::deque<Ball> stack(1);
  stack[0] = initial_ball;
  while (stack.size() > 0) {
    Ball b = stack.back();
    stack.pop_back();
    if ( (depth > 0 && b.word_len >= depth) || 
         (depth < 0 && b.radius < pixel_width/2.0) ) {
      int x = int(num_pixels*((b.center.real() - region_ll.real()) / drawing_width));
      int y = int(num_pixels*((b.center.imag() - region_ll.imag()) / drawing_width));
      if (0 <= x && x < num_pixels && 0 <= y && y < num_pixels) {
        bmp[x][y] = (b.last_gen_index() == 0 ? Point3d<unsigned char>(0xFF, 0xAA, 0x00) :
                                               Point3d<unsigned char>(0x00, 0xAA, 0xFF));
      }
    } else {
      Ball bz = act_on_right(0, b);
      Ball bw = act_on_right(1, b);
      if (!is_ball_disjoint(bz, region_ll, region_ur)) stack.push_front(bz);
      if (!is_ball_disjoint(bw, region_ll, region_ur)) stack.push_front(bw);
    }
  }
    
}


void ifs::draw_mand_to_array(std::vector<std::vector<Point3d<unsigned char> > >&bmp,
                             const cpx& region_ll, const cpx& region_ur,
                             int connected_depth, int contains_half_depth) {
  double drawing_width = region_ur.real() - region_ll.real();
  int num_pixels = bmp.size();
  double pixel_width = drawing_width/double(num_pixels);
  int num_pixels_width = bmp.size();
  
  std::vector<std::vector<Point2d<int> > > raw_data(num_pixels_width);
  for (int i=0; i<num_pixels_width; ++i) {
    raw_data[i] = std::vector<Point2d<int> >(num_pixels, Point2d<int>(-1,-1));
  }
  
  ifs temp_IFS;
  
  std::cout << "Writing picture to " << num_pixels_width << " x " << num_pixels << " array.\n";
  std::cout << "Connected depth: " << connected_depth << "\n";
  std::cout << "Contains half dept: " << contains_half_depth << "\n";
  std::cout << "Computing array\n";
  
  for (int i=0; i<num_pixels_width; ++i) {
    double x = region_ll.real() + (i+0.5)*pixel_width;
    //double theta = PI/2.0 - (i+0.5)*pixel_theta_width;
    std::cout << "\rColumn " << i << " of " << num_pixels_width;  std::cout.flush();
    
    for (int j=0; j<num_pixels; ++j) {
      double y = region_ll.imag() + (j+0.5)*pixel_width;
      //double r = radius_initial + (j+0.5)*pixel_radius_width; 
      //double x = r*cos(theta);
      //double y = r*sin(theta);
      temp_IFS.set_params( cpx(x,y), cpx(x,y) );
      
      int diff;
        
      if ((connected_depth>0) && temp_IFS.is_connected(connected_depth, diff)) {
        raw_data[i][j].x = diff;
      }
      
      if ((contains_half_depth>0) && temp_IFS.contains_half(contains_half_depth, diff)) {
        raw_data[i][j].y = diff;
        //std::cout << "half depth: " << raw_data[i][j].y << "\n";
      }
    }
  }
  std::cout << "\n";
  
  int min_conn_diff=-1;
  int max_conn_diff=-1;  
  int min_half_diff=-1;
  int max_half_diff=-1;
  for (int i=0; i<num_pixels_width; ++i) {
    for (int j=0; j<num_pixels; ++j) { 
      if (raw_data[i][j].x > -1 && (min_conn_diff == -1 || raw_data[i][j].x < min_conn_diff)) {
        min_conn_diff = raw_data[i][j].x;
      }
      if (raw_data[i][j].x > -1 && (max_conn_diff == -1 || raw_data[i][j].x > max_conn_diff)) {
        max_conn_diff = raw_data[i][j].x;
      }
      if (raw_data[i][j].y > -1 && (min_half_diff == -1 || raw_data[i][j].y < min_half_diff)) {
        min_half_diff = raw_data[i][j].y;
      }
      if (raw_data[i][j].y > -1 && (max_half_diff == -1 || raw_data[i][j].y > max_half_diff)) {
        max_half_diff = raw_data[i][j].y;
      }
    }
  }
  
  std::cout << "Computed connectedness range: " << min_conn_diff << " " << max_conn_diff << "\n";
  std::cout << "Computed half range: " << min_half_diff << " " << max_half_diff << "\n";
  
  int conn_offset = min_conn_diff;
  int conn_range = max_conn_diff - min_conn_diff;
  double conn_scaling_factor = (max_conn_diff == min_conn_diff ? 1 : (256.0 / double(max_conn_diff - min_conn_diff)));
  int half_offset = min_half_diff;
  int half_range = max_half_diff - min_half_diff;
  double half_scaling_factor = (max_half_diff == min_half_diff ? 1 : (256.0 / double(max_half_diff - min_half_diff)));
  
  std::cout << "Computed connectedness scaling: +" << conn_offset << " *" << conn_scaling_factor << "\n";
  std::cout << "Computed half scaling: +" << half_offset << " *" << half_scaling_factor << "\n";
  
  for (int i=0; i<num_pixels_width; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      bmp[i][j].x = bmp[i][j].y = bmp[i][j].z = 255;
      if (raw_data[i][j].x != -1) {
        bmp[i][j].z = int(56+200.0*pow( (double(raw_data[i][j].x - conn_offset)/double(conn_range)), 0.6));
        bmp[i][j].x = bmp[i][j].y = 0;
      } if (raw_data[i][j].y != -1) {
        bmp[i][j].x = int(56+200.0*pow( (double(raw_data[i][j].y - half_offset)/double(half_range)), 0.1));
      }
    }
  }
  
  std::cout << "Wrote to array\n"; std::cout.flush();
}


void ifs::draw_mand_to_array_radial(std::vector<std::vector<Point3d<unsigned char> > >&bmp,
                                    int num_pixels, //height
                                    double radius_start, double radius_end,
                                    double theta_start, double theta_end,
                                    int connected_depth, int contains_half_depth) {
  
  
  //radial picture
  double radius_length = radius_end-radius_start;
  double theta_length = theta_start - theta_end;
  int num_pixels_width = int((theta_length/radius_length)*num_pixels);
  double pixel_theta_width = theta_length/num_pixels_width;
  double pixel_radius_width = radius_length/num_pixels;
  bmp.resize(num_pixels_width);
  for (int i=0; i<num_pixels_width; ++i) {
    bmp[i].resize(num_pixels);
  }
  std::vector<std::vector<Point2d<int> > > raw_data(num_pixels_width);
  for (int i=0; i<num_pixels_width; ++i) {
    raw_data[i] = std::vector<Point2d<int> >(num_pixels, Point2d<int>(-1,-1));
  }
  
  ifs temp_IFS;
  
  std::cout << "Writing picture to " << num_pixels_width << " x " << num_pixels << " array.\n";
  std::cout << "Connected depth: " << connected_depth << "\n";
  std::cout << "Contains half dept: " << contains_half_depth << "\n";
  std::cout << "Computing array\n";
  
  for (int i=0; i<num_pixels_width; ++i) {
    double theta = theta_start - (i+0.5)*pixel_theta_width;
    std::cout << "\rColumn " << i << " of " << num_pixels_width;  std::cout.flush();
    
    for (int j=0; j<num_pixels; ++j) {
      double r = radius_start + (j+0.5)*pixel_radius_width; 
      double x = r*cos(theta);
      double y = r*sin(theta);
      temp_IFS.set_params( cpx(x,y), cpx(x,y) );
      
      int diff;
        
      if (connected_depth>0) { 
        if (temp_IFS.is_connected(connected_depth, diff)) {
          raw_data[i][j].x = diff;
        } else {
          raw_data[i][j].x = -diff;
        }
      }
      
      if (contains_half_depth>0 && temp_IFS.contains_half(contains_half_depth, diff)) {
        raw_data[i][j].y = diff;
      }
    }
  }
  std::cout << "\n";
  
  for (int i=0; i<num_pixels_width; ++i) {
    int lower_x = (i-50 < 0 ? 0 : i-50);
    int upper_x = (i+50 >= num_pixels_width ? num_pixels_width : i+50);
    std::cout << "\rCalibrating column " << i << " of " << num_pixels_width;
    std::cout.flush();
    for (int j=0; j<num_pixels; ++j) {
      if (raw_data[i][j].x == -1 && raw_data[i][j].y == -1) {
        bmp[i][j].x = bmp[i][j].y = bmp[i][j].z = 255;
        continue;
      }
      int lower_y = (j-20 < 0 ? 0 : j-20);
      int upper_y = (j+20 >= num_pixels ? num_pixels : j+20);
      Point2d<int> conn_extremes(-1,-1);
      Point2d<int> half_extremes(-1,-1);
      for (int k=lower_x; k<upper_x; ++k) {
        for (int m=lower_y; m<upper_y; ++m) {
          if (raw_data[k][m].x > -1) {
            if (conn_extremes.x == -1 || abs(raw_data[k][m].x) < conn_extremes.x) 
              conn_extremes.x = abs(raw_data[k][m].x);
            if (conn_extremes.y == -1 || abs(raw_data[k][m].x) > conn_extremes.y)
              conn_extremes.y = abs(raw_data[k][m].x);
          }
          if (raw_data[k][m].y > -1) {
            if (half_extremes.x == -1 || raw_data[k][m].y < half_extremes.x)
              half_extremes.x = raw_data[k][m].y;
            if (half_extremes.y == -1 || raw_data[k][m].y > half_extremes.y)
              half_extremes.y = raw_data[k][m].y;
          }
        }
      }
      int conn_offset = conn_extremes.x;
      int conn_range = conn_extremes.y - conn_extremes.x;
      double conn_scaling_factor = (conn_range==0 ? 1 : (256.0 / double(conn_range)));
      int half_offset = half_extremes.x;
      int half_range = half_extremes.y - half_extremes.x;
      double half_scaling_factor = (half_range==0 ? 1 : (256.0 / double(half_range)));
      //std::cout << "Got range : " << conn_extremes << " " << conn_offset <<  " " << conn_scaling_factor << "\n";
      bmp[i][j].x = bmp[i][j].y = bmp[i][j].z = 255;
      if (connected_depth > 0) {
        if (raw_data[i][j].x > 0) { 
          bmp[i][j].z = int(56+200.0*pow( (double(raw_data[i][j].x - conn_offset)/double(conn_range)), 0.6));
          bmp[i][j].x = bmp[i][j].y = 0;
        } else {
          //bmp[i][j].y = int(56+200.0*pow( (double(-raw_data[i][j].x - conn_offset)/double(conn_range)), 0.6));
          bmp[i][j].x = bmp[i][j].z = bmp[i][j].y = 255;
        }
      } if (raw_data[i][j].y != -1) {
        bmp[i][j].x = int(56+200.0*pow( (double(raw_data[i][j].y - half_offset)/double(half_range)), 0.1));
      }
    }
  }
  
  std::cout << "\nWrote to array\n"; std::cout.flush();
}

