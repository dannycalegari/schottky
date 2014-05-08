


bool ifs::hole_boundary_containing_point(std::vector<cpx>& path, bool& closed, 
                                         cpx p, int verbose) {
  if (verbose>0) {
    std::cout << "Rerunning mandlebrot to get the pixel grid...\n";
  }
  int num_pixels = drawing_width/mesh;
  std::vector<std::vector<bool> > grid(num_pixels, std::vector<bool>(num_pixels, false));
  cpx ll = center - cpx(wind, wind);
  cpx ur = center + cpx(wind, wind);
  double wid = 2*wind;
  double pixel_width = wid/double(num_pixels);
  
  //make a new ifs so we don't screw stuff up
  ifs IFS2(z,w,drawing_width,mode);
  IFS2.depth = depth;
  
  //find the center pixel
  double min_r;
  Point2d<int> marked_pixel(int( (p.real() - ll.real()) / pixel_width ),
                            int( (p.imag() - ll.imag()) / pixel_width ) );
  IFS2.set_params(p,p);
  if (verbose>0) {
    std::cout << "Marked pixel: " << marked_pixel << " p = " << p << "\n";
    std::cout << "Center: " << center << "\n";
    std::cout << "wind: " << wind << "\n";
  }
  if (!IFS2.minimal_enclosing_radius(min_r) || IFS2.circ_connected(min_r)) {
    if (verbose>0) std::cout << "Connected or weird\n";
    return false;
  }
  
  //recompute everything
  for (int i=0; i<num_pixels; ++i) {
    for (int j=0; j<num_pixels; ++j) {
      cpx Z(ll.real() + i*pixel_width, ll.imag() + j*pixel_width);
      IFS2.set_params(Z,Z);
      if (!IFS2.minimal_enclosing_radius(min_r)) {
        grid[i][j] = true;
      }
      grid[i][j] = IFS2.circ_connected(min_r);
    }
  }
  
  return hole_boundary_containing_point_from_grid(path, closed, grid, marked_pixel, ll, ur, verbose);
}




//return a path which is the boundary of the hole containing 
//the point p (false if p isn't *in* a hole)
//if the hole touches the boundary of the current window, 
//the path is not a loop, otherwise, it is
bool ifs::hole_boundary_containing_point_from_grid(std::vector<cpx>& path, bool& closed, 
                                                   const std::vector<std::vector<bool> >& grid,
                                                   const Point2d<int>& marked_pixel, 
                                                   cpx ll, cpx ur, int verbose) {
  int num_pixels = grid.size();
  double w = (ur.real() - ll.real());
  double pixel_width = w/double(grid.size());
  //follow the boundary
  //the direction records which way we're looking.  0 is down, going ccw around
  std::vector<Point2d<int> > pixel_boundary(0);
  int offset_i[4] = {0,1,0,-1};
  int offset_j[4] = {-1,0,1,0};
  int i=marked_pixel.x;
  int j=marked_pixel.y;
  
  if (grid[i][j]) {
    std::cout << "Marked pixel is not in a hole?\n";
    path.resize(0);
    closed = false;
    return false;
  }
  
  
  while (0 <= j && j < num_pixels && grid[i][j] == false) {
    j += 1;
  }
  --j; //step back so we're on a good pixel
  int start_i = i;
  int start_j = j;
  int start_dir = 0; 
  //make the start direction correct -- the direction we would have come in from
  int ii = i + offset_i[start_dir];
  int jj = j + offset_j[start_dir];
  while ( 0 <= ii && ii < num_pixels && 0 <= jj && jj < num_pixels && grid[ii][jj] == false ) {
    start_dir = (start_dir == 0 ? 3 : start_dir-1);
    ii = i + offset_i[start_dir];
    jj = j + offset_j[start_dir];
  }
  start_dir =  (start_dir == 3 ? 0 : start_dir+1);
  int current_dir = start_dir;
  
  if (verbose>0) {
    std::cout << "Marked pixel: " << marked_pixel << "\n";
    std::cout << "Found the start pixel, dir " << i << " " << j << " " << start_dir << "\n";
  }
  
  do {
    //push the current pixel
    pixel_boundary.push_back(Point2d<int>(i,j));
    
    //scan to find the next pixel
    do{ 
      current_dir = (current_dir == 0 ? 3 : current_dir-1);
      ii = i + offset_i[current_dir];
      jj = j + offset_j[current_dir];
    } while (ii < 0 || ii >= num_pixels || 0 > jj || jj >= num_pixels || grid[ii][jj] == true);
    
    //make that the current pixel
    i = ii;
    j = jj;
    current_dir = (current_dir + 2)%4;
  
  } while (i != start_i || j != start_j || current_dir != start_dir);
  
  
  //get the real boundary along set A
  int real_boundary_start = -1;
  for (int i=(int)pixel_boundary.size()-1; i>=0; --i) {
    int im1 = (i==0 ? pixel_boundary.size()-1 : i-1);
    if ( (pixel_boundary[im1].x == 0 || pixel_boundary[im1].x == num_pixels-1 ||
          pixel_boundary[im1].y == 0 || pixel_boundary[im1].y == num_pixels-1))  {
      real_boundary_start = i;
      break;
    }
  }
  if (real_boundary_start == -1) { //there's no window boundary
    path.resize(pixel_boundary.size());
    for (int i=0; i<(int)pixel_boundary.size(); ++i) {
      path[i] = cpx(ll.real() + (pixel_boundary[i].x+0.5)*pixel_width, 
                    ll.imag() + (pixel_boundary[i].y+0.5)*pixel_width);
    }
    closed = true;
  } else {
    path.resize(0);
    int i=real_boundary_start;
    while (pixel_boundary[i].x > 0 && pixel_boundary[i].x < num_pixels-1 &&
           pixel_boundary[i].y > 0 && pixel_boundary[i].y < num_pixels-1) {
      path.push_back( cpx(ll.real() + (pixel_boundary[i].x+0.5)*pixel_width, 
                          ll.imag() + (pixel_boundary[i].y+0.5)*pixel_width) );
      i = (i == (int)pixel_boundary.size()-1 ? 0 : i+1);
    }
    closed = false;
  }
  
  //display the boundary
  if (verbose>0) {
    std::cout << "Boundary has length " << pixel_boundary.size() << "\n";
  
    int pixel_group_width = (num_pixels > 512 ? 1 : 512/num_pixels);
    int X2_num_pixels = (num_pixels > 512 ? num_pixels : (512/num_pixels)*num_pixels);
    XGraphics X2(X2_num_pixels, X2_num_pixels, 1, Point2d<float>(0,0));
    int wcol = X2.get_rgb_color(1,1,1);
    int bcol = X2.get_rgb_color(0,0,0);
    int rcol = X2.get_rgb_color(1,0.2,0);
    int gcol = X2.get_rgb_color(0,1,0.2);
    for (int i=0; i<num_pixels; ++i) {
      for (int j=0; j<num_pixels; ++j) {
        if (num_pixels > 512) {
          X2.draw_point(Point2d<int>(i,j), (grid[i][j] ? bcol : wcol));
        } else {
          X2.draw_box(Point2d<int>(i*pixel_group_width,j*pixel_group_width), pixel_group_width, (grid[i][j] ? bcol : wcol));
        }
      }
    }
    for (int i=0; i<(int)pixel_boundary.size(); ++i) {
      int x = pixel_boundary[i].x;
      int y = pixel_boundary[i].y;
      if (num_pixels > 512) {
        X2.draw_point(Point2d<int>(x,y), rcol);
      } else {
        X2.draw_box(Point2d<int>(x*pixel_group_width,y*pixel_group_width), pixel_group_width, rcol);
      }
    }
    for (int i=0; i<(int)path.size(); ++i) {
      int x = int( (path[i].real()-ll.real()) / pixel_width);
      int y = int( (path[i].imag()-ll.imag()) / pixel_width);
      if (num_pixels > 512) {
        X2.draw_point(Point2d<int>(x,y), gcol);
      } else {
        X2.draw_box(Point2d<int>(x*pixel_group_width,y*pixel_group_width), pixel_group_width, gcol);
      }
    }
    if (num_pixels > 512) {
      X2.draw_point(marked_pixel, rcol);
    } else {
      X2.draw_box(Point2d<int>(marked_pixel.x*pixel_group_width,marked_pixel.y*pixel_group_width), pixel_group_width, rcol);
    }
    X2.wait_for_key();
  }
    
  return true;
}





















