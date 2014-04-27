


//return a path which is the boundary of the hole containing 
//the point p (false if p isn't *in* a hole)
//if the hole touches the boundary of the current window, 
//the path is not a loop, otherwise, it is
bool ifs::hole_boundary_containing_point(std::vector<cpx>& path, bool& closed, 
                                         cpx p, int verbose) {
  if (verbose>0) {
    std::cout << "Rerunning mandlebrot to get the pixel grid...\n";
  }
  std::vector<std::vector<bool> > grid(drawing_width, std::vector<bool>(drawing_width, false));
  cpx ll = center - cpx(wind, wind);
  double w = 2*wind;
  double pixel_width = w/double(drawing_width);
  
  //make a new ifs so we don't screw stuff up
  ifs IFS2(z,w,drawing_width,mode);
  IFS2.depth = depth;
  
  //recompute everything
  double min_r;
  for (int i=0; i<drawing_width; ++i) {
    for (int j=0; j<drawing_width; ++j) {
      cpx Z(ll.real() + i*pixel_width, ll.imag() + j*pixel_width);
      IFS2.set_params(Z,Z);
      if (!IFS2.minimal_enclosing_radius(min_r)) {
        grid[i][j] = true;
      }
      grid[i][j] = IFS2.circ_connected(min_r);
    }
  }
  
  //follow the boundary
  path.resize(0);
  closed = true;
      
      
  
  
  return true;
}