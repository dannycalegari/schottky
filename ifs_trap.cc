//Function to check for a trap

bool ifs::find_trap() {
  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double z_restriction = abs(0.5*z-0.5)/(1.0-abs(z));
  double w_restriction = abs(0.5*w+0.5)/(1.0-abs(w));
  double smallest_radius = (z_restriction < w_restriction ? z_restriction : w_restriction);
  
  //create the set of image points
  //each image point is indexed by the binary digits
  //so 1011 means fgff, and it's a left action
  //of course we need to know the word length to parse how many g's are in front
  //starting depth will always be 8?
  std::vector<cpx> points;
  
  //starting mesh size will always be 512?
  
  
  
  return true;
}

void ifs::draw_trap() {
}
