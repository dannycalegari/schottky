//Functions to check for a trap


//this records the trap computation grid
struct TrapGrid {
  int width; //it's a square
  
  std::vector<int> f_ball_hit;
  std::vector<int> g_ball_hit;
  TrapGrid() {
    width = 0;
    f_ball_hit = std::vector<int>(0);
    g_ball_hit = std::vector<int>(0);
  }
};


//this records a ball
struct Ball {
  Ball() { 
    center = 0.5;
    radius = 1.0;
  }
  Ball(cpx c, double r) {
    center = c;
    radius = r;
  }
  cpx center;
  double radius;
};

//compute the image ball
//(left action)
Ball ifs::act_on_left(int index, const Ball& b) {
  if (index == 0) {
    return Ball( z*b.center, az*b.radius );
  } else {
    return Ball( (w*(b.center - 1.0)) + 1.0, aw*b.radius );
  }
}

//take a list of balls and compute one more level of depth 
//(on the left)
void ifs::compute_next_ball_depth(std::vector<Ball>& balls, int current_depth) {
  std::vector<Ball> balls_temp;
  balls_temp.swap(balls);
  int L_cur = 1<<current_depth;
  int L_new = 1<<(current_depth+1);
  balls.resize(L_new);
  for (int j=0; j<L_cur; ++j) {
    //for each j, we want to append on the left both a 0 and a 1
    //(the highest bit position is the last function we applied)
    Ball parent_ball = balls_temp[j];
    balls[j] = act_on_left(0, parent_ball);
    balls[(j | (1<<current_depth))] = act_on_left(1, parent_ball);
  }
}

//take a seed ball and compute all the image balls
//create a list of image points
//each image point is indexed by the binary digits
//so 1011 means fgff, and it's a left action
//of course we need to know the word length to parse how many g's are in front
void ifs::compute_balls(std::vector<Ball>& balls, const Ball& ball_seed, int compute_depth) {
  std::vector<Ball> balls_temp;
  balls.resize(2);
  balls[0] = act_on_left(0, ball_seed);
  balls[1] = act_on_left(1, ball_seed);
  for (int i=1; i<compute_depth; ++i) {
    compute_next_ball_depth(balls, i);
  }
}


bool ifs::find_trap() {
  //find the radius of the smallest closed ball about 1/2 which 
  //is mapped inside itself under both f and g
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
	double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  double min_initial_radius = (z_restriction > w_restriction 
                                             ? z_restriction 
                                             : w_restriction);
  
  
  //starting depth will always be 8?
  int current_depth = 8;
  std::vector<Ball> balls;
  compute_balls(balls, Ball(0.5, min_initial_radius), current_depth);
  
  //now we need to effectively increase the radius of the 
  //*previous* step until the balls of the *current* step 
  //form a connected set for both f and g
  //we'll start by doubling the radius (it's probably pretty small)
  double more_dramatic_action = (az < aw ? az : aw);
  double smallest_previous_radius = pow( more_dramatic_action, current_depth-1 ) * min_initial_radius;
  //double epsilon = 
  
  
  
  
  
  //figure out a radius boost such that adding this much on to the 
  //radius of all image balls makes the f ball union connected and 
  //the g ball union also connected
  double larger_action = (az < aw ? az : aw);
  double max_current_radius = pow( larger_action, current_depth ) * min_initial_radius;
  //this needs to be some positive number; this is kind of random
  double current_radius_boost = 0.1*max_current_radius;
  //while (true) {
    //fill out the grid with the current balls and radius boost
    
  //}
  
  //starting mesh size will always be 512?
  
  
  
  return true;
}

void ifs::draw_trap() {
}















