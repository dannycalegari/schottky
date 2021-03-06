#ifndef __IFS__
#define __IFS__

#include <ostream>
#include <bitset>

#include "graphics.h"
#include "cpx.h"

#include "point.h"

/*************************************************************************
 * convex hull computational geometry stuff
 * ***********************************************************************/
//this does dot product as if they were vectors
double cpx_dot(cpx a, cpx b);

//return a complex number b such that {b,a} is an oriented basis for R^2
cpx perp_to(cpx a);

//the halfspace is all x such that v.(x-p) <= 0 
struct halfspace {
  cpx v;
  cpx p;
  halfspace() {v = p = 0.0;}
  halfspace(cpx V, cpx P) {
    v = V;
    p = P;
  }
  bool contains(cpx x) {
    return cpx_dot(v, x-p) <= 0;
  }
  bool strictly_contains(cpx x) {
    return cpx_dot(v, x-p) < 0;
  }
  double val(cpx x) {
    return cpx_dot(v, x-p);
  }
};

//returns the halfspace which is on the left when walking x1->x2
halfspace halfspace_on_left(cpx x1, cpx x2);

//returns a list of indices of the points that are in the convex hull
void convex_hull(std::vector<int>& ch, 
                 const std::vector<cpx>& X);


/*************************************************************************
 * A ball is a center and radius, and it knows the word that was applied
 * to the original ball to get it
 * ***********************************************************************/
struct Ball {
  cpx center;
  double radius;
  cpx to_z;     //this points to where the center would go if we applied z
  cpx to_w;     //this points to where the center would go under w
  std::bitset<64> word;     //the z and w word we applied to get here (highest bit is most recent application)
  int word_len; //the number of significant bits
  Ball();
  Ball(cpx c, double r);
  Ball(cpx c, cpx to_z, cpx to_w, double r);
  Ball(cpx c, cpx to_z, cpx to_w, double r, const std::bitset<64>& w, int wl);
  int last_gen_index() const;
  bool is_disjoint(const Ball& b);
  bool is_disjoint(const cpx& ll, const cpx& ur);
  bool is_contained(const cpx& ll, const cpx& ur);
};

std::ostream& operator<<(std::ostream& os, const Ball& b);

void ball_convex_hull(std::vector<int>& ch,
                      std::vector<cpx>& boundary_points,
                      std::vector<halfspace>& H,
                      const std::vector<Ball>& balls);

void box_containing_points(const std::vector<cpx>& points, cpx& ll, cpx& ur);



/**************************************************************************
 * a bit word
 **************************************************************************/
struct Bitword {
  std::bitset<64> w;
  int len;
  Bitword () {
    w = std::bitset<64>(0);
    len = 0;
  }
  Bitword (const std::bitset<64>& W, int L) {
    w = W;
    len = L;
  }
  bool operator<(const Bitword& b) const;
  Bitword prefix(int n) const;
  Bitword suffix(int n) const;
  std::string str() const;
  int reverse_get(int n) const;
  cpx apply(cpx z, cpx x) const;
  Bitword pow(int n) const;
  Bitword append(int n, int L) const;
  bool operator==(const Bitword& b) const;
  bool operator!=(const Bitword& b) const;
};

std::ostream& operator<<(std::ostream& os, const Bitword& b);


/**************************************************************************
 * data defining an IFS
 * IFS is generated by two similarities f and g; f scales by z and 
 * is centered at 0, g scales by w and is centered at 1.  When referred to
 * by number, f is 0 and g is 1
 * ************************************************************************/
class ifs{
  public:
    //main data
    cpx z,w;               // parameters for IFS
    double az, aw;         //absolute values of them
    int sync;              // sync=0 w is arbitrary, sync=1 w=z, sync=2 w=-z
    
    
    // these are the construction functions
    ifs();
    ifs(cpx a, cpx b, int width, int mode);
    void initialize(cpx a, cpx b, int width, int mode);
    void reinitialize(cpx a, cpx b);
    
    
    //general computation
    
    cpx iterate(int index, cpx u);
    Ball act_on_left(int index, const Ball& b) const;
    Ball act_on_right(int index, const Ball& b) const;
    cpx apply_bitword(const Bitword& b, cpx x) const;
    void compute_next_ball_depth(std::vector<Ball>& balls, int current_depth);
    void compute_balls(std::vector<Ball>& balls, const Ball& ball_seed, int compute_depth);
    void compute_next_ball_depth_right(std::vector<Ball>& balls, int current_depth);
    void compute_balls_right(std::vector<Ball>& balls, const Ball& ball_seed, int compute_depth);
    void box_containing_balls(const std::vector<Ball>& balls, cpx& ll, cpx& ur);
    void refine_balls_into_box(std::vector<Ball>& balls, const cpx& ll, const cpx& ur);
    bool is_ball_disjoint(const Ball& b, const cpx& ll, const cpx& ur);
    void find_close_images_with_distinct_first_letters(const Ball& b, int length, Ball& zb, Ball& wb);
    void find_aligned_images_with_distinct_first_letters(const Ball& initial_ball, 
                                                         cpx p1, cpx p2, int search_depth,
                                                         Ball& zb, Ball& wb,
                                                         double ratio_goal, double ratio_lower_limit);
    cpx center_of_mass(const std::vector<Ball>& balls);
    bool minimal_enclosing_radius(double& r);
    double distance_from_balls(cpx p, const std::vector<Ball>& balls);
    double when_ray_hits_ball(cpx p, cpx v, const Ball& b);
    double when_ray_hits_ball(cpx p, cpx v, const std::vector<Ball>& balls);
    void find_closest_uv_words(std::vector<std::pair<Bitword,Bitword> >& words, 
                               int uv_depth,
                               double last_step_tolerance=0.0,
                               int list_size_max=-1);
    void find_closest_uv_words_along_path(const std::vector<cpx>& path, 
                                          bool closed_path, 
                                          int word_len);
    bool close_to_set_C(int n_depth, double epsilon);
    void compute_uv_graph(std::vector<Point3d<int> >& uv_graph, 
                          std::vector<Ball>& balls, 
                          int uv_depth, 
                          int verbose);
    void set_params(cpx Z, cpx W);    
    void half_balls(std::vector<Bitword>& half_words, 
                    int first_depth, int second_depth);
    std::vector<int> coefficient_list(const Bitword& u); 
    void word_deriv(const Bitword& u, const cpx& z0, cpx& deriv, double& err);
    
    
    //IFS drawing
    bool color_ifs;        // draw fL and gL in different colors
    bool chunky_ifs;       //draw the ifs with chunky balls
    cpx seed;              // initial seed for IFS in ifs mode
    double step;           // step size to adjust z or w in ifs mode
    bool find_close_uv_words; //whether to find close uv words
    
    void draw_dots(int depth, cpx u);
    void draw_color_dots(int d, cpx u, long color);
    void draw_color_chunky_dots(int d, cpx u, long color, double radius);
    void draw_limit_set();
    
    //IFS coordinate computation
    bool compute_coordinates(double* theta, double* lambda, int n_depth);
    
    //mandelbrot mode
    cpx center;            // center of screen in mandelbrot mode
    double wind;           // size of window in mandelbrot mode
    int mesh;              // size of mesh in mandelbrot mode in pixels
    int mode;              // draw mode: 0 for limit set, 1 for mandelbrot set
    bool disconnection_depth; // whether to draw the depth for disconnected sets
    bool draw_contains_half; //whether to color the points that contain 1/2
    
    void zoom(const Point2d<int>& p);
    void draw_mandelbrot_set();
		
		
    //connectedness testing
    int depth;             // depth to recurse to draw or detect connectedness
    int exit_depth;        // size of tree to detect connectedness
    
    bool old_circles_intersect(cpx c1, cpx a1, cpx c2, cpx a2, double R, int d); // recursive test; do circles intersect?
    bool circles_intersect(cpx center_1, cpx z_img_1, cpx w_img_1, double R1,
                            cpx center_2, cpx z_img_2, cpx w_img_2, double R2, int d);
    bool circ_connected(double r=-1);                                        // circle algorithm to test for connectedness
    bool is_connected(int d, int& difficulty); //automates it
    bool contains_point(cpx pt, double r=-1);
    bool contains_point_recurse(const cpx& pt, const Ball& b, int d);
		bool contains_half(int d, int& difficulty);
		
    //trap construction
    bool draw_trap_mode;      //whether to check for a trap and draw it in limit set mode
    int trap_depth;           //maximal depth to look for traps
    //Trap current_trap;        //the last trap we created
    
    bool find_trap_given_balls_old(const std::vector<Ball>& initial_balls, 
                                int max_refinements,
                                int max_pixels,
                                bool far_trap_points,
                                double* minimum_trap_distance,
                                int verbose);
    bool find_trap_given_balls(const std::vector<Ball>& balls,
                               int max_pixels,
                               double* min_trap_distance,
                               int verbose);
    bool find_trap(int max_uv_depth, int max_n_depth, int max_pixels, double Cz, double* epsilon, double* difficulty, int verbose=0);
    bool find_traps_along_loop(const std::vector<cpx>& loop, 
                                bool draw_it, 
                                int verbose);
    void draw_trap();	
    
    //trap-like vectors
    void get_TLB_constants(cpx ll, cpx ur, double& K, double& C, double& A, double& Z);
    bool trap_like_balls(std::vector<Ball>& TLB, 
                          double initial_radius, 
                          double radius_subtraction,
                          int n_depth, 
                          int verbose);
    void trap_like_balls_from_balls(std::vector<Ball>& TLB, 
                        int num_TL_balls, 
                        int num_ball_trials,
                        const std::vector<Ball>& balls,
                        double radius_subtraction,
                        int verbose);
    bool TLB_for_region(std::vector<Ball>& TLB,
                        cpx ll, cpx ur, int n_depth, double* TLB_C, double* TLB_Z, int verbose);
    int check_TLB_and_uv_words(const std::vector<Ball>& TLB, 
                               const std::vector<std::pair<Bitword,Bitword> >& words);
    int check_TLB(const std::vector<Ball>& TLB, 
                  double* TLB_C, double* TLB_Z, 
                  double& trap_radius, 
                  std::vector<std::pair<Bitword,Bitword> >* trap_w, 
                  int uv_depth) ;
    bool find_TLB_along_loop(const std::vector<cpx>& loop, 
                             bool draw_it, 
                             int verbose);
    //hardcoded spiral traps
    int check_limit_TLB(const std::vector<Ball>& TLB, 
                        double* TLB_C, double* TLB_Z, 
                        double& trap_radius, 
                        std::vector<std::pair<Bitword,Bitword> >* trap_w,
                        int n_limit);
    int check_limit_TLB_recursive(const std::vector<Ball>& TLB, 
                                  double* TLB_C, double* TLB_Z, 
                                  double& trap_radius, 
                                  std::vector<std::pair<Bitword,Bitword> >* trap_w, 
                                  int uv_depth);
    bool find_trap_like_vectors;
    
    
    //set B
    void deriv_bounds_around_ball(const Bitword& u, 
                                  cpx z0, 
                                  double r, 
                                  double& max,
                                  double& min);
    bool certify_set_B_point(const Bitword& u, bool certify_all, double& within);
    std::vector<Bitword> get_half_balls_along_path(const std::vector<cpx>& path,
                                                   int d,
                                                   int verbose);
    std::vector<Bitword> get_certified_half_balls_along_path(const std::vector<cpx>& path,
                                                             int d,
                                                             int verbose);
    bool certify_set_B_path(const std::vector<cpx>& path, int initial_depth, int verbose);
    std::vector<Ball> subdivide_half_prefix(const Bitword& u, 
                                             cpx start_z,
                                             int d,
                                             cpx ll, cpx ur);
    
    
    //balls
    void draw_set_B_balls(const std::vector<Bitword>& balls, 
                          cpx initial_point, 
                          int d, 
                         int verbose);
    cpx solve_for_half(const Bitword& u, cpx start, double tol);
    
    //hole boundary finding functions
    bool hole_boundary_containing_point(std::vector<cpx>& path, bool& closed, 
                                        cpx p, int verbose);
    bool hole_boundary_containing_point_from_grid(std::vector<cpx>& path, bool& closed, 
                                                  const std::vector<std::vector<bool> >& grid,
                                                  const Point2d<int>& marked_pixel, 
                                                  cpx ll, cpx ur, int verbose) ;
    
    //picture creation
    void draw_ifs_to_array(std::vector<std::vector<Point3d<unsigned char> > >& bmp, 
                           const cpx& region_ll, const cpx& region_ur, 
                           int depth); 
    void draw_mand_to_array(std::vector<std::vector<Point3d<unsigned char> > >&bmp,
                           const cpx& region_ll, const cpx& region_ur,
                           int connected_depth, int contains_half_depth); 
    void draw_mand_to_array_radial(std::vector<std::vector<Point3d<unsigned char> > >&bmp,
                                    int num_pixels, //height
                                    double radius_start, double radius_end,
                                    double theta_start, double theta_end,
                                    int connected_depth, int contains_half_depth);
    //Main interface and drawing functions
    XGraphics X;
    int drawing_width;
    int drawing_radius; //drawing_width/2
    
    Point2d<int> cpx_to_point(cpx w); //this is for ifs mode
    int cpx_to_radius(cpx w);
    cpx point_to_cpx(const Point2d<int>& p);  //this is for mandelbrot mode
    Point2d<int> cpx_to_point_mandelbrot(cpx w); //this is for mandelbrot mode
    void draw();
    void user_interface();
    void input_loop(std::vector<cpx>& loop);
};



/****************************************************************************
 * A more general ball with less information
 ****************************************************************************/
struct nBall {
  cpx center;
  double radius;
  cpx one; // this is the vector which points where 1 is in the original ball
  
  nBall();
  nBall(cpx c, double r, cpx o);
  bool is_disjoint(const nBall& other) const;
  bool is_disjoint(const cpx& ll, const cpx& ur) const;
  bool is_contained(const cpx& ll, const cpx& ur) const;
};

std::ostream& operator<<(std::ostream& os, const nBall& b);

//the boolean means whether the ball is completely contained in the window
//the ints are (last gen, depth)
struct nBall_stuff { 
    bool contained;
    int last_gen;
    int depth;
    nBall ball;
    nBall_stuff() {}
    nBall_stuff(bool c, int l, int d, nBall b) { 
      contained = c; last_gen = l; depth = d; ball = b;
    }
  };
  
std::ostream& operator<<(std::ostream& os, const nBall_stuff& b);
/***************************************************************************
 * a more general IFS
 ***************************************************************************/
struct nIFS {
  cpx c;
  std::vector<cpx> centers;
  
  //initializing
  nIFS();
  nIFS(int n, cpx C);
  void set_param(cpx C);
  void set_centers(const std::vector<cpx>& new_centers);
  void set_std_centers(int n);
  
  //return the minimal radius of a ball sent inside itself under all 
  //functions
  double minimal_initial_radius(); 
  
  //acting
  nBall act_on_left(int i, const nBall& b);
  nBall act_on_right(int i, const nBall& b);
  
  //computation
  bool is_connected(int depth, int& difficulty);
};

/**************************************************************************
 * A most general ball
 **************************************************************************/
struct gBall {
  cpx center;
  double radius;
  std::vector<cpx> ifs_centers; //this records where the IFS centers 
                               //are upon subdividing
  
  gBall();
  gBall(cpx c, double r, const std::vector<cpx>& i_c);
  bool is_disjoint(const gBall& other) const;
  bool is_disjoint(const cpx& ll, const cpx& ur) const;
  bool is_contained(const cpx& ll, const cpx& ur) const;
};

std::ostream& operator<<(std::ostream& os, const gBall& b);

//the boolean means whether the ball is completely contained in the window
//the ints are (last gen, depth)
struct gBall_stuff { 
  bool contained;
  int last_gen;
  int depth;
  gBall ball;
  gBall_stuff() {}
  gBall_stuff(bool c, int l, int d, gBall b) { 
    contained = c; last_gen = l; depth = d; ball = b;
  }
};
  
std::ostream& operator<<(std::ostream& os, const gBall_stuff& b);

/***************************************************************************
 * a most general IFS
 ***************************************************************************/
struct gIFS {
  std::vector<cpx> factors;
  std::vector<cpx> centers;
  
  //initializing
  gIFS();
  gIFS(const std::vector<cpx>& facs, const std::vector<cpx>& cents);
  void set_params(const std::vector<cpx>& new_factors);
  void set_centers(const std::vector<cpx>& new_centers);
  
  //return a ball centered at the average of the centers and sent 
  //inside itself under all maps
  bool minimal_ball(gBall& b);
  
  //acting
  gBall act_on_left(int i, const gBall& b);
  gBall act_on_right(int i, const gBall& b);
  
  //computation
  bool is_connected(int depth, int& difficulty);
};
  
/***************************************************************************
 * a 2d affine map 
 * *************************************************************************/
struct LinearMap {
  double a,b,c,d;
  LinearMap();
  LinearMap(double A, double B, double C, double D);
  LinearMap operator*(const LinearMap& other) const;
  Point2d<double> operator()(const Point2d<double>& X) const;
};

struct AffineMap {
  LinearMap A;
  Point2d<double> t;
  AffineMap();
  AffineMap(const LinearMap& a, const Point2d<double>& T);
  AffineMap(double a, double b, double c, double d, double x, double y);
  AffineMap operator*(const AffineMap& other) const;
  Point2d<double> operator()(const Point2d<double>& X) const;
};


bool is_point_contained(const cpx& c, const cpx& ll, const cpx& ur);

struct Box_Stuff {
  bool contained;
  int last_gen;
  int word;
  int depth;
  std::vector< Point2d<double> > box;
  Box_Stuff() {}
  Box_Stuff(bool C, int L, int W, int D, const std::vector< Point2d<double> >& B) {
    contained = C;
    last_gen = L;
    word = W;
    depth = D;
    box = B;
  }
  bool is_disjoint(cpx ll, cpx ur) {
    //check if the box is separated from the square by a square side
    if ( (box[0].x < ll.real() && box[1].x < ll.real() && box[2].x < ll.real() && box[3].x < ll.real()) ||
         (box[0].x > ur.real() && box[1].x > ur.real() && box[2].x > ur.real() && box[3].x > ur.real()) ||
         (box[0].y < ll.imag() && box[1].y < ll.imag() && box[2].y < ll.imag() && box[3].y < ll.imag()) ||
         (box[0].y > ur.imag() && box[1].y > ur.imag() && box[2].y > ur.imag() && box[3].y > ur.imag()) ) {
      return true;
    }
    //check if the box is separated from the square by a box side
    std::vector<Point2d<double> > square(4);
    square[0] = Point2d<double>(ll.real(), ll.imag());
    square[1] = Point2d<double>(ur.real(), ll.imag());
    square[2] = Point2d<double>(ur.real(), ur.imag());
    square[3] = Point2d<double>(ll.real(), ur.imag());
    for (int i=0; i<(int)box.size(); ++i) {
      Point2d<double> v = box[(i+1)%4] - box[i];
      Point2d<double> vp(-v.y, v.x);
      if ( dot(vp, square[0]-box[i]) < 0 && dot(vp, square[1]-box[i]) < 0 &&
           dot(vp, square[2]-box[i]) < 0 && dot(vp, square[3]-box[i]) < 0    ) {
        return true;
      }
    }
    return false;
  }
  bool is_contained(cpx ll, cpx ur) {
    for (int i=0; i<(int)box.size(); ++i) {
      if (!is_point_contained(cpx(box[i].x, box[i].y), ll, ur)) {
        return false;
      }
    }
    return true;
  }
};

Point4d<double> point_as_weighted_average_in_box(const Point2d<double>& point_target, 
                                                 const Point2d<double>& center, 
                                                 double radius);

/***************************************************************************
 * a 2d affine IFS      
 * *************************************************************************/
struct ifs2d {
  std::vector<AffineMap> gens;
  AffineMap semigroup_element(const std::vector<int>& gen_word);
  AffineMap semigroup_element(int list, int n);
};






























#endif