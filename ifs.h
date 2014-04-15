#ifndef __IFS__
#define __IFS__

#include <ostream>

#include "graphics.h"
#include "cpx.h"

/*************************************************************************
 * A ball is a center and radius, and it knows the word that was applied
 * to the original ball to get it
 * ***********************************************************************/
struct Ball {
  cpx center;
  double radius;
  cpx to_z;     //this points to where the center would go if we applied z
  cpx to_w;     //this points to where the center would go under w
  int word;     //the z and w word we applied to get here (highest bit is most recent application)
  int word_len; //the number of significant bits
  Ball();
  Ball(cpx c, cpx to_z, cpx to_w, double r);
  Ball(cpx c, cpx to_z, cpx to_w, double r, int w, int wl);
  int last_gen_index() const;
  bool is_disjoint(const Ball& b);
};

std::ostream& operator<<(std::ostream& os, const Ball& b);

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
		Ball act_on_left(int index, const Ball& b);
		Ball act_on_right(int index, const Ball& b);
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
                                                         double ratio_goal);
    cpx center_of_mass(const std::vector<Ball>& balls);
		bool minimal_enclosing_radius(double& r);
		
		
		//IFS drawing
		bool color_ifs;        // draw fL and gL in different colors
		bool chunky_ifs;       //draw the ifs with chunky balls
		cpx seed;              // initial seed for IFS in ifs mode
		double step;           // step size to adjust z or w in ifs mode
		
		void draw_dots(int depth, cpx u);
		void draw_color_dots(int d, cpx u, long color);
		void draw_color_chunky_dots(int d, cpx u, long color, double radius);
		void draw_limit_set();
		
		
		//mandlebrot mode
		cpx center;            // center of screen in mandelbrot mode
		double wind;           // size of window in mandelbrot mode
		int mesh;              // size of mesh in mandelbrot mode in pixels
		int mode;              // draw mode: 0 for limit set, 1 for mandelbrot set
		bool disconnection_depth; // whether to draw the depth for disconnected sets
		
		void zoom(const Point2d<int>& p);
		void draw_mandelbrot_set();
		
		
		//connectedness testing
		int depth;             // depth to recurse to draw or detect connectedness
		int exit_depth;        // size of tree to detect connectedness
		
		bool old_circles_intersect(cpx c1, cpx a1, cpx c2, cpx a2, double R, int d); // recursive test; do circles intersect?
                bool circles_intersect(cpx center_1, cpx z_img_1, cpx w_img_1, double R1,
                                       cpx center_2, cpx z_img_2, cpx w_img_2, double R2, int d);
		bool circ_connected(double r=-1);                                        // circle algorithm to test for connectedness
		
		
		//trap construction
		bool draw_trap_mode;      //whether to check for a trap and draw it in limit set mode
		int trap_depth;           //maximal depth to look for traps
		//Trap current_trap;        //the last trap we created
		
		bool find_trap_given_balls(const std::vector<Ball>& initial_balls, 
                                           int max_refinements,
                                           int max_pixels,
                                           int verbose);
		bool find_trap(double* epsilon, int verbose=0);
		bool find_traps_along_loop(const std::vector<cpx>& loop, 
		                           bool draw_it, 
		                           int verbose);
		void draw_trap();	
		
		
		//Main interface and drawing functions
		XGraphics X;
		int drawing_width;
		int drawing_radius; //drawing_width/2
		
		Point2d<int> cpx_to_point(cpx w); //this is for ifs mode
		int cpx_to_radius(cpx w);
		cpx point_to_cpx(const Point2d<int>& p);  //this is for mandlebrot mode
		Point2d<int> cpx_to_point_mandlebrot(cpx w); //this is for mandlebrot mode
		void draw();
		void user_interface();
		void input_loop(std::vector<cpx>& loop);
};

#endif