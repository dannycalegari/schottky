#include "ifs.h"


/**************************************************************************
 * There are a bunch of big functions, so the ifs class is in 
 * multiple files.  This file contains only basic computational functions
 * ************************************************************************/

//these are the other files
#include "ifs_draw.cc"         //everything about drawing it to the xwindow
#include "ifs_interface.cc"    //everything about interacting
#include "ifs_connected.cc"    //the function to detect if the IFS is connected 
#include "ifs_trap.cc"         //functions to build a trap


ifs::ifs() {
  initialize(0,0, 1024);
}

ifs::ifs(cpx a, cpx b, int width) {
  initialize(a,b, width);
}


void ifs::initialize(cpx a, cpx b, int width){
	// initialize z to a and w to b
	z=a;
	w=b;
	az = abs(z);
	aw = abs(w);
	sync=0;
	color_ifs=true;
	chunky_ifs=false;
	disconnection_depth=false;
	draw_trap_mode = false;
	step=0.01;	// size of adjustments to z and w
	seed=0.0;	// initial seed point for IFS
	center=0.0;	// in mandelbrot mode; center of screen, size of window, and mesh of accuracy
	wind=1.0;
	mesh=1;
	depth=12;	  // depth to iterate IFS or detect connectedness to
	trap_depth = depth;  //depth to search for traps 
	drawing_width = width;
	drawing_radius = drawing_width/2;
	
	//initialize the graphics window
	//extra 580 for the right hand text
	X = XGraphics(drawing_width+580, drawing_width, 1, Point2d<float>(0,0));
};




//compute the image ball
Ball ifs::act_on_left(int index, const Ball& b) {
  int word = b.word;
  int word_len = b.word_len;
  if (index == 0) {
    return Ball( z*b.center, az*b.radius, word, word_len+1 );
  } else {
    return Ball( (w*(b.center - 1.0)) + 1.0, aw*b.radius, word | (1 << word_len), word_len+1 );
  }
}

double ifs::minimal_enclosing_radius() {
  //initialize the chunky radius to contain the whole set
  double z_restriction = abs(0.5*z-0.5)/(1.0-az);
  double w_restriction = abs(0.5-0.5*w)/(1.0-aw);
  return (z_restriction > w_restriction ? z_restriction : w_restriction);
  //cout << "z: " << z << " az: " << az << " w: " << w << " aw: " << aw << "\n"; 
  //cout << "Z restriction: " << z_restriction << "\n";
  //cout << "W restriction: " << w_restriction << "\n";
  //cout << "Computed minimal radius as " << chunky_radius << "\n";
}

cpx ifs::iterate(int index, cpx u){
	// apply generator to u
	if(index==0){
		return(u*z);
	} else {
		return(((u-1.0)*w)+1.0);
	};
};


//convert a complex to the point in the drawing
//0 is the center, and the range is [-1,1]x[-1,1]
Point2d<int> ifs::cpx_to_point(cpx w) {
  Point2d<int> p();
  p.x = drawing_radius + int(drawing_radius*w.real());
  p.y = drawing_radius + int(drawing_radius*w.imag());
  return p;
}

int ifs::cpx_to_radius(cpx w) {
  int r;
  r= int(drawing_radius*abs(w));
  return r;
}


cpx ifs::point_to_cpx(const Point2d<int>& p) {
  cpx w;
  w = cpx( double(p.x-drawing_radius)/drawing_radius,
	   double(p.y-drawing_radius)/drawing_radius );
  return w;
};



