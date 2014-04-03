/*	schottky.cc

	Version 0.01 March 15 2014

	Copyright Danny Calegari

	released under the terms of the GNU GPL, version 3.0
	see attached licence for more details
*/

// standard libraries to include

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>

using namespace std;

// preprocessor definitions

#define SMALLSCREEN		0	// for small resolutions, set to 1

#define PI 		3.14159265358979
#define TWOPI	6.28318530717959
#define cpx		complex<double>				// complex number
#define cvec 	vector<cpx >	// vector of complex numbers
#define ivec	vector<int >	// vector of integers
#define cmat	vector< vector<cpx > >		// matrix of complex

struct point{	// used for X windows graphics
	int x,y;
};

// global constants

cpx I (0.0,1.0);

#include "points.cc"
#include "graphics.cc"
#include "trap.cc"
#include "ifs.cc"
#include "connected.cc"
#include "interface.cc"
#include "ifs_trap.cc"


int main(int argc, char *argv[]){ 
	cpx z,w;
	ifs G;
	char c;
	
	G.initialize(cos(TWOPI/3.0)+I*sin(TWOPI/3.0),0.5);	// default value : Sierpinski triangle

	while(1){
		cout << "enter 'i' for IFS or 'm' for mandelbrot:";
		cin >> c;
		if(c=='i'){
			setup_graphics();
			G.mode=0;
			break;
		} else if(c=='m'){
			setup_graphics();
			G.mode=1;
			break;
		} else {
			cout << "Just type a single character ('i' or 'm') then [return].\n";
		};
	};
	G.draw();
	G.user_interface();
	
	return(0);
};
