/*	schottky.cc

	Version 0.01 March 15 2014

	Copyright Danny Calegari

	released under the terms of the GNU GPL, version 3.0
	see attached licence for more details
*/

// standard libraries to include

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "ifs.h"
#include "ifs_gui.h"

#define TWOPI 6.28318530718


int main(int argc, char *argv[]) { 
  char c;
  int w = 1024;
  cpx I(0,1);
  int mode = 0;
  
  if (argc > 1) {
    w = atoi(argv[1]);
  }
  
  //std::cout << sizeof(long long int) <<"\n";
  
//  std::cout << "enter 'i' for IFS or 'm' for mandelbrot (n for new GUI):";
//  std::cin >> c;
//   if (c == 'n') {
//     IFSGui G;
//     G.launch(BOTH, cpx(-0.58, 0.33));
//   } else {
//     mode = (c=='i' ? 0 : 1);
//     ifs G(cos(TWOPI/3.0)+I*sin(TWOPI/3.0), 0.5, w, mode);      // default value : Sierpinski triangle
//     G.user_interface();
//   }
  
  IFSGui G;
  G.launch(BOTH, cpx(-0.58, 0.33));
  
  return 0;
}
