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
  cpx I(0,1);
  int mode = 0;

  //the default start point, used when argv[1] is absent or does not parse
  cpx start(-0.58, 0.33);

  //argv[1], if given, is the start parameter -- "re im", "re+imi", a bare
  //"deg" (degrees on |s| = 1/sqrt2, the convention certify_arc and funddom
  //use), or "deg@r".  This used to be `w = atoi(argv[1])`, assigned to a
  //variable nothing else in this file reads -- dead on arrival.
  if (argc > 1) {
    cpx parsed;
    if (parse_parameter(argv[1], parsed)) {
      start = parsed;
    } else {
      std::cerr << "schottky: cannot parse start parameter \"" << argv[1]
                << "\"; expected \"re im\", \"re+imi\", \"deg\" (degrees on "
                << "|s|=1/sqrt2), or \"deg@r\".  Using the default instead.\n";
    }
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
  G.launch(BOTH, start);

  return 0;
}
