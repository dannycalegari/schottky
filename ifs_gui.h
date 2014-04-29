#ifndef __IFS_GUI__
#define __IFS_GUI__

extern "C" {
#include <X11/Xlib.h>
#include <X11/Xutil.h>
}


#include "ifs.h"

enum IFSWindowMode { LIMIT, MANDLEBROT, BOTH };

struct IFSGui {
  
  //which mode
  IFSWindowMode window_mode;
  
  //limit and mandlebrot settings
  ifs IFS;
  cpx limit_ll;
  cpx limit_ur;
  int limit_depth;
  cpx mand_ll;
  cpx mand_ur;
  int mand_depth;
  int mand_pixel_group_size;
  
  //graphics stuff
  Display* display;
  int screen;
  Window main_window;
  Pixmap limit_pixmap;
  Pixmap mand_pixmap;
  
  
  bool main_window_initialized;
  int main_window_height;
  int main_window_width;
  bool limit_pixmap_intialized;
  bool mand_pixmap_initialized;
  int limit_sidebar_size;
  int mand_sidebar_size;
  

  void launch(IFSWindowMode m = BOTH);
  void reset_and_pack_window();
  void main_loop();

};

#endif