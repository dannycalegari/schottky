#include "ifs_gui.h"



void IFSGui::reset_and_pack_window() {

}

void IFSGui::main_loop() {

}




void IFSGui::launch(IFSWindowMode m) {
  
  //set the initial window mode
  window_mode = m;
  
  //set the initial settings for limit and mandlebrot
  limit_ll = cpx(-1, -1.5);
  limit_ur = cpx(2, 1.5);
  limit_depth = 12;
  mand_ll = cpx(-1,-1);
  mand_ur = cpx(1,1);
  mand_depth = 12;
  mand_pixel_group_size = 4;
  
  
  //set up the graphics
  display = XOpenDisplay(NULL);
  screen = DefaultScreen(display);
  if (display == NULL) {
    std::cout << "Failed to open display\n";
    return;
  }
  initialized = false;
  
  //reset (set) the window
  reset_and_pack_window();
  
  //go for it
  main_loop();
  
  //clean up at the end
  XCloseDisplay(display);
  
}