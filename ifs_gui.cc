#include "ifs_gui.h"



void IFSGui::reset_and_pack_window() {
  //destroy the main window, if it exists
  if (main_window_initialized) {
    XDestroyWindow(main_window);
  }
  //figure out how big to make the window (here ss = sidebar size)
  //for each of limit and mandlebrot, we need height x+ss and width x+ss,
  //if if it's just one of them, we can set the height to be the display height
  //minus 200 (say), as long as x+ss is smaller than the width
  //otherwise, we can make it x+ss, as long as 2x+ss < width
  int display_width = XDisplayWidth(display, screen);
  int display_height = XDisplayHeight(display, screen);
  int ss = (limit_sidebar_size > mand_sidebar_size ? limit_sidebar_size : mand_sidebar_size);
  int width_rest = (window_mode == BOTH ? (display_width-2*ss)/2 :
                                                 display_width - ss );
  int height_rest = display_height - 200;
  int x = (width_rest > height_rest ? height_rest : width_rest);
  
  if (m == MANDLEBROT) {
    window_height = x + 200;
    window_width = x + mand_sidebar_size;
  } else if (m == LIMIT) {
    window_height = x + 200;
    window_width = x + limit_sidebar_size;
  } else {
    window_height = x + 200;
    window_width = 2*x + mand_sidebar_size + limit_sidebar_size;
  }
  
  //create the window
  main_window = XCreateSimpleWindow(display, 
                                    RootWindow(display, screen), 20, 20,
                                    window_width, window_height,
                                    BlackPixel(display, screen), WhitePixel(display, screen));
  XSelectInput(display, main_window, ExposureMask |
                                     PointerMotionMask |
                                     KeyPressMask |
                                     ButtonPressMask |
                                     ButtonReleaseMask );
  XMapWindow(display, main_window);
  //wait until the window is actually mapped
  while (true) {  
    XEvent e;
    XNextEvent(display, &e);
    if (e.type == Expose) break;
  }
  main_window_initialized = true;
  
  //create the limit and mandlebrot stuff
  limit_widget = Widget(x,x);
  mand_widget = Widget(x,x);
  
  //create the buttons and stuff
  widgets.resize(0);
  
  
  
  
  
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
  main_window_initialized = false;
  
  //reset (set) the window
  limit_sidebar_size = 100;
  mand_sidebar_size = 100;
  reset_and_pack_window();
  
  //go for it
  main_loop();
  
  //clean up at the end
  XCloseDisplay(display);
  
}