#include "ifs_gui.h"










/***************************************************************************
 * widget functions
 ***************************************************************************/
Widget::Widget() {
  //don't do anything
}

WidgetDraw::WidgetDraw(IFSGui* i, int w, int h) {
  width = w;
  height = h;
  ifsg = i;
}

WidgetButton::WidgetButton(IFSGui* i, const std::string& t, int w, int h, void (IFSGui::*f)()) {
  ifsg = i;
  width = w;
  height = h;
  text = t;
  clicker = f;
}

WidgetText::WidgetText(IFSGui* i, const std::string& t, int w, int h) {
  ifsg = i;
  text = t;
  width = w; 
  height = h;
}

WidgetCheck::WidgetCheck(IFSGui* i, const std::string& t, int w, int h, bool c, void (IFSGui::*f)()) {
  width = w;
  height = h;
  ifsg = i;
  text = t;
  checked = c;
  clicker = f;
} 

WidgetLeftArrow::WidgetLeftArrow(IFSGui* i, int w, int h, void (IFSGui::*f)()) {
  ifsg = i;
  width = w;
  height = h;
  clicker = f;
}

WidgeRightArrow::WidgetRightArrow(IFSGui* i, int w, int h, void (IFSGui::*f)()) {
  ifsg = i;
  width = w;
  height = h;
  clicker = f;
}

















/****************************************************************************
 * signal handlers
 ****************************************************************************/
 
 
 
void IFSGui::S_limit_increase_depth() {
}

void IFSGui::S_limit_decrease_depth() {
}

void IFSGui::S_limit_switch_chunky() {
}

void IFSGui::S_limit_zoom_in() {
}

void IFSGui::S_limit_zoom_out() {
}

void IFSGui::S_limit_recenter() {
}


/****************************************************************************
 * main gui functions
 ****************************************************************************/

void IFSGui::pack_widget_upper_right(const Widget* w1, Widget* w2) {
  //figure out where it can go
  int desired_x;
  if (w1 != NULL) {
    desired_x = w1->ul.x + w1->width;
  } else {
    desired_x = 0;
  }
  
  //go through and check the other widgets to see how 
  //far down they obstruct this one
  int greatest_y_obstruction = 0;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w1) continue;
    if (widgets[i]->ul.x == desired_x && 
        widgets[i]->ul.y + widgets[i]->height > greatest_y_obstruction) {
      greatest_y_obstruction = widgets[i]->ul.y + widgets[i]->height;
    }
  }
  if (greatest_y_obstruction + w2->height > main_window_height) {
    std::cout << "Cannot pack widget -- too tall\n";
    return;
  }
  int y = greatest_y_obstruction;
  
  //now determine whether we have to shove it over to make room
  int greatest_x_obstruction = desired_x;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i]->ul.y + widgets[i]->height > y && 
        widgets[i]->ul.y < y + w2->height && 
        widgets[i]->ul.x + widgets[i]->width > greatest_x_obstruction) {
      greatest_x_obstruction = widgets[i]->ul.x + widgets[i]->width;
    }
  }
  int x = greatest_x_obstruction;
  
  //find the position
  w2->ul = Point2d<int>(x, y);
  
  //record it in the list of widgets
  widgets.push_back(w2);
}




void IFSGui::reset_and_pack_window() {
  //destroy the main window, if it exists
  if (main_window_initialized) {
    XDestroyWindow(display, main_window);
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
  
  if (window_mode == MANDLEBROT) {
    main_window_height = x + 200;
    main_window_width = x + mand_sidebar_size;
  } else if (window_mode == LIMIT) {
    main_window_height = x + 200;
    main_window_width = x + limit_sidebar_size;
  } else {
    main_window_height = x + 200;
    main_window_width = 2*x + mand_sidebar_size + limit_sidebar_size;
  }
  
  //create the window
  main_window = XCreateSimpleWindow(display, 
                                    RootWindow(display, screen), 20, 20,
                                    main_window_width, main_window_height, 4,
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
  
  //create the buttons and stuff and pack them
  widgets.resize(0);
  
  //stuff that shows up everywhere
  W_switch_to_limit = WidgetButton(this, "Switch to limit", 50, 20);
  W_switch_to_mandlebrot = WidgetButton(this, "Switch to mandlebrot", 50,20);
  W_switch_to_combined = WidgetButton(this, "Switch to combined", 50, 20);
  
  W_point_title = WidgetText(this, "Highlighted IFS options:", 100, 20);
  W_point_connected_check = WidgetCheck(this, "Connectedness", 100, 20, (point_connected ? 1 : 0));
  W_point_contains_half_check = WidgetCheck(this, "Contains 1/2", 100, 20, (point_contains_half ? 1 : 0));
  W_point_trap_check = WidgetCheck(this, "Find trap", 50, 20, (point_trap ? 1 : 0));
  W_point_uv_word_check = WidgetCheck(this, "Find uv words", 50, 20, (point_uv_words ? 1 : 0));
  
  
  if (window_mode != MANDLEBROT) {
    W_limit_plot = WidgetDraw(this, x,x);
    W_limit_depth_title = WidgetText(this, "Depth: ", 50, 20);
    W_limit_depth_leftarrow = WidgetLeftArrow(this, 20,20);
    std::stringstream T; T.str(""); T << limit_depth;
    W_limit_depth_label = WidgetText(this, T.str().c_str(), 30, 20);
    W_limit_depth_rightarrow = WidgetRightArrow(this, 20,20);
    W_limit_chunky_title = WidgetText(this, "Chunky: ", 50, 20);
    W_limit_chunky_on = WidgetCheck(this, "on", 20, 20, (limit_chunky ? 1 : 0));
    W_limit_chunky_off = WidgetCheck(this, "off", 30, 20, (limit_chunky ? 0 : 1));
    W_limit_zoom_title = WidgetText(this, "Zoom: ", 40, 20);
    W_limit_zoom_in = WidgetButton(this, "in", 30, 20);
    W_limit_zoom_out = WidgetButton(this, "out", 30, 20);
    W_limit_center_title = WidgetText(this, "(Click to center)", 40, 20);
    
    pack_widget_upper_right(NULL, &W_limit_plot);
    if (window_mode == LIMIT) {
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_mandlebrot);
      pack_widget_upper_right(&W_switch_to_mandlebrot, &W_switch_to_combined);
    }
    pack_widget_upper_right(&W_limit_plot, &W_limit_depth_title);
    pack_widget_upper_right(&W_limit_depth_title, &W_limit_depth_leftarrow);
    pack_widget_upper_right(&W_limit_depth_leftarrow, &W_limit_depth_label);
    pack_widget_upper_right(&W_limit_depth_label, &W_limit_depth_rightarrow);
    pack_widget_upper_right(&W_limit_plot, &W_limit_chunky_title);
    pack_widget_upper_right(&W_limit_chunky_title, &W_limit_chunky_on);
    pack_widget_upper_right(&W_limit_chunky_on, &W_limit_chunky_off);
    pack_widget_upper_right(&W_limit_plot, &W_limit_zoom_title);
    pack_widget_upper_right(&W_limit_zoom_title, &W_limit_zoom_in);
    pack_widget_upper_right(&W_limit_zoom_in, &W_limit_zoom_out);
    pack_widget_upper_right(&W_limit_plot, &W_limit_center_title);
    
    if (window_mode == LIMIT) {
      pack_widget_upper_right(&W_limit_plot, &W_point_title);
      pack_widget_upper_right(&W_limit_plot, &W_point_connected_check);
      pack_widget_upper_right(&W_limit_plot, &W_point_contains_half_check);
      pack_widget_upper_right(&W_limit_plot, &W_point_trap_check);
      pack_widget_upper_right(&W_limit_plot, &W_point_uv_word_check);
    }
  }
  
  if (window_mode != LIMIT) {
    W_mand_plot = WidgetDraw(this, x,x);
    W_mand_connected_check = WidgetCheck(this, "Connectedness - depth:", 100, 20, (mand_connected ? 1 : 0));
    W_mand_connected_depth_leftarrow = WidgetLeftArrow(this, 20,20);
    std::stringstream T;  T.str("");  T << mand_connected_depth;
    W_mand_connected_depth_label = WidgetText(this, T.str().c_str(), 30, 20);
    W_mand_connected_depth_rightarrow = WidgetRightArrow(this, 20,20);
    W_mand_contains_half_check = WidgetCheck(this, "Contains 1/2 - depth:", 100, 20, (mand_contains_half ? 1 : 0));
    W_mand_contains_half_depth_leftarrow = WidgetLeftArrow(this, 20,20);
    T.str("");  T << mand_contains_half_depth;
    W_mand_contains_half_depth_label = WidgetText(this, T.str().c_str(), 30, 20);
    W_mand_contains_half_depth_rightarrow = WidgetRightArrow(this, 20,20);
    W_mand_trap_check = WidgetCheck(this, "Traps - depth:", 100, 20, (mand_trap ? 1 : 0));;
    W_mand_trap_depth_leftarrow = WidgetLeftArrow(this, 20,20);
    T.str("");  T << mand_trap_depth;
    W_mand_trap_depth_label = WidgetText(this, T.str().c_str(), 30, 20);
    W_mand_trap_depth_rightarrow = WidgetRightArrow(this, 20,20);
    
    if (window_mode == MANDLEBROT) {
      pack_widget_upper_right(NULL, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_switch_to_limit, &W_switch_to_combined);
    } else {
      pack_widget_upper_right(&W_limit_depth_title, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_switch_to_limit, &W_switch_to_mandlebrot);
    }
    pack_widget_upper_right(&W_mand_plot, &W_mand_connected_check);
    pack_widget_upper_right(&W_mand_connected_check, &W_mand_connected_depth_leftarrow);
    pack_widget_upper_right(&W_mand_connected_depth_leftarrow, &W_mand_connected_depth_label);
    pack_widget_upper_right(&W_mand_connected_depth_label, &W_mand_connected_depth_rightarrow);
    pack_widget_upper_right(&W_mand_plot, &W_mand_contains_half_check);
    pack_widget_upper_right(&W_mand_contains_half_check, &W_mand_contains_half_depth_leftarrow);
    pack_widget_upper_right(&W_mand_contains_half_depth_leftarrow, &W_mand_contains_half_depth_label);
    pack_widget_upper_right(&W_mand_contains_half_depth_label, &W_mand_contains_half_depth_rightarrow);
    pack_widget_upper_right(&W_mand_plot, &W_mand_trap_check);
    pack_widget_upper_right(&W_mand_trap_check, &W_mand_trap_depth_leftarrow);
    pack_widget_upper_right(&W_mand_trap_depth_leftarrow, &W_mand_trap_depth_label);
    pack_widget_upper_right(&W_mand_trap_depth_label, &W_mand_trap_depth_rightarrow);
    
    pack_widget_upper_right(&W_mand_plot, &W_limit_depth_title);
    pack_widget_upper_right(&W_mand_plot, &W_point_title);
    pack_widget_upper_right(&W_mand_plot, &W_point_connected_check);
    pack_widget_upper_right(&W_mand_plot, &W_point_contains_half_check);
    pack_widget_upper_right(&W_mand_plot, &W_point_trap_check);
    pack_widget_upper_right(&W_mand_plot, &W_point_uv_word_check);
  }
  
  //draw all the widgets
  for (int i=0; i<(int)widgets.size(); ++i) {
    widgets[i]->initial_draw();
  }
  
  
}

void IFSGui::main_loop() {
  XEvent e;
  XNextEvent(display, &e);
}




void IFSGui::launch(IFSWindowMode m) {
  
  //set the initial window mode
  window_mode = m;
  
  //set the initial settings for limit and mandlebrot
  limit_ll = cpx(-1, -1.5);
  limit_ur = cpx(2, 1.5);
  limit_depth = 12;
  limit_chunky = true;
  
  mand_ll = cpx(-1,-1);
  mand_ur = cpx(1,1);
  mand_pixel_group_size = 4;
  mand_connected = true;
  mand_connected_depth = 12;
  mand_contains_half = false;
  mand_contains_half_depth = 12;
  mand_trap = false;
  mand_trap_depth = 20;
  
  point_connected = true;
  point_contains_half = false;
  point_trap = false;
  point_uv_words = false;
  
  
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