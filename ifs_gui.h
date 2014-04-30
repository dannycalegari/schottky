#ifndef __IFS_GUI__
#define __IFS_GUI__

extern "C" {
#include <X11/Xlib.h>
#include <X11/Xutil.h>
}

#include <string>
#include <vector>

#include "ifs.h"
#include "point.h"

struct IFSGui;


struct Widget {
  Point2d<int> ul; //the upper left have corner
  int height;
  int width;
  IFSGui* ifsg; //the gui, so we can call the member functions
  Pixmap p;
  
  virtual void initial_draw();
  virtual void redraw();
  virtual bool contains_pixel();
  Widget() {
    //don't do anything
  }
};

struct WidgetDraw : Widget {
  WidgetDraw(IFSGui* i, int w, int h) {
    width = w;
    height = h;
    ifsg = i;
  }
};

struct WidgetButton : Widget {
  std::string text;
  void (IFSGui::*clicker)();
  WidgetButton(IFSGui* i, const std::string& t, int w, int h) {
    ifsg = i;
    width = w;
    height = h;
    text = t;
  }
};

struct WidgetText : Widget {
std::string text;
  WidgetText(IFSGui* i, const std::string& t, int w, int h) {
    ifsg = i;
    text = t;
    width = w; 
    height = h;
  }
};

struct WidgetCheck : Widget {
  std::string text;
  bool checked;
  void (IFSGui::*clicker)();
  
  WidgetCheck(IFSGui* i, const std::string& t, int w, int h, bool c) {
    width = w;
    height = h;
    ifsg = i;
    text = t;
    checked = c;
  } 
};

struct WidgetLeftArrow : Widget {
  void (IFSGui::*clicker)();
  WidgetLeftArrow(IFSGui* i, int w, int h) {
    ifsg = i;
    width = w;
    height = h;
  }
};

struct WidgetRightArrow : Widget {
  void (IFSGui::*clicker)();
  WidgetRightArrow(IFSGui* i, int w, int h) {
    ifsg = i;
    width = w;
    height = h;
  }
};

















enum IFSWindowMode { LIMIT, MANDLEBROT, BOTH };

struct IFSGui {
  
  //which mode
  IFSWindowMode window_mode;
  
  //limit and mandlebrot settings
  //limit
  ifs IFS;
  cpx limit_ll;
  cpx limit_ur;
  int limit_depth;
  bool limit_chunky;
  
  //mandlebrot
  cpx mand_ll;
  cpx mand_ur;
  int mand_pixel_group_size;
  bool mand_connected;
  int mand_connected_depth;
  bool mand_contains_half;
  int mand_contains_half_depth;
  bool mand_trap;
  int mand_trap_depth;
  bool point_connected;
  bool point_contains_half;
  bool point_trap;
  bool point_uv_words;
  
  //graphics stuff
  Display* display;
  int screen;
  Window main_window;
  
  std::vector<Widget*> widgets;
  
  //widgets:
  WidgetButton W_switch_to_limit;
  WidgetButton W_switch_to_mandlebrot;
  WidgetButton W_switch_to_combined;
  WidgetText W_point_title;
  WidgetCheck W_point_connected_check;
  WidgetCheck W_point_contains_half_check;
  WidgetCheck W_point_trap_check;
  WidgetCheck W_point_uv_word_check;
  
  WidgetDraw W_limit_plot;
  WidgetText W_limit_depth_title;
  WidgetLeftArrow W_limit_depth_leftarrow;
  WidgetText W_limit_depth_label;
  WidgetRightArrow W_limit_depth_rightarrow;
  WidgetText W_limit_chunky_title;
  WidgetCheck W_limit_chunky_on;
  WidgetCheck W_limit_chunky_off;
  WidgetText W_limit_zoom_title;
  WidgetButton W_limit_zoom_in;
  WidgetButton W_limit_zoom_out;
  WidgetText W_limit_center_title;
  
  WidgetDraw W_mand_plot;
  WidgetCheck W_mand_connected_check;
  WidgetLeftArrow W_mand_connected_depth_leftarrow;
  WidgetText W_mand_connected_depth_label;
  WidgetRightArrow W_mand_connected_depth_rightarrow;
  WidgetCheck W_mand_contains_half_check;
  WidgetLeftArrow W_mand_contains_half_depth_leftarrow;
  WidgetText W_mand_contains_half_depth_label;
  WidgetRightArrow W_mand_contains_half_depth_rightarrow;
  WidgetCheck W_mand_trap_check;
  WidgetLeftArrow W_mand_trap_depth_leftarrow;
  WidgetText W_mand_trap_depth_label;
  WidgetRightArrow W_mand_trap_depth_rightarrow;
  
  bool main_window_initialized;
  int main_window_height;
  int main_window_width;
  int limit_sidebar_size;
  int mand_sidebar_size;
  

  void pack_widget_upper_right(const Widget* w1, Widget* w2);
  void launch(IFSWindowMode m = BOTH);
  void reset_and_pack_window();
  void main_loop();

};

#endif