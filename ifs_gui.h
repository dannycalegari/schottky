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
  Point2d<int> ul; //the upper left hand corner
  int height;
  int width;
  IFSGui* ifsg; //the gui, so we can call the member functions
  Pixmap p;
  GC gc;
  void (IFSGui::*click_signal)(XEvent*);
  
  virtual void initial_draw() {}
  virtual void redraw() {}
  virtual bool contains_pixel(int x, int y);
  Widget() {}
};

struct WidgetDraw : Widget {
  WidgetDraw() {}
  WidgetDraw(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*));
  void redraw();
  void initial_draw();
};

struct WidgetButton : Widget {
  std::string text;
  Point2d<int> text_position;
  
  WidgetButton() {}
  WidgetButton(IFSGui* i, const std::string& t, int w, int h, void (IFSGui::*f)(XEvent*));
  void initial_draw();
};

struct WidgetText : Widget {
  std::string text;
  Point2d<int> text_position;
  
  WidgetText() {}
  WidgetText(IFSGui* i, const std::string& t, int w, int h);
  void update_text(const std::string& s);
  void redraw();
  void initial_draw();
};

struct WidgetCheck : Widget {
  std::string text;
  bool checked;
  Point2d<int> text_position;
  
  WidgetCheck() {}
  WidgetCheck(IFSGui* i, const std::string& t, int w, int h, bool c, void (IFSGui::*f)(XEvent*));
  void redraw();
  void initial_draw();
};

struct WidgetLeftArrow : Widget {
  WidgetLeftArrow() {}
  WidgetLeftArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*));
  void initial_draw();
};

struct WidgetRightArrow : Widget {
  WidgetRightArrow() {}
  WidgetRightArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*));
  void initial_draw();
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
  bool limit_colors;
  double limit_pixel_width;
  Point2d<int> limit_cpx_to_pixel(const cpx& c);
  cpx limit_pixel_to_cpx(const Point2d<int>& p);
  
  //mandlebrot
  cpx mand_ll;
  cpx mand_ur;
  int mand_pixel_group_size;
  int mand_num_pixel_groups;
  double mand_pixel_group_width;
  double mand_pixel_width;
  bool mand_connected;
  int mand_connected_depth;
  bool mand_contains_half;
  int mand_contains_half_depth;
  bool mand_trap;
  int mand_trap_depth;
  cpx mand_pixel_group_to_cpx(const Point2d<int>& p); 
  cpx mand_pixel_to_cpx(const Point2d<int>& p);
  Point2d<int> mand_cpx_to_pixel(const cpx& c);
  int mand_get_color(const Point3d<int>& p);
  

  //data for mandlebrot
  std::vector<std::vector<Point3d<int> > > mand_data_grid;
  bool mand_grid_connected_valid;
  bool mand_grid_contains_half_valid;
  bool mand_grid_trap_valid;
  
  //data about highlighted point
  bool point_connected_check;
  int point_connected_depth;
  bool point_is_connected;
  bool point_contains_half_check;
  int point_contains_half_depth;
  bool point_is_contains_half;
  //bool point_trap_check;
  bool point_uv_words_check;
  int point_uv_words_depth;
  //bool point_is_trap;
  std::vector<std::pair<Bitword,Bitword> > point_uv_words;
  
  
  //computation functions
  void draw_limit();
  void recenter_limit(cpx c);
  void draw_mand();
  void change_highlighted_ifs(cpx c);
  void mand_zoom(double radius_multiplier);
  void mand_reset_mesh();
  void mand_recenter();
  void recompute_point_data();
  
  //graphics stuff
  Display* display;
  int screen;
  Window main_window;
  Colormap col_map;
  std::vector<Widget*> widgets;
  
  int get_rgb_color(double r, double g, double b);
  
  //widgets:
  WidgetButton W_switch_to_limit;
  WidgetButton W_switch_to_mandlebrot;
  WidgetButton W_switch_to_combined;
  
  WidgetText W_point_title;
  WidgetCheck W_point_connected_check;
  WidgetLeftArrow W_point_connected_leftarrow;
  WidgetText W_point_connected_depth_label;
  WidgetRightArrow W_point_connected_rightarrow;
  WidgetText W_point_connected_status;
  
  WidgetCheck W_point_contains_half_check;
  WidgetLeftArrow W_point_contains_half_leftarrow;
  WidgetText W_point_contains_half_depth_label;
  WidgetRightArrow W_point_contains_half_rightarrow;
  WidgetText W_point_contains_half_status;
  
  //WidgetCheck W_point_trap_check;
  //WidgetLeftArrow W_point_trap_decrease_depth;
  //WidgetText W_point_trap_depth_label;
  //WidgetRightArrow W_point_trap_increase_depth;
  //WidgetText W_point_trap_status;
  
  WidgetCheck W_point_uv_words_check;
  WidgetLeftArrow W_point_uv_words_leftarrow;
  WidgetText W_point_uv_words_depth_label;
  WidgetRightArrow W_point_uv_words_rightarrow;
  WidgetText W_point_uv_words_status;
  
  WidgetDraw W_limit_plot;
  WidgetText W_limit_depth_title;
  WidgetLeftArrow W_limit_depth_leftarrow;
  WidgetText W_limit_depth_label;
  WidgetRightArrow W_limit_depth_rightarrow;
  WidgetCheck W_limit_chunky;
  WidgetCheck W_limit_colors;
  WidgetText W_limit_zoom_title;
  WidgetButton W_limit_zoom_in;
  WidgetButton W_limit_zoom_out;
  WidgetText W_limit_center_title;
  
  
  WidgetDraw W_mand_plot;
  WidgetText W_mand_options_title;
  WidgetButton W_mand_recenter;
  WidgetText W_mand_zoom_title;
  WidgetButton W_mand_zoom_in;
  WidgetButton W_mand_zoom_out;
  WidgetText W_mand_mesh_title;
  WidgetLeftArrow  W_mand_mesh_leftarrow;
  WidgetText W_mand_mesh_label;
  WidgetRightArrow W_mand_mesh_rightarrow;
  
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
  
  //signal functions
  void S_switch_to_limit(XEvent* e);
  void S_switch_to_mandlebrot(XEvent* e);
  void S_switch_to_combined(XEvent* e);
  
  void S_limit_draw(XEvent* e);
  void S_limit_increase_depth(XEvent* e);
  void S_limit_decrease_depth(XEvent* e);
  void S_limit_switch_chunky(XEvent* e);
  void S_limit_switch_colors(XEvent* e);
  void S_limit_zoom_in(XEvent* e);
  void S_limit_zoom_out(XEvent* e);
  
  void S_mand_draw(XEvent* e);
  void S_mand_recenter(XEvent* e);
  void S_mand_zoom_in(XEvent* e);
  void S_mand_zoom_out(XEvent* e);
  void S_mand_decrease_mesh(XEvent* e);
  void S_mand_increase_mesh(XEvent* e);
  void S_mand_connected(XEvent* e);
  void S_mand_connected_increase_depth(XEvent* e);
  void S_mand_connected_decrease_depth(XEvent* e);
  void S_mand_contains_half(XEvent* e);
  void S_mand_contains_half_increase_depth(XEvent* e);
  void S_mand_contains_half_decrease_depth(XEvent* e);
  void S_mand_trap(XEvent* e);
  void S_mand_trap_increase_depth(XEvent* e);
  void S_mand_trap_decrease_depth(XEvent* e);
  
  void S_point_connected(XEvent* e);
  void S_point_connected_increase_depth(XEvent* e);
  void S_point_connected_decrease_depth(XEvent* e);
  void S_point_contains_half(XEvent* e);
  void S_point_contains_half_increase_depth(XEvent* e);
  void S_point_contains_half_decrease_depth(XEvent* e);
  //void S_point_trap(XEvent* e);
  //void S_point_connected_increase_depth(XEvent* e);
  //void S_point_connected_decrease_depth(XEvent* e);
  void S_point_uv_words(XEvent* e);
  void S_point_uv_words_increase_depth(XEvent* e);
  void S_point_uv_words_decrease_depth(XEvent* e);
  
  
  
  bool main_window_initialized;
  int main_window_height;
  int main_window_width;
  int limit_sidebar_size;
  int mand_sidebar_size;
  

  void pack_widget_upper_right(const Widget* w1, Widget* w2);
  void launch(IFSWindowMode m = BOTH, const cpx& c = cpx(0.5,0.5));
  void reset_and_pack_window();
  void main_loop();

};

#endif