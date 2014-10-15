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
  virtual void clear();
  bool contains_pixel(int x, int y);
  bool intersects_rectangle(const Point2d<int>& ul, int w, int h);
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


struct IFSPath {
  bool is_valid;
  std::vector<cpx> path;
  bool closed;
  bool has_traps;
  std::vector<Ball> traps;
  std::vector<int> trap_colors;
  bool has_coordinates;
  std::vector<std::pair<double, double> > coordinates;
  bool has_uv_words;
  std::vector<std::pair<Bitword,Bitword> > uv_words;
  bool has_half_words;
  std::vector<Bitword> half_words;
  int half_depth;
  int half_start;
  int half_end;
  int movie_length;
  int movie_fps;
  bool movie_with_mandlebrot;
  IFSPath() {
    is_valid = false;
    path.resize(0);
    closed = has_traps = has_coordinates = has_uv_words = has_half_words = false;
    traps.resize(0);
    trap_colors.resize(0);
    coordinates.resize(0);
    uv_words.resize(0);
    movie_length = 10;
    movie_fps = 30;
    movie_with_mandlebrot = false;
    half_depth = 10;
    half_start = 0;
    half_end = 0;
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
  bool limit_auto_depth;
  bool limit_chunky;
  bool limit_colors;
  double limit_pixel_width;
  bool limit_uv_graph;
  int limit_uv_graph_depth;
  bool limit_nifs;
  bool limit_2d;
  std::vector<cpx> limit_marked_points;
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
  bool mand_limit_trap;
  bool mand_dirichlet;
  int mand_dirichlet_depth;
  bool mand_set_C;
  int mand_set_C_depth;
  bool mand_theta;
  int mand_theta_depth;
  cpx mand_pixel_group_to_cpx(const Point2d<int>& p);
  Point2d<int> mand_cpx_to_pixel_group(const cpx& c);
  cpx mand_pixel_to_cpx(const Point2d<int>& p);
  Point2d<int> mand_cpx_to_pixel(const cpx& c);
  int mand_get_color(PointNd<6,int>& p);
  

  //data for mandlebrot
  std::vector<std::vector<PointNd<6,int> > > mand_data_grid;
  bool mand_grid_connected_valid;
  bool mand_grid_contains_half_valid;
  bool mand_grid_trap_valid;
  bool mand_grid_dirichlet_valid;
  bool mand_grid_set_C_valid;
  bool mand_grid_theta_valid;
  
  //data about highlighted point
  bool point_connected_check;
  int point_connected_depth;
  bool point_is_connected;
  bool point_contains_half_check;
  int point_contains_half_depth;
  bool point_is_contains_half;
  bool point_uv_words_check;
  int point_uv_words_depth;
  std::vector<std::pair<Bitword,Bitword> > point_uv_words;
  bool point_trap_check;
  int point_trap_depth;
  std::vector<std::pair<Bitword,Bitword> > point_trap_words;
  bool point_coordinates_check;
  int point_coordinates_depth;
  double point_coordinates_theta;
  double point_coordinates_lambda;
  
  //data about path
  IFSPath path;
  bool currently_drawing_path;
  void make_path_drawing_buttons();
  void make_path_task_buttons(bool created_by_drawing);
  void make_path_creation_buttons(bool cancelling);
  
  //computation functions
  void draw_limit();
  void draw_nifs_limit();
  void draw_2d_limit();
  void recenter_limit(cpx c);
  void draw_mand();
  void change_highlighted_ifs(cpx c);
  void mand_zoom(double radius_multiplier);
  void mand_reset_mesh();
  void mand_recenter();
  void mand_draw_ball(const Ball& b, int col);
  void recompute_point_data();
  void find_traps_along_path(int verbose);
  void find_coordinates_along_path(int verbose);
  
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
  WidgetText W_point_point;
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
  
  WidgetCheck W_point_trap_check;
  WidgetLeftArrow W_point_trap_leftarrow;
  WidgetText W_point_trap_depth_label;
  WidgetRightArrow W_point_trap_rightarrow;
  WidgetText W_point_trap_status;
  
  WidgetCheck W_point_uv_words_check;
  WidgetLeftArrow W_point_uv_words_leftarrow;
  WidgetText W_point_uv_words_depth_label;
  WidgetRightArrow W_point_uv_words_rightarrow;
  WidgetText W_point_uv_words_status;
  
  WidgetCheck W_point_coordinates_check;
  WidgetLeftArrow W_point_coordinates_leftarrow;
  WidgetText W_point_coordinates_depth_label;
  WidgetRightArrow W_point_coordinates_rightarrow;
  WidgetText W_point_coordinates_status;
  
  WidgetDraw W_limit_plot;
  WidgetText W_limit_depth_title;
  WidgetLeftArrow W_limit_depth_leftarrow;
  WidgetText W_limit_depth_label;
  WidgetRightArrow W_limit_depth_rightarrow;
  WidgetCheck W_limit_depth_auto;
  WidgetCheck W_limit_chunky;
  WidgetCheck W_limit_colors;
  WidgetText W_limit_zoom_title;
  WidgetButton W_limit_zoom_in;
  WidgetButton W_limit_zoom_out;
  WidgetText W_limit_center_title;
  WidgetCheck W_limit_uv_graph;
  WidgetText W_limit_uv_graph_depth_title;
  WidgetLeftArrow W_limit_uv_graph_depth_leftarrow;
  WidgetText W_limit_uv_graph_depth_label;
  WidgetRightArrow W_limit_uv_graph_depth_rightarrow;
  WidgetCheck W_limit_nifs;
  WidgetCheck W_limit_2d;
  
  
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
  WidgetCheck W_mand_limit_trap_check;
  WidgetCheck W_mand_dirichlet_check;
  WidgetLeftArrow W_mand_dirichlet_depth_leftarrow;
  WidgetText W_mand_dirichlet_depth_label;
  WidgetRightArrow W_mand_dirichlet_depth_rightarrow;
  WidgetCheck W_mand_set_C_check;
  WidgetLeftArrow W_mand_set_C_depth_leftarrow;
  WidgetText W_mand_set_C_depth_label;
  WidgetRightArrow W_mand_set_C_depth_rightarrow;
  WidgetCheck W_mand_theta_check;
  WidgetLeftArrow W_mand_theta_depth_leftarrow;
  WidgetText W_mand_theta_depth_label;
  WidgetRightArrow W_mand_theta_depth_rightarrow;
  WidgetText W_mand_mouse_label;
  WidgetText W_mand_mouse_X;
  WidgetText W_mand_mouse_Y;
  WidgetButton W_mand_output_window;
  
  WidgetText W_mand_path_drawing_title;
  WidgetButton W_mand_path_create_by_drawing_button;
  WidgetButton W_mand_path_create_by_boundary_button;
  WidgetButton W_mand_path_finish_cancel_button;
  WidgetButton W_mand_path_finish_path_button;
  WidgetButton W_mand_path_finish_loop_button;
  WidgetText W_mand_path_tasks_title;
  WidgetButton W_mand_path_delete_button;
  WidgetButton W_mand_path_find_traps_button;
  WidgetButton W_mand_path_find_coordinates_button;
  WidgetButton W_mand_path_create_movie_button;
  WidgetText W_mand_path_movie_length_title;
  WidgetLeftArrow W_mand_path_movie_decrease_length;
  WidgetText W_mand_path_movie_length_label;
  WidgetRightArrow W_mand_path_movie_increase_length;
  WidgetCheck W_mand_path_movie_with_mandlebrot;
  WidgetButton W_mand_path_find_uv_words_button;
  WidgetButton W_mand_path_find_half_words_button;
  WidgetText W_mand_path_half_depth_title;
  WidgetLeftArrow W_mand_path_half_depth_leftarrow;
  WidgetText W_mand_path_half_depth_label;
  WidgetRightArrow W_mand_path_half_depth_rightarrow;
  WidgetText W_mand_path_half_start_stop_title;
  WidgetLeftArrow W_mand_path_half_start_leftarrow;
  WidgetText W_mand_path_half_start_label;
  WidgetRightArrow W_mand_path_half_start_rightarrow;
  WidgetLeftArrow W_mand_path_half_end_leftarrow;
  WidgetText W_mand_path_half_end_label;
  WidgetRightArrow W_mand_path_half_end_rightarrow;
  
  
  //signal functions
  void S_switch_to_limit(XEvent* e);
  void S_switch_to_mandlebrot(XEvent* e);
  void S_switch_to_combined(XEvent* e);
  
  void S_limit_draw(XEvent* e);
  void S_limit_increase_depth(XEvent* e);
  void S_limit_decrease_depth(XEvent* e);
  void S_limit_auto_depth(XEvent* e);
  void S_limit_switch_chunky(XEvent* e);
  void S_limit_switch_colors(XEvent* e);
  void S_limit_zoom_in(XEvent* e);
  void S_limit_zoom_out(XEvent* e);
  void S_limit_uv_graph(XEvent* e);
  void S_limit_uv_graph_decrease_depth(XEvent* e);
  void S_limit_uv_graph_increase_depth(XEvent* e);
  void S_limit_nifs(XEvent* e);
  void S_limit_2d(XEvent* e);
  
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
  void S_mand_limit_trap(XEvent* e);
  void S_mand_dirichlet(XEvent* e);
  void S_mand_dirichlet_decrease_depth(XEvent* e);
  void S_mand_dirichlet_increase_depth(XEvent* e);
  void S_mand_set_C(XEvent* e);
  void S_mand_set_C_decrease_depth(XEvent* e);
  void S_mand_set_C_increase_depth(XEvent* e);
  void S_mand_theta(XEvent* e);
  void S_mand_theta_decrease_depth(XEvent* e);
  void S_mand_theta_increase_depth(XEvent* e);
  void S_mand_output_window(XEvent* e);
  
  
  void S_point_connected(XEvent* e);
  void S_point_connected_increase_depth(XEvent* e);
  void S_point_connected_decrease_depth(XEvent* e);
  void S_point_contains_half(XEvent* e);
  void S_point_contains_half_increase_depth(XEvent* e);
  void S_point_contains_half_decrease_depth(XEvent* e);
  void S_point_trap(XEvent* e);
  void S_point_trap_increase_depth(XEvent* e);
  void S_point_trap_decrease_depth(XEvent* e);
  void S_point_uv_words(XEvent* e);
  void S_point_uv_words_increase_depth(XEvent* e);
  void S_point_uv_words_decrease_depth(XEvent* e);
  void S_point_coordinates(XEvent* e);
  void S_point_coordinates_increase_depth(XEvent* e);
  void S_point_coordinates_decrease_depth(XEvent* e);
  
  void S_mand_path_create_by_drawing_button(XEvent* e);
  void S_mand_path_create_by_boundary(XEvent* e);
  void S_mand_path_finish_cancel(XEvent* e);
  void S_mand_path_finish_path(XEvent* e);
  void S_mand_path_finish_loop(XEvent* e);
  void S_mand_path_delete(XEvent* e);
  void S_mand_path_find_traps(XEvent* e);
  void S_mand_path_find_coordinates(XEvent* e);
  void S_mand_path_create_movie(XEvent* e);
  void S_mand_path_movie_decrease_length(XEvent* e);
  void S_mand_path_movie_increase_length(XEvent* e);
  void S_mand_path_movie_with_mandlebrot(XEvent* e);
  void S_mand_path_find_uv_words(XEvent* e);
  void S_mand_path_find_half_words(XEvent* e);
  void S_mand_path_half_increase_depth(XEvent* e);
  void S_mand_path_half_decrease_depth(XEvent* e);
  void S_mand_path_half_increase_start(XEvent* e);
  void S_mand_path_half_decrease_start(XEvent* e);
  void S_mand_path_half_increase_end(XEvent* e);
  void S_mand_path_half_decrease_end(XEvent* e);
  
  bool main_window_initialized;
  int main_window_height;
  int main_window_width;
  int limit_sidebar_size;
  int mand_sidebar_size;
  
  

  void detach_widget(Widget* w);
  void pack_widget_upper_right(const Widget* w1, Widget* w2);
  void launch(IFSWindowMode m = BOTH, const cpx& c = cpx(0.5,0.5));
  void reset_and_pack_window();
  void main_loop();

};

#endif