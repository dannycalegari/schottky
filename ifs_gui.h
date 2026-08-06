#ifndef __IFS_GUI__
#define __IFS_GUI__

extern "C" {
#include <X11/Xlib.h>
#include <X11/Xutil.h>
//named explicitly rather than relied on through Xutil.h, since we test for a
//dozen keysyms now and not just XK_q
#include <X11/keysym.h>
}

#include <string>
#include <vector>

#include "ifs.h"
#include "funddom_core.h"
#include "point.h"

struct IFSGui;

//Parses the text typed into W_point_param_entry -- "re im", "re+imi", a bare
//"deg" (degrees on the circle |s| = 1/sqrt(2), the certify_arc/funddom
//convention), or "deg@r" -- into a parameter value.  Returns false (leaving
//out untouched) if raw does not match any of those forms.  Declared here
//(rather than kept file-static in ifs_gui.cc) so schottky.cc can parse
//argv[1] with the exact same rules the GUI entry uses.
bool parse_parameter(const std::string& raw, cpx& out);

//Parses a pair of eventually-periodic words, "fg(fffffggggg) gf(gggggfffff)" or the same in
//bits, into the sign strings A and B of the landmark parameterisation; `why` gets a reason
//on failure.  See the comment on the definition in ifs_gui.cc.
bool parse_uv_pair(const std::string& raw, std::string& A, std::string& B, std::string& why);

//Parses the FINITE form, two equal-length words with no brackets ("fgff gfgg" or "0100 1011"),
//into the f/g strings the coincidence solver takes.  See the definition in ifs_gui.cc.
bool parse_uv_finite(const std::string& raw, std::string& u, std::string& v, std::string& why);

struct Widget {
  Point2d<int> ul; //the upper left hand corner
  int height;
  int width;
  IFSGui* ifsg; //the gui, so we can call the member functions
  Pixmap p;
  GC gc;
  void (IFSGui::*click_signal)(XEvent*);
  /* One line shown in the status bar while the pointer is over this widget.  Empty
     means "say nothing".  This is the cheapest way to make an opaque control
     discoverable: it costs the user no navigation and gets found by accident. */
  std::string help;
  
  virtual void initial_draw() = 0;
  virtual void redraw() = 0;
  //neither is pure: every widget except WidgetEntry is silent on keys, so
  //giving both a default here means only WidgetEntry has to say anything.
  virtual void handle_key(XEvent* e) {}
  virtual bool wants_keys() const { return false; }
  void clear();
  bool contains_pixel(int x, int y);
  bool intersects_rectangle(const Point2d<int>& p, int w, int h);
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
  void redraw();
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

//A single-line editable text field with a static label prefix (e.g. "s = ").
//Nothing else in the toolkit takes keyboard input; see the focus model on
//IFSGui (focus_widget) and the three main_loop edits that route KeyPress
//and ButtonPress there.
struct WidgetEntry : Widget {
  std::string label;          //static prefix, e.g. "s = "
  std::string text;           //editable content
  std::string committed;      //snapshot for Escape-revert
  int  caret;                 //0..text.size()
  int  max_chars;
  bool focused;
  Point2d<int> text_position;  //as WidgetText: (5, height/2 + (ascent-descent)/2)
  int  char_w, ascent, descent; //"fixed" is monospaced -> exact caret math
  int  label_px;
  void (IFSGui::*enter_signal)(XEvent*);   //fired on Return

  WidgetEntry() {}
  WidgetEntry(IFSGui* i, const std::string& lab, const std::string& init,
              int w, int h, void (IFSGui::*on_enter)(XEvent*));
  void initial_draw();
  void redraw();
  void handle_key(XEvent* e);
  bool wants_keys() const { return true; }
};

struct WidgetLeftArrow : Widget {
  WidgetLeftArrow() {}
  WidgetLeftArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*));
  void initial_draw();
  void redraw();
};

struct WidgetRightArrow : Widget {
  WidgetRightArrow() {}
  WidgetRightArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*));
  void initial_draw();
  void redraw();
};


struct IFSPath {
  bool is_valid;
  std::vector<cpx> path;
  bool closed;
  bool has_traps;
  std::vector<Ball> traps;
  std::vector<int> trap_colors;
  bool has_uv_words;
  std::vector<std::pair<Bitword,Bitword> > uv_words;
  bool has_half_words;
  std::vector<Bitword> half_words;
  int half_depth;
  int half_start;
  int half_end;
  int movie_length;
  int movie_fps;
  bool movie_with_mandelbrot;
  IFSPath() {
    is_valid = false;
    path.resize(0);
    closed = has_traps = has_uv_words = has_half_words = false;
    traps.resize(0);
    trap_colors.resize(0);
    uv_words.resize(0);
    movie_length = 10;
    movie_fps = 30;
    movie_with_mandelbrot = false;
    half_depth = 10;
    half_start = 0;
    half_end = 0;
  }
};



enum IFSWindowMode { LIMIT, MANDELBROT, BOTH };

struct IFSGui {
  
  //which mode
  IFSWindowMode window_mode;
  
  //limit and mandelbrot settings
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
  /* Draw the trap that certifies the current parameter, over the limit set.  A trap
     is what makes s a member of int(M), and until now the program could compute one
     and print its two words but never show it. */
  bool limit_trap;
  int limit_trap_depth;
  /* find_trap_words costs about 0.08s, and the limit-set overlay would otherwise redo
     it on every zoom and pan.  Cached on the exact parameter and depth it was computed
     for; the parameter is set exactly, so comparing cpx for equality is right here. */
  /* where the last drawn trap was, so the view can be taken to it: at a typical
     parameter the certifying balls are a few pixels across in the default view, and
     hunting for them by hand is most of the work of using the feature */
  bool limit_trap_located;
  cpx limit_trap_center;
  double limit_trap_radius;
  bool trap_cache_valid, trap_cache_ok;
  cpx trap_cache_z;
  int trap_cache_depth;
  std::vector<std::pair<Bitword,Bitword> > trap_cache_words;
  std::string trap_cache_why;
  bool limit_nifs;
  //the limit-set view to come back to when the difference-set layer is switched off: it
  //reframes the pane on 0, since Lambda - Lambda is centred there and Lambda is not
  cpx nifs_saved_ll, nifs_saved_ur;
  bool limit_gifs;
  bool limit_2d;
  std::vector<cpx> limit_marked_points;
  Point2d<int> limit_cpx_to_pixel(const cpx& c);
  cpx limit_pixel_to_cpx(const Point2d<int>& p);
  
  //mandelbrot
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
  /* Landmark points: the renormalization points, i.e. the sigma at which the
     limit-trap mechanism of CKW Lemma 9.2.5 applies.  Enumerated by
     fd_landmarks() in funddom_core.h; they depend only on the complexity bound
     mand_landmarks_N, not on the window, so the list is cached until N changes.
     They are drawn as an OVERLAY of dots rather than as a mand_data_grid layer,
     because they are a discrete set of points, not a property of every pixel. */
  /* Optional guide circles on the parameter pane.  |s| = 1/2 and |s| = 1/sqrt2
     bracket the annulus in which dM lies: below 1/2 the attractor is a Cantor set,
     at or above 1/sqrt2 it is robustly connected.  Drawing them on demand explains
     at a glance why several controls do nothing outside that band. */
  bool mand_circle_half;
  bool mand_circle_sqrt2;
  bool mand_landmarks;
  bool mand_roots;                 //show the finite-coincidence roots as well
  int  mand_roots_deg;             //largest polynomial degree enumerated
  int  mand_root_list_deg;         //what the cached list was built for (-1 = none)
  bool mand_root_list_on;
  int  mand_landmarks_N;            /* the complexity bound a+b <= N */
  std::vector<fd_landmark> mand_landmark_list;
  int  mand_landmark_list_N;        /* the N the cached list was built for; -1 = empty */
  int  mand_landmark_selected;      /* index into the list, or -1 */
  int  mand_annulus_res;            /* raster size of the fundamental-annulus figure */
  int  mand_annulus_cmax;           /* largest tail length c tried per parameter */
  int  mand_annulus_ball_depth;     /* depth for the trap-like ball computation */
  std::string mand_annulus_last;    /* summary of the last annulus run, for the readout */
  cpx mand_pixel_group_to_cpx(const Point2d<int>& p);
  Point2d<int> mand_cpx_to_pixel_group(const cpx& c);
  cpx mand_pixel_to_cpx(const Point2d<int>& p);
  Point2d<int> mand_cpx_to_pixel(const cpx& c);
  int mand_get_color(PointNd<4,int>& p);
  int mand_output_picture_size;
  

  //data for mandelbrot
  std::vector<std::vector<PointNd<4,int> > > mand_data_grid;
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
  bool point_uv_words_check;
  int point_uv_words_depth;
  std::vector<std::pair<Bitword,Bitword> > point_uv_words;
  bool point_trap_check;
  int point_trap_depth;
  std::vector<std::pair<Bitword,Bitword> > point_trap_words;

  //data about path
  IFSPath path;
  bool currently_drawing_path;
  void make_path_drawing_buttons();
  void make_path_task_buttons(bool created_by_drawing);
  void make_path_creation_buttons(bool cancelling);
  
  //computation functions
  void draw_limit();
  void draw_nifs_limit();
  void draw_gifs_limit();
  void draw_2d_limit();
  void recenter_limit(cpx c);
  void draw_mand();
  void change_highlighted_ifs(cpx c);
  void mand_zoom(double radius_multiplier);
  void mand_reset_mesh();
  void mand_recenter();
  bool window_frame_origin(int& ox, int& oy);
  void mand_draw_ball(const Ball& b, int col);
  //paints one mesh cell into draw_mand's client-side raster (NULL = straight to the pixmap)
  void mand_put_cell(XImage* img, int i, int j, unsigned long col);
  void recompute_point_data();
  void find_traps_along_path(int verbose);
  
  //graphics stuff
  Display* display;
  int screen;
  Window main_window;
  Colormap col_map;
  std::vector<Widget*> widgets;
  //the widget currently receiving KeyPress events, or NULL if none.  Reset to
  //NULL in launch, in reset_and_pack_window (at the widgets.resize(0) that
  //throws away everyone's ul), and in detach_widget when the widget being
  //detached is the one focused -- otherwise typing could go to a widget with
  //a stale ul or a pixmap that no longer belongs to the current mode.
  Widget* focus_widget;
  Widget* hover_widget;      /* last widget whose help we reported; avoids re-flushing */
  //the atom the window manager sends in a ClientMessage when the close button is
  //clicked; without answering it the server just kills us with an X error
  Atom wm_delete_window;

  int get_rgb_color(double r, double g, double b);

  //the one line at the very bottom of the window that says what the program is
  //doing.  set_status flushes, which is what makes a message appear during a
  //long computation instead of after it.
  void set_status(const std::string& s);
  void set_progress(double frac, const std::string& s);   /* status text plus a bar */
  bool status_widget_ready;

  //widgets:
  WidgetText W_status;
  WidgetButton W_switch_to_limit;
  WidgetButton W_switch_to_mandelbrot;
  WidgetButton W_switch_to_combined;

  WidgetText W_point_title;
  WidgetText W_point_point;
  //type "re im", "re+imi", a bare "deg" (degrees on |s|=1/sqrt2, the
  //certify_arc/funddom convention), or "deg@r", then press Enter to move
  //the highlighted parameter there
  WidgetEntry W_point_param_entry;
  WidgetEntry W_point_uv_entry;

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
  WidgetCheck W_limit_trap;
  WidgetLeftArrow W_limit_trap_depth_leftarrow;
  WidgetText W_limit_trap_depth_label;
  WidgetRightArrow W_limit_trap_depth_rightarrow;
  WidgetButton W_limit_trap_zoom;
  WidgetCheck W_limit_nifs;
  WidgetCheck W_limit_gifs;
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
  WidgetCheck W_mand_scale_check;
  WidgetCheck W_mand_circle_half_check;
  WidgetCheck W_mand_circle_sqrt2_check;
  /* Roots: every finite coincidence u(0) = v(0) of degree <= mand_roots_deg, marked and
     selectable exactly as the landmarks are.  Shares the marked-point list with them, so
     both can be shown at once and are told apart by colour and by their spec. */
  WidgetCheck      W_mand_roots_check;
  WidgetLeftArrow  W_mand_roots_leftarrow;
  WidgetText       W_mand_roots_label;
  WidgetRightArrow W_mand_roots_rightarrow;
  WidgetCheck W_mand_landmarks_check;
  WidgetLeftArrow W_mand_landmarks_leftarrow;
  WidgetText W_mand_landmarks_label;
  WidgetRightArrow W_mand_landmarks_rightarrow;
  WidgetButton W_mand_annulus_button;
  WidgetText W_mand_annulus_result;   /* the last annulus result, kept on screen */
  WidgetText W_mand_legend;           /* which layer wins where; see mand_get_color */
  WidgetText W_mand_mouse_label;
  WidgetText W_mand_mouse_X;
  WidgetText W_mand_mouse_Y;
  WidgetButton W_mand_output_window;
  WidgetButton W_mand_output_picture;
  WidgetLeftArrow W_mand_output_picture_size_leftarrow;
  WidgetText W_mand_output_picture_size_label;
  WidgetRightArrow W_mand_output_picture_size_rightarrow;
  
  WidgetText W_mand_path_drawing_title;
  WidgetButton W_mand_path_create_by_drawing_button;
  WidgetButton W_mand_path_create_by_boundary_button;
  WidgetButton W_mand_path_finish_cancel_button;
  WidgetButton W_mand_path_finish_path_button;
  WidgetButton W_mand_path_finish_loop_button;
  WidgetText W_mand_path_tasks_title;
  WidgetButton W_mand_path_delete_button;
  WidgetButton W_mand_path_find_traps_button;
  WidgetButton W_mand_path_create_movie_button;
  WidgetText W_mand_path_movie_length_title;
  WidgetLeftArrow W_mand_path_movie_decrease_length;
  WidgetText W_mand_path_movie_length_label;
  WidgetRightArrow W_mand_path_movie_increase_length;
  WidgetCheck W_mand_path_movie_with_mandelbrot;
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

  //the click_signal shared by every WidgetEntry, since a click on a text
  //entry means "put the caret here", not any entry-specific action
  void S_entry_click(XEvent* e);

  void S_switch_to_limit(XEvent* e);
  void S_switch_to_mandelbrot(XEvent* e);
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
  void S_limit_trap(XEvent* e);
  void S_limit_trap_increase_depth(XEvent* e);
  void S_limit_trap_decrease_depth(XEvent* e);
  void S_limit_trap_zoom(XEvent* e);
  /* The single trap search, shared by the point panel and the limit-set overlay, so
     that the words the one reports are the balls the other draws. */
  bool find_trap_words(int uv_depth, std::vector<std::pair<Bitword,Bitword> >& out,
                       std::string& why);
  void draw_limit_trap();
  void S_limit_uv_graph_decrease_depth(XEvent* e);
  void S_limit_uv_graph_increase_depth(XEvent* e);
  void S_limit_nifs(XEvent* e);
  void S_limit_gifs(XEvent* e);
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
  void S_mand_circle_half(XEvent* e);
  void S_mand_circle_sqrt2(XEvent* e);
  void mand_draw_guide_circles();
  /* Scale bar and axes.  A deep zoom in parameter space is otherwise disorienting: the
     picture is self-similar, so nothing in it says how far in you are. */
  void mand_draw_scale();
  bool mand_scale;
  void S_mand_scale(XEvent* e);
  void mand_update_legend();
  void S_mand_landmarks(XEvent* e);
  void S_mand_landmarks_decrease(XEvent* e);
  void S_mand_landmarks_increase(XEvent* e);
  void S_mand_annulus(XEvent* e);      /* the limit-trap figure at the selected landmark */
  void mand_rebuild_landmarks();       /* refresh the cache if N changed */
  void mand_draw_landmarks();          /* overlay the dots on W_mand_plot */
  int  mand_landmark_near(int wx, int wy, int tol_px);  /* index, or -1 */
  void S_mand_output_window(XEvent* e);
  void S_mand_output_picture(XEvent* e);
  void S_mand_output_picture_increase_size(XEvent* e);
  void S_mand_output_picture_decrease_size(XEvent* e);
  

  //enter_signal for W_point_param_entry: parses its text and, on success,
  //moves the highlighted parameter there
  void S_point_param_entered(XEvent* e);
  void S_point_uv_entered(XEvent* e);
  void S_point_roots(XEvent* e);
  void S_mand_roots(XEvent* e);
  void S_mand_roots_increase(XEvent* e);
  void S_mand_roots_decrease(XEvent* e);
  void do_uv_roots_search();
  void uv_show_results(const std::vector<fd_landmark>& found,
                       const std::string& lead, const char* what);

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
  
  void S_mand_path_create_by_drawing_button(XEvent* e);
  void S_mand_path_create_by_boundary(XEvent* e);
  void S_mand_path_finish_cancel(XEvent* e);
  void S_mand_path_finish_path(XEvent* e);
  void S_mand_path_finish_loop(XEvent* e);
  void S_mand_path_delete(XEvent* e);
  void S_mand_path_find_traps(XEvent* e);
  void S_mand_path_create_movie(XEvent* e);
  void S_mand_path_movie_decrease_length(XEvent* e);
  void S_mand_path_movie_increase_length(XEvent* e);
  void S_mand_path_movie_with_mandelbrot(XEvent* e);
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
  //the window the cached landmark list was built for, and whether it came from the
  //targeted search (which depends on the window) or the exhaustive one (which does not)
  bool mand_landmark_list_targeted;
  //the list came from a typed u,v pair rather than from enumeration or search, so it does
  //not depend on the window and must not be rebuilt when the window changes
  bool mand_landmark_list_from_uv;
  cpx  mand_landmark_list_ll, mand_landmark_list_ur;
  //where the WM frame was last seen, so a mode switch (which destroys and recreates the
  //window) puts it back where the user left it instead of at the top left corner
  int main_window_x;
  int main_window_y;
  /* One font for the whole interface, loaded on demand and kept.  See gui_font(). */
  XFontStruct* the_gui_font;
  XFontStruct* gui_font();
  int limit_sidebar_size;
  int mand_sidebar_size;
  
  

  void detach_widget(Widget* w);
  void pack_widget_upper_right(const Widget* w1, Widget* w2);
  //Tab from an entry: walks widgets (in packing order) for the next one
  //with wants_keys() == true, wrapping around; a no-op if none exists.
  void focus_next_entry();
  void launch(IFSWindowMode m = BOTH, const cpx& c = cpx(0.5,0.5));
  void reset_and_pack_window();
  void main_loop();

};

#endif