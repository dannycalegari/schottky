#include "ifs_gui.h"


#include "movie.h"







/***************************************************************************
 * widget functions
 ***************************************************************************/
bool Widget::contains_pixel(int x, int y) {
  return (ul.x <= x) && (x < ul.x + width) && (ul.y <= y) && (y < ul.y + height);
}

void Widget::clear() {
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, ifsg->main_window, gc, ul.x, ul.y, width, height);
}


WidgetDraw::WidgetDraw(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*)) {
  width = w;
  height = h;
  ifsg = i;
  click_signal = f;
  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XDrawRectangle(ifsg->display, p, gc, 0, 0, width-1, height-1);
}

void WidgetDraw::redraw() {
  initial_draw();
}

void WidgetDraw::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}



WidgetButton::WidgetButton(IFSGui* i, const std::string& t, int w, int h, void (IFSGui::*f)(XEvent*)) {
  ifsg = i;
  width = w;
  height = h;
  text = t;
  click_signal = f;
  
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  XFontStruct* font = XLoadQueryFont(ifsg->display, "fixed");
  XSetFont(ifsg->display, gc, font->fid);  
  XCharStruct te;
  int fdir, fdescent, fascent;
  XTextExtents(font, text.c_str(), text.size(), &fdir, &fascent, &fdescent, &te);
  
  int desired_width = te.rbearing - te.lbearing;
  int width_offset = desired_width/2;
  int height_offset = (te.ascent - te.descent)/2; 
  if (w > 0) {
    width = w; 
  } else {
    width = desired_width + 10;
  }
  text_position = Point2d<int>(width/2 - width_offset, height/2 + height_offset); 

  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  
  //clear the pixmap
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  
  //draw the square
  XSetLineAttributes(ifsg->display, gc, 0.1, LineSolid, CapButt, JoinMiter);
  XDrawLine(ifsg->display, p, gc, 1, 1, 1, height-2);
  XDrawLine(ifsg->display, p, gc, 1, height-2, width-2, height-2);
  XDrawLine(ifsg->display, p, gc, width-2, height-2, width-2, 1);
  XDrawLine(ifsg->display, p, gc, width-2, 1, 1, 1);
  //draw the label
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size()); 
}

void WidgetButton::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}


WidgetText::WidgetText(IFSGui* i, const std::string& t, int w, int h) {
  ifsg = i;
  text = t;
  height = h;
  click_signal = NULL;
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  XFontStruct* font = XLoadQueryFont(ifsg->display, "fixed");
  XSetFont(ifsg->display, gc, font->fid);  
  XCharStruct te;
  int fdir, fdescent, fascent;
  XTextExtents(font, text.c_str(), text.size(), &fdir, &fascent, &fdescent, &te);
  
  int desired_width = te.rbearing - te.lbearing;
  //int width_offset = desired_width/2;
  int height_offset = (te.ascent - te.descent)/2; 
  if (w > 0) {
    width = w; 
  } else {
    width = desired_width + 10;
  }
  text_position = Point2d<int>(5, height/2 + height_offset); 
  //text_position = Point2d<int>(5, height/2);
  
  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  
  //clear the pixmap
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  //draw the label
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size()); 
}

void WidgetText::update_text(const std::string& s) {
  text = s;
  //clear it
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size()); 
  redraw();
} 

void WidgetText::redraw() {
  initial_draw();
}

void WidgetText::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
  //std::cout << "Drawing string: " << text << "\n";
  //std::cout << "At position: " << ul.x << " " << ul.y << "\n";
}

WidgetCheck::WidgetCheck(IFSGui* i, const std::string& t, int w, int h, bool c, void (IFSGui::*f)(XEvent*)) {
  width = w;
  height = h;
  ifsg = i;
  text = t;
  checked = c;
  click_signal = f;
  
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  XFontStruct* font = XLoadQueryFont(ifsg->display, "fixed");
  XSetFont(ifsg->display, gc, font->fid);  
  XCharStruct te;
  int fdir, fdescent, fascent;
  XTextExtents(font, text.c_str(), text.size(), &fdir, &fascent, &fdescent, &te);
  
  int desired_width = te.rbearing - te.lbearing;
  //int width_offset = desired_width/2;
  int height_offset = (te.ascent - te.descent)/2; 
  if (w > 0) {
    width = w; 
  } else {
    width = desired_width + 10 + 10 + 5;
  }
  text_position = Point2d<int>(20, height/2 + height_offset+1); 
  
  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  
  //clear the pixmap
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  
  //draw the check box
  if (checked) {
    XFillRectangle(ifsg->display, p, gc, 5, height/2-5, 11, 11);
  } else {
    XDrawRectangle(ifsg->display, p, gc, 5, height/2-5, 10, 10);
  }
  
  //draw the label
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size()); 
} 

void WidgetCheck::redraw() {
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  if (checked) {
    XFillRectangle(ifsg->display, p, gc, 5, height/2-5, 11, 11);
  } else {
    XDrawRectangle(ifsg->display, p, gc, 5, height/2-5, 10, 10);
  }
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size());  
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}
  

void WidgetCheck::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}


WidgetLeftArrow::WidgetLeftArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*)) {
  ifsg = i;
  width = w;
  height = h;
  click_signal = f;
  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  //clear the pixmap
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  
  //draw the left arrow box
  XPoint tri[3];
  tri[0].x = width-5; tri[0].y = 5;
  tri[1].x = 5; tri[1].y = height/2;
  tri[2].x = width-5; tri[2].y = height-5;
  XFillPolygon(ifsg->display, p, gc, tri, 3, Convex, CoordModeOrigin);
}

void WidgetLeftArrow::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}



WidgetRightArrow::WidgetRightArrow(IFSGui* i, int w, int h, void (IFSGui::*f)(XEvent*)) {
  ifsg = i;
  width = w;
  height = h;
  click_signal = f;
  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  //clear the pixmap
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XSetBackground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  
  //draw the right arrow box
  XPoint tri[3];
  tri[0].x = 5; tri[0].y = 5;
  tri[1].x = 5; tri[1].y = height-5;
  tri[2].x = width-5; tri[2].y = height/2;
  XFillPolygon(ifsg->display, p, gc, tri, 3, Convex, CoordModeOrigin);
}

void WidgetRightArrow::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}
















/****************************************************************************
 * signal handlers
 ****************************************************************************/
 
void IFSGui::S_switch_to_limit(XEvent* e) {
  if (e->type != ButtonPress) return;
  window_mode = LIMIT;
  reset_and_pack_window();
}

void IFSGui::S_switch_to_mandlebrot(XEvent* e) {
  if (e->type != ButtonPress) return;
  window_mode = MANDLEBROT;
  reset_and_pack_window();
}

void IFSGui::S_switch_to_combined(XEvent* e) {
  if (e->type != ButtonPress) return;
  window_mode = BOTH;
  reset_and_pack_window();
}
 
//limit set
void IFSGui::S_limit_draw(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  
  //if we clicked, center the window
  //on the mouse pointer
  if (e->type == ButtonPress) {
    int widget_x = e->xbutton.x - W_limit_plot.ul.x;
    int widget_y = e->xbutton.y - W_limit_plot.ul.y;
    cpx c = limit_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    recenter_limit(c);
  }
}

void IFSGui::S_limit_increase_depth(XEvent* e) {
  if (e->type == MotionNotify) return;
  ++limit_depth;
  std::stringstream T; T.str(""); T << limit_depth;
  W_limit_depth_label.update_text(T.str());
  draw_limit();
}

void IFSGui::S_limit_decrease_depth(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  --limit_depth;
  std::stringstream T; T.str(""); T << limit_depth;
  W_limit_depth_label.update_text(T.str());
  draw_limit();
}


void IFSGui::S_limit_auto_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_auto_depth = !limit_auto_depth;
  W_limit_depth_auto.checked = limit_auto_depth;
  W_limit_depth_auto.redraw();
  draw_limit();
}


void IFSGui::S_limit_switch_chunky(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  limit_chunky = !limit_chunky;
  W_limit_chunky.checked = limit_chunky;
  W_limit_chunky.redraw();
  draw_limit();
}

void IFSGui::S_limit_switch_colors(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  limit_colors = !limit_colors;
  W_limit_colors.checked = limit_colors;
  W_limit_colors.redraw();
  draw_limit();
}

void IFSGui::S_limit_zoom_in(XEvent* e) {
  if (e->type == ButtonPress) {
    double radius = (limit_ur.real() - limit_ll.real())/2.0;
    radius /= 1.5;
    cpx center = (limit_ll + limit_ur) / 2.0;
    limit_ll = center - cpx(radius, radius);
    limit_ur = center + cpx(radius, radius);
    limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(W_limit_plot.width);
    draw_limit();
  }
}

void IFSGui::S_limit_zoom_out(XEvent* e) {
  if (e->type == ButtonPress) {
    double radius = (limit_ur.real() - limit_ll.real())/2.0;
    radius *= 1.5;
    cpx center = (limit_ll + limit_ur) / 2.0;
    limit_ll = center - cpx(radius, radius);
    limit_ur = center + cpx(radius, radius);
    limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(W_limit_plot.width);
    draw_limit();
  }
}

void IFSGui::S_limit_uv_graph(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_uv_graph = !limit_uv_graph;
  W_limit_uv_graph.checked = limit_uv_graph;
  W_limit_uv_graph.redraw();
  draw_limit();
}

void IFSGui::S_limit_uv_graph_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  if (limit_uv_graph_depth == 1) return;
  --limit_uv_graph_depth;
  std::stringstream T; T.str(""); T << limit_uv_graph_depth;
  W_limit_uv_graph_depth_label.update_text(T.str());
  draw_limit();
}

void IFSGui::S_limit_uv_graph_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  ++limit_uv_graph_depth;
  std::stringstream T; T.str(""); T << limit_uv_graph_depth;
  W_limit_uv_graph_depth_label.update_text(T.str());
  draw_limit();
}



//mandlebrot
void IFSGui::S_mand_draw(XEvent* e) {
  if (e->type == KeyPress) return;
        
  //if we're drawing a loop, and we click, we need to add that to 
  //the loop
  if (e->type == ButtonPress && currently_drawing_path) {
    int widget_x = e->xbutton.x - W_mand_plot.ul.x;
    int widget_y = e->xbutton.y - W_mand_plot.ul.y;
    cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    path.path.push_back(c);
    //draw the next line
    int rcol = get_rgb_color(1,0.2,1);
    Widget& MW = W_mand_plot;
    XSetForeground(display, MW.gc, rcol);
    if (path.path.size() == 1) {
      XDrawPoint(display, MW.p, MW.gc, widget_x, widget_y);
      XCopyArea(display, MW.p, main_window, MW.gc, widget_x, widget_y, 
                                                   1,1,
                                                   MW.ul.x + widget_x, MW.ul.y + widget_y); 
    } else {
      Point2d<int> p = mand_cpx_to_pixel(path.path[path.path.size()-2]);
      XDrawLine(display, MW.p, MW.gc, p.x, p.y, widget_x, widget_y);
      XCopyArea(display, MW.p, main_window, MW.gc, 0,0, 
                                                   MW.width, MW.height,
                                                   MW.ul.x, MW.ul.y); 
    }
  
  //the following is run if the button is pressed or if there is
  //motion where the button is down
  } else if ( (e->type == ButtonPress && e->xbutton.button == Button1) ||
              (e->type == MotionNotify && ((e->xmotion.state >> 8)&1)) ) {
    int widget_x = e->xbutton.x - W_mand_plot.ul.x;
    int widget_y = e->xbutton.y - W_mand_plot.ul.y;
    cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    //std::cout << "Moving IFS to " << widget_x << " " << widget_y << "=" << c << "\n";
    change_highlighted_ifs(c);
  
  //it's a right mouse click -- zoom in
  } else if (e->type == ButtonPress && e->xbutton.button == Button3) { 
    int widget_x = e->xbutton.x - W_mand_plot.ul.x;
    int widget_y = e->xbutton.y - W_mand_plot.ul.y;
    cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    IFS.set_params(c,c);
    mand_zoom(0.5);
    recompute_point_data();
  }
  
  //additionally, if the mouse is moved, we need to update the 
  //text
  if (e->type == MotionNotify) {
    int widget_x = e->xbutton.x - W_mand_plot.ul.x;
    int widget_y = e->xbutton.y - W_mand_plot.ul.y;
    cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    std::stringstream T; T.str("");
    T << "Mouse: " << c;
    W_mand_mouse_label.update_text(T.str());
  }
  
}


void IFSGui::S_mand_recenter(XEvent* e) {
  if (e->type == ButtonPress) mand_recenter();
}

void IFSGui::S_mand_zoom_in(XEvent* e) {
  if (e->type == ButtonPress) mand_zoom(0.5);
}

void IFSGui::S_mand_zoom_out(XEvent* e) {
  if (e->type == ButtonPress) mand_zoom(2.0);
}

void IFSGui::S_mand_decrease_mesh(XEvent* e) {
  if (e->type == ButtonPress) {
    if (mand_pixel_group_size == 1) return;
    mand_pixel_group_size /= 2;
    std::stringstream T; T.str(""); T << mand_pixel_group_size;
    W_mand_mesh_label.update_text(T.str());
    mand_reset_mesh();
    draw_mand();
  }
}

void IFSGui::S_mand_increase_mesh(XEvent* e) {
  if (e->type == ButtonPress) {
    mand_pixel_group_size *= 2;
    std::stringstream T; T.str(""); T << mand_pixel_group_size;
    W_mand_mesh_label.update_text(T.str());
    mand_reset_mesh();
    draw_mand();
  }
}

void IFSGui::S_mand_connected(XEvent* e) {
  if (e->type == ButtonPress) {
    mand_connected = !mand_connected;
    W_mand_connected_check.checked = mand_connected;
    W_mand_connected_check.redraw();
    draw_mand();
  }
}

void IFSGui::S_mand_connected_increase_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    ++mand_connected_depth;
    std::stringstream T; T.str(""); T << mand_connected_depth;
    W_mand_connected_depth_label.update_text(T.str());
    mand_grid_connected_valid = false;
    draw_mand();
  }
}

void IFSGui::S_mand_connected_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    --mand_connected_depth;
    std::stringstream T; T.str(""); T << mand_connected_depth;
    W_mand_connected_depth_label.update_text(T.str());
    mand_grid_connected_valid = false;
    draw_mand();
  }
}

void IFSGui::S_mand_contains_half(XEvent* e) {
  if (e->type == ButtonPress) {
    mand_contains_half = !mand_contains_half;
    W_mand_contains_half_check.checked = mand_contains_half;
    W_mand_contains_half_check.redraw();
    draw_mand();
  }
}

void IFSGui::S_mand_contains_half_increase_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    ++mand_contains_half_depth;
    std::stringstream T; T.str(""); T << mand_contains_half_depth;
    W_mand_contains_half_depth_label.update_text(T.str());
    mand_grid_contains_half_valid = false;
    draw_mand();
  }
}

void IFSGui::S_mand_contains_half_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    --mand_contains_half_depth;
    std::stringstream T; T.str(""); T << mand_contains_half_depth;
    W_mand_contains_half_depth_label.update_text(T.str());
    mand_grid_contains_half_valid = false;
    draw_mand();
  }
}

void IFSGui::S_mand_trap(XEvent* e) {
  if (e->type == ButtonPress) {
    mand_trap = !mand_trap;
    W_mand_trap_check.checked = mand_trap;
    W_mand_trap_check.redraw();
    draw_mand();
  }
}
void IFSGui::S_mand_trap_increase_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    ++mand_trap_depth;
    std::stringstream T; T.str(""); T << mand_trap_depth;
    W_mand_trap_depth_label.update_text(T.str());
    mand_grid_trap_valid = false;
    draw_mand();
  }
}
void IFSGui::S_mand_trap_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    --mand_trap_depth;
    std::stringstream T; T.str(""); T << mand_trap_depth;
    W_mand_trap_depth_label.update_text(T.str());
    mand_grid_trap_valid = false;
    draw_mand();
  }
}

//point
void IFSGui::S_point_connected(XEvent* e) {
  if (e->type != ButtonPress) return;
  point_connected_check = !point_connected_check;
  W_point_connected_check.checked = point_connected_check;
  W_point_connected_check.redraw();
  recompute_point_data();
}

void IFSGui::S_point_connected_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  ++point_connected_depth;
  std::stringstream T; T.str(""); T << point_connected_depth;
  W_point_connected_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_connected_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  --point_connected_depth;
  std::stringstream T; T.str(""); T << point_connected_depth;
  W_point_connected_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_contains_half(XEvent* e) {
  if (e->type != ButtonPress) return;
  point_contains_half_check = !point_contains_half_check;
  W_point_contains_half_check.checked = point_contains_half_check;
  W_point_contains_half_check.redraw();
  recompute_point_data();
}

void IFSGui::S_point_contains_half_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  ++point_contains_half_depth;
  std::stringstream T; T.str(""); T << point_contains_half_depth;
  W_point_contains_half_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_contains_half_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  --point_contains_half_depth;
  std::stringstream T; T.str(""); T << point_contains_half_depth;
  W_point_contains_half_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_uv_words(XEvent* e) {
  if (e->type != ButtonPress) return;
  point_uv_words_check = !point_uv_words_check;
  W_point_uv_words_check.checked = point_uv_words_check;
  W_point_uv_words_check.redraw();
  recompute_point_data();
}
void IFSGui::S_point_uv_words_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  ++point_uv_words_depth;
  std::stringstream T; T.str(""); T << point_uv_words_depth;
  W_point_uv_words_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_uv_words_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  --point_uv_words_depth;
  std::stringstream T; T.str(""); T << point_uv_words_depth;
  W_point_uv_words_depth_label.update_text(T.str());
  recompute_point_data();
}





void IFSGui::S_mand_path_create_by_drawing_button(XEvent* e) {
  if (e->type != ButtonPress) return;
  currently_drawing_path = true;
  path = IFSPath();
  path.is_valid = true;
  make_path_drawing_buttons();
}
  
void IFSGui::S_mand_path_create_by_boundary(XEvent* e) {
  if (e->type != ButtonPress) return;
  //make a grid of what's connected
  std::vector<std::vector<bool> > grid(mand_num_pixel_groups, 
                                       std::vector<bool>(mand_num_pixel_groups, false));
  for (int i=0; i<mand_num_pixel_groups; ++i) {
    for (int j=0; j<mand_num_pixel_groups; ++j) {
      grid[i][j] = (mand_data_grid[i][mand_num_pixel_groups-j-1].x > 0);
    }
  }
  Point2d<int> p = mand_cpx_to_pixel_group(IFS.z);
  p.y = mand_num_pixel_groups - p.y -1;
  path = IFSPath();
  IFS.hole_boundary_containing_point_from_grid(path.path, path.closed, grid, p, mand_ll, mand_ur, 0); 
  path.is_valid = true;
  make_path_task_buttons(false);
  draw_mand();
}

void IFSGui::S_mand_path_finish_cancel(XEvent* e) {
  if (e->type != ButtonPress) return;
  currently_drawing_path = false;
  path = IFSPath();
  make_path_creation_buttons(true);
  draw_mand();
}

void IFSGui::S_mand_path_finish_path(XEvent* e) {
  if (e->type != ButtonPress) return;
  currently_drawing_path = false;
  path.closed = false;
  make_path_task_buttons(true);
}

void IFSGui::S_mand_path_finish_loop(XEvent* e) {
  if (e->type != ButtonPress) return;
  currently_drawing_path = false;
  path.closed = true;
  make_path_task_buttons(true);
  Point2d<int> p1 = mand_cpx_to_pixel(path.path[path.path.size()-1]);
  Point2d<int> p2 = mand_cpx_to_pixel(path.path[0]);
  Widget& MW = W_mand_plot;
  XSetForeground(display, MW.gc, get_rgb_color(1,0.2,1));
  XDrawLine(display, MW.p, MW.gc, p1.x, p1.y, p2.x, p2.y);
  XCopyArea(display, MW.p, main_window, MW.gc, 0,0, 
                                               MW.width, MW.height,
                                               MW.ul.x, MW.ul.y); 
}

void IFSGui::S_mand_path_delete(XEvent* e) {
  if (e->type != ButtonPress) return;
  path = IFSPath();
  make_path_creation_buttons(false);
  draw_mand();
}

void IFSGui::S_mand_path_find_traps(XEvent* e) {
  if (e->type != ButtonPress) return;
  find_traps_along_path(0);
}

void IFSGui::S_mand_path_create_movie(XEvent* e) {
  if (e->type != ButtonPress || !path.is_valid) return;
  //create the mandlebrot connectedness grid
  std::vector<std::vector<bool> > mand_connected_grid;
  if (path.movie_with_mandlebrot) {
    mand_connected_grid = std::vector<std::vector<bool> >(mand_data_grid.size(), std::vector<bool>(mand_data_grid.size(), false));
    for (int i=0; i<(int)mand_data_grid.size(); ++i) {
      for (int j=0; j<(int)mand_data_grid.size(); ++j) {
        mand_connected_grid[i][mand_data_grid.size()-j-1] = (mand_data_grid[i][j].x >= 0);
      }
    }
  }
  (void)ifs_movie_from_path(IFS, path.path, path.closed, "ifs_movie",
                            limit_ll, limit_ur, limit_depth, 
                            path.movie_with_mandlebrot, &mand_ll, &mand_ur, &mand_connected_grid, 
                            W_limit_plot.width, W_limit_plot.height, 
                            path.movie_fps, path.movie_length, 1);
}

void IFSGui::S_mand_path_movie_decrease_length(XEvent* e) {
  if (e->type != ButtonPress) return;
  if (path.movie_length >= 1) {
    --path.movie_length;
    std::stringstream T; T.str(""); T << path.movie_length;
    W_mand_path_movie_length_label.update_text(T.str());
  }
}
void IFSGui::S_mand_path_movie_increase_length(XEvent* e) {
  if (e->type != ButtonPress) return;
  ++path.movie_length;
  std::stringstream T; T.str(""); T << path.movie_length;
  W_mand_path_movie_length_label.update_text(T.str());
}

void IFSGui::S_mand_path_movie_with_mandlebrot(XEvent* e) {
  if (e->type != ButtonPress) return;
  path.movie_with_mandlebrot = !path.movie_with_mandlebrot;
  W_mand_path_movie_with_mandlebrot.checked = path.movie_with_mandlebrot;
  W_mand_path_movie_with_mandlebrot.redraw();
}

void IFSGui::S_mand_path_find_uv_words(XEvent* e) {
  if (e->type != ButtonPress || !path.is_valid) return;
  ifs temp_IFS;
  temp_IFS.find_closest_uv_words_along_path(path.path, path.closed, point_uv_words_depth); 
}

void IFSGui::make_path_drawing_buttons() {
  detach_widget(&W_mand_path_create_by_drawing_button);
  detach_widget(&W_mand_path_create_by_boundary_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_drawing_title);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_cancel_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_path_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_loop_button);
  W_mand_path_drawing_title.initial_draw();
  W_mand_path_finish_cancel_button.initial_draw();
  W_mand_path_finish_path_button.initial_draw();
  W_mand_path_finish_loop_button.initial_draw();
}

void IFSGui::make_path_task_buttons(bool created_by_drawing) {
  //if it was created by drawing, we need to take those away
  if (created_by_drawing) {
    detach_widget(&W_mand_path_drawing_title);
    detach_widget(&W_mand_path_finish_cancel_button);
    detach_widget(&W_mand_path_finish_path_button);
    detach_widget(&W_mand_path_finish_loop_button);
  } else {
    detach_widget(&W_mand_path_create_by_drawing_button);
    detach_widget(&W_mand_path_create_by_boundary_button);
  } 
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_tasks_title);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_delete_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_traps_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_movie_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_length_title);
  pack_widget_upper_right(&W_mand_path_movie_length_title, &W_mand_path_movie_decrease_length);
  pack_widget_upper_right(&W_mand_path_movie_decrease_length, &W_mand_path_movie_length_label);
  pack_widget_upper_right(&W_mand_path_movie_length_label, &W_mand_path_movie_increase_length);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_with_mandlebrot);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_uv_words_button);
  W_mand_path_tasks_title.initial_draw();
  W_mand_path_delete_button.initial_draw();
  W_mand_path_find_traps_button.initial_draw();
  W_mand_path_create_movie_button.initial_draw();
  W_mand_path_movie_length_title.initial_draw();
  W_mand_path_movie_decrease_length.initial_draw();
  std::stringstream T; T.str(""); T << path.movie_length;
  W_mand_path_movie_length_label.update_text(T.str());
  W_mand_path_movie_with_mandlebrot.checked = path.movie_with_mandlebrot;
  W_mand_path_movie_with_mandlebrot.initial_draw();
  W_mand_path_movie_with_mandlebrot.redraw();
  W_mand_path_movie_increase_length.initial_draw();
  W_mand_path_find_uv_words_button.initial_draw();
}

void IFSGui::make_path_creation_buttons(bool cancelling) {
  if (cancelling) {
    detach_widget(&W_mand_path_drawing_title);
    detach_widget(&W_mand_path_finish_cancel_button);
    detach_widget(&W_mand_path_finish_path_button);
    detach_widget(&W_mand_path_finish_loop_button);
  } else {
    detach_widget(&W_mand_path_tasks_title);
    detach_widget(&W_mand_path_delete_button);
    detach_widget(&W_mand_path_find_traps_button);
    detach_widget(&W_mand_path_create_movie_button);
    detach_widget(&W_mand_path_movie_length_title);
    detach_widget(&W_mand_path_movie_decrease_length);
    detach_widget(&W_mand_path_movie_length_label);
    detach_widget(&W_mand_path_movie_increase_length);
    detach_widget(&W_mand_path_movie_with_mandlebrot);
    detach_widget(&W_mand_path_find_uv_words_button);
  }
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_by_drawing_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_by_boundary_button);
  W_mand_path_create_by_drawing_button.initial_draw();
  W_mand_path_create_by_boundary_button.initial_draw();
}


/**************************************************************************
 * back end functions
 * ************************************************************************/
Point2d<int> IFSGui::limit_cpx_to_pixel(const cpx& c) {
  int x = int( (c.real() - limit_ll.real()) / limit_pixel_width );
  int real_y = int( (c.imag() - limit_ll.imag()) / limit_pixel_width );
  int x11_y = W_limit_plot.height - real_y;
  return Point2d<int>(x, x11_y);
}

cpx IFSGui::limit_pixel_to_cpx(const Point2d<int>& p) {
  double r = (p.x+0.5)*limit_pixel_width + limit_ll.real();
  double i = limit_ur.imag() - (p.y+0.5)*limit_pixel_width;
  return cpx(r,i);
}

void IFSGui::draw_limit() {
  double min_r;
  if (!IFS.minimal_enclosing_radius(min_r)) {
    return;
  }
  
  Widget& LW = W_limit_plot;
  
  //clear the limit widget
  XSetForeground(display, LW.gc, WhitePixel(display, screen));
  XFillRectangle(display, LW.p, LW.gc, 0, 0, LW.width, LW.height);
  XSetForeground(display, LW.gc, BlackPixel(display, screen));
  XDrawRectangle(display, LW.p, LW.gc, 0, 0, LW.width-1, LW.height-1);
  XSetFillStyle(display, LW.gc, FillSolid);
  
  int blue_color = (limit_uv_graph ? get_rgb_color(0.5,0.9,1) : get_rgb_color(0,0.6,1.0));
  int yellow_color = (limit_uv_graph ? get_rgb_color(1,0.8,0.5) : get_rgb_color(1.0,0.6,0));
  
  Ball initial_ball(0.5,(IFS.z-1.0)/2.0,(1.0-IFS.w)/2.0,min_r);
  std::vector<std::pair<bool,Ball> > stack(0);
  stack.push_back(std::make_pair(false, initial_ball));
  while (stack.size() > 0) {
    std::pair<bool,Ball> b = stack.back();
    stack.pop_back();
    //if the ball is disjoint from the window, we might as well 
    //get rid of it
    if (!b.first && b.second.is_disjoint(limit_ll, limit_ur)) continue;
    if ( (!limit_auto_depth && b.second.word_len >= limit_depth) ||
         (limit_auto_depth && b.second.radius < limit_pixel_width/2.0) ) {
      Point2d<int> p = limit_cpx_to_pixel(b.second.center);
      double r = b.second.radius / limit_pixel_width;
      if (r <= 1.0) r = 1.0;
      if (limit_colors) {
        XSetForeground(display, LW.gc, (b.second.last_gen_index()==0 ? blue_color : yellow_color));
      }
      if (limit_chunky) {
        XFillArc(display, LW.p, LW.gc, p.x-r, p.y-r, int(2.0*r), int(2.0*r), 23040, 23040);
      } else {
        XDrawPoint(display, LW.p, LW.gc, p.x, p.y);
      }
      continue;
    }
    //if the ball isn't disjoint from the window, maybe it is contained in it?
    Ball bz = IFS.act_on_right(0, b.second);
    Ball bw = IFS.act_on_right(1, b.second);
    if (b.first) {
      stack.push_back(std::make_pair(true, bz));
      stack.push_back(std::make_pair(true, bw));
    } else {
      if (b.second.is_contained(limit_ll, limit_ur)) {
        stack.push_back(std::make_pair(true, bz));
        stack.push_back(std::make_pair(true, bw));
      } else {
        stack.push_back(std::make_pair(false, bz));
        stack.push_back(std::make_pair(false, bw));
      }
    }
  }
  
  //draw the marked points 0, 1/2, 1
  for (int i=0; i<(int)limit_marked_points.size(); ++i) {
    cpx& c = limit_marked_points[i];
    int rcol = get_rgb_color(1,0.1,0);
    if (limit_ll.real() < c.real() && c.real() < limit_ur.real() &&
        limit_ll.imag() < c.imag() && c.imag() < limit_ur.imag()) {
      Point2d<int> p = limit_cpx_to_pixel(c);
      XSetForeground(display, LW.gc, rcol);
      XFillArc(display, LW.p, LW.gc, p.x-3, p.y-3, 3, 3, 23040, 23040);
    }
  }
  
  //draw the uv graph
  if (limit_uv_graph) {
    //generate all the balls of the given depth
    std::stringstream T; T.str("");
    std::vector<Ball> uv_graph_balls;
    std::vector<Point3d<int> > uv_graph_edges; //(i,j,k) records an edge between i and j which swaps a suffix of length k
    IFS.compute_uv_graph(uv_graph_edges, uv_graph_balls, limit_uv_graph_depth, 0);
  
    //get the size of a text string of this length
    XFontStruct* font = XLoadQueryFont(display, "fixed");
    XCharStruct te;
    int fdir, fdescent, fascent;
    T << Bitword(uv_graph_balls[0].word, uv_graph_balls[0].word_len);
    XTextExtents(font, T.str().c_str(), T.str().size(), &fdir, &fascent, &fdescent, &te);
    int text_width_offset = (te.rbearing - te.lbearing)/2;
    int text_height_offset = (te.ascent - te.descent)/2; 
  
    //draw the balls
    XSetForeground(display, LW.gc, BlackPixel(display, screen));
    double r = 15;
    for (int i=0; i<(int)uv_graph_balls.size(); ++i) {
      Point2d<int> p = limit_cpx_to_pixel(uv_graph_balls[i].center);
      //double r = uv_graph_balls[i].radius / limit_pixel_width;
      T.str(""); T << Bitword(uv_graph_balls[i].word, uv_graph_balls[i].word_len);
      XDrawArc(display, LW.p, LW.gc, p.x-r, p.y-r, int(2*r), int(2*r), 23040, 23040);
      XDrawString(display, LW.p, LW.gc, p.x-text_width_offset, p.y+text_height_offset, T.str().c_str(), T.str().size()); 
    }
    
    //draw the edges
    T.str("8");
    XTextExtents(font, T.str().c_str(), T.str().size(), &fdir, &fascent, &fdescent, &te);
    text_width_offset = (te.rbearing - te.lbearing)/2;
    text_height_offset = (te.ascent - te.descent)/2; 
    XSetLineAttributes(display, LW.gc, 1.5, LineSolid, 1, 1);
    for (int i=0; i<(int)uv_graph_edges.size(); ++i){
      Point3d<int>& e = uv_graph_edges[i];
      Point2d<int> c1 = limit_cpx_to_pixel(uv_graph_balls[e.x].center);
      Point2d<int> c2 = limit_cpx_to_pixel(uv_graph_balls[e.y].center);
      Point2d<float> v(c2.x-c1.x, c2.y-c1.y);
      v = v/(float)sqrt(dot(v,v));
      Point2d<float> c1p(c1.x + r*v.x, c1.y + r*v.y);
      Point2d<float> c2p(c2.x - r*v.x, c2.y - r*v.y);
      Point2d<float> perp(-10*v.y, 10*v.x);
      Point2d<float> text_center = Point2d<float>((c1p.x + c2p.x)/2.0, (c1p.y + c2p.y)/2.03);
      text_center = text_center + perp;
      T.str(""); T << e.z;
      XDrawString(display, LW.p, LW.gc, text_center.x-text_width_offset, text_center.y+text_height_offset, T.str().c_str(), T.str().size());
      XDrawLine(display, LW.p, LW.gc, c1p.x, c1p.y, c2p.x, c2p.y);
    }
  }
  
  LW.redraw();
}

void IFSGui::recenter_limit(cpx c) {
  double radius = (limit_ur.real() - limit_ll.real())/2.0;
  limit_ll = c - cpx(radius, radius);
  limit_ur = c + cpx(radius, radius);
  draw_limit();
}




cpx IFSGui::mand_pixel_group_to_cpx(const Point2d<int>& p) {
  double r = (p.x + 0.5)*mand_pixel_group_width + mand_ll.real();
  double i = mand_ur.imag() - (p.y + 0.5)*mand_pixel_group_width;
  return cpx(r,i);
}

Point2d<int> IFSGui::mand_cpx_to_pixel_group(const cpx& c) {
  return Point2d<int>( (c.real() - mand_ll.real()) / mand_pixel_group_width,
                       mand_num_pixel_groups - ((c.imag() - mand_ll.imag())/ mand_pixel_group_width) );
}

cpx IFSGui::mand_pixel_to_cpx(const Point2d<int>& p) {
  double r = (p.x + 0.5)*mand_pixel_width + mand_ll.real();
  double i = mand_ur.imag() - (p.y + 0.5)*mand_pixel_width;
  return cpx(r,i);
}

Point2d<int> IFSGui::mand_cpx_to_pixel(const cpx& c) {
  return Point2d<int>( (c.real() - mand_ll.real()) / mand_pixel_width,
                       W_mand_plot.height - ((c.imag() - mand_ll.imag()) / mand_pixel_width) );
}

int IFSGui::mand_get_color(const Point3d<int>& p) {
  if (mand_trap && p.z > 0) { //use the trap color
    return get_rgb_color(0, double(p.z)/100, 1.0);
  } else if (mand_contains_half && p.y > 0) {
    return get_rgb_color(0.5, double(p.y)/100, 0.5);
  } else if (mand_connected && p.x >= 0) {
    return p.x*0x000001;
  } else {
    return WhitePixel(display, screen);
  }
}

void IFSGui::mand_draw_ball(const Ball& b, int col) {
  Point2d<int> p = mand_cpx_to_pixel(b.center);
  Widget& MP = W_mand_plot;
  XSetForeground(display, MP.gc, col);
  int r = int(b.radius / mand_pixel_width);
  if (r < 2) r = 2;
  XFillArc(display, MP.p, MP.gc, p.x-r, p.y-r, 2*r, 2*r, 23040, 23040);
  XCopyArea(display, MP.p, main_window, MP.gc, p.x-r, p.y-r, 2*r,2*r, MP.ul.x+p.x-r, MP.ul.y+p.y-r);
}
  
  
  
  

void IFSGui::draw_mand() {
  ifs temp_IFS;
  Widget& MW = W_mand_plot;

  std::vector<Ball> TLB;
  bool found_TLB = false;
  double TLB_neighborhood;
  if (mand_trap && !mand_grid_trap_valid) {
    //std::cout << "About to find TLB\n";
    temp_IFS.set_params(IFS.z, IFS.z);
    temp_IFS.TLB_for_region(TLB, TLB_neighborhood, mand_ll, mand_ur, 15, 0);
    found_TLB = (TLB.size() != 0);
  }
  
  for (int i=0; i<(int)mand_num_pixel_groups; ++i) {
    for (int j=0; j<(int)mand_num_pixel_groups; ++j) {
      
      //do the necessary computations
      cpx c = mand_pixel_group_to_cpx(Point2d<int>(i,j));
      temp_IFS.set_params(c,c);
      if (mand_connected && !mand_grid_connected_valid) {
        if (!temp_IFS.is_connected(mand_connected_depth, mand_data_grid[i][j].x) ) {
          mand_data_grid[i][j].x = -1;
        }
      }
      if (mand_contains_half && !mand_grid_contains_half_valid) {
        if (!temp_IFS.contains_half(mand_contains_half_depth, mand_data_grid[i][j].y)) {
          mand_data_grid[i][j].y = -1;
        }
      }
      if (mand_trap && !mand_grid_trap_valid && found_TLB) {
        double trap_radius;
        int multiplier = 100/mand_trap_depth;
        mand_data_grid[i][j].z = multiplier*temp_IFS.check_TLB(TLB,trap_radius,TLB_neighborhood,mand_trap_depth);
      }
      
      //draw the pixel for the impatient
      int col = mand_get_color(mand_data_grid[i][j]);
      XSetForeground(display, MW.gc, col);
      XFillRectangle(display, MW.p, MW.gc, i*mand_pixel_group_size, 
                                           j*mand_pixel_group_size, 
                                           mand_pixel_group_size, 
                                           mand_pixel_group_size);
      XCopyArea(display, MW.p, main_window, MW.gc, i*mand_pixel_group_size, 
                                                   j*mand_pixel_group_size, 
                                                   mand_pixel_group_size, 
                                                   mand_pixel_group_size, 
                                                   MW.ul.x + i*mand_pixel_group_size,
                                                   MW.ul.y + j*mand_pixel_group_size);
    }
  }
  if (mand_connected && !mand_grid_connected_valid) mand_grid_connected_valid = true;
  if (mand_contains_half && !mand_grid_contains_half_valid) mand_grid_contains_half_valid = true;
  if (mand_trap && !mand_grid_trap_valid) mand_grid_trap_valid = true;
  
  //now draw the highlighted point
  Point2d<int> h = mand_cpx_to_pixel(IFS.z);
  int rcol = get_rgb_color(1.0,0.1,0.0);
  XSetForeground(display, MW.gc, rcol);
  XFillArc(display, MW.p, MW.gc, h.x-2, h.y-2, 4, 4, 23040, 23040);
  
  //now draw the path, if there is one
  if (path.is_valid) {
    int rcol = get_rgb_color(1,0.2,1);
    XSetForeground(display, MW.gc, rcol);
    for (int i=0; i<(int)path.path.size()-1; ++i) {
      Point2d<int> p1 = mand_cpx_to_pixel(path.path[i]);
      Point2d<int> p2 = mand_cpx_to_pixel(path.path[i+1]);
      XDrawLine(display, MW.p, MW.gc, p1.x, p1.y, p2.x, p2.y);
    }
    if (path.closed) {
      Point2d<int> p1 = mand_cpx_to_pixel(path.path[path.path.size()-1]);
      Point2d<int> p2 = mand_cpx_to_pixel(path.path[0]);
      XDrawLine(display, MW.p, MW.gc, p1.x, p1.y, p2.x, p2.y);
    }
    //if the path has traps, draw them
    if (path.has_traps) {
      for (int i=0; i<(int)path.traps.size(); ++i) {
        XSetForeground(display, MW.gc, path.trap_colors[i]);
        Point2d<int> p = mand_cpx_to_pixel(path.traps[i].center);
        int r = int(path.traps[i].radius / mand_pixel_width);
        if (r < 2) r = 2;
        XFillArc(display, MW.p, MW.gc, p.x-r, p.y-r, 2*r, 2*r, 23040, 23040);
      }
    }
  }
  
  MW.redraw();
}


  
void IFSGui::change_highlighted_ifs(cpx c) {
  //redraw the current red dot with the saved stuff
  Widget& MW = W_mand_plot;
  Point2d<int> p = mand_cpx_to_pixel(IFS.z);
  int pg_i = (p.x/mand_pixel_group_size) - 2;
  if (pg_i < 0) pg_i = 0;
  int upper_limit_i = pg_i + 6;
  if (upper_limit_i >= mand_num_pixel_groups) upper_limit_i = mand_num_pixel_groups-1;
  int pg_j = (p.y/mand_pixel_group_size) - 2;
  if (pg_j < 0) pg_j = 0;
  int upper_limit_j = pg_j + 6;
  if (upper_limit_j >= mand_num_pixel_groups) upper_limit_j = mand_num_pixel_groups-1;
  
  for (int i=pg_i; i<upper_limit_i; ++i) {
    for (int j=pg_j; j<upper_limit_j; ++j) {
      int col = mand_get_color(mand_data_grid[i][j]);
      XSetForeground(display, MW.gc, col);
      XFillRectangle(display, MW.p, MW.gc, i*mand_pixel_group_size, 
                                           j*mand_pixel_group_size, 
                                           mand_pixel_group_size, 
                                           mand_pixel_group_size);
    }
  }
  XCopyArea(display, MW.p, main_window, MW.gc, pg_i*mand_pixel_group_size, 
                                               pg_j*mand_pixel_group_size, 
                                               6*mand_pixel_group_size, 
                                               6*mand_pixel_group_size, 
                                               MW.ul.x + pg_i*mand_pixel_group_size,
                                               MW.ul.y + pg_j*mand_pixel_group_size);
  //now actually switch the point
  IFS.set_params(c,c);
  Point2d<int> h = mand_cpx_to_pixel(c);
  int rcol = get_rgb_color(1.0,0.1,0.0);
  XSetForeground(display, MW.gc, rcol);
  XFillArc(display, MW.p, MW.gc, h.x-2, h.y-2, 5, 5, 23040, 23040);
  XCopyArea(display, MW.p, main_window, MW.gc,h.x-2, h.y-2, 6, 6, 
                                              MW.ul.x + h.x-2, MW.ul.y + h.y-2);
  
  if (window_mode != MANDLEBROT) {
    draw_limit();
  }
  recompute_point_data();
  
}


//this zooms on whatever the highlighted ifs is
void IFSGui::mand_zoom(double radius_multiplier) {
  double radius = (mand_ur.real() - mand_ll.real())/2.0;
  cpx c = IFS.z;
  radius *= radius_multiplier;
  mand_ll = c - cpx(radius, radius);
  mand_ur = c + cpx(radius, radius);
  mand_reset_mesh();
  draw_mand();
}

//recenters on the highlighted ifs
void IFSGui::mand_recenter() {
  double radius = (mand_ur.real() - mand_ll.real())/2.0;
  cpx c = IFS.z;
  mand_ll = c - cpx(radius, radius);
  mand_ur = c + cpx(radius, radius);
  mand_grid_connected_valid = false;
  mand_grid_contains_half_valid = false;
  mand_grid_trap_valid = false;
  draw_mand();
}

void IFSGui::mand_reset_mesh() {
  mand_pixel_width = (mand_ur.real() - mand_ll.real()) / double(W_mand_plot.width);
  mand_pixel_group_width = mand_pixel_group_size * mand_pixel_width;
  mand_num_pixel_groups = W_mand_plot.width / mand_pixel_group_size;
  mand_data_grid.resize(mand_num_pixel_groups);
  for (int i=0; i<mand_num_pixel_groups; ++i) {
    mand_data_grid[i] = std::vector<Point3d<int> >(mand_num_pixel_groups, Point3d<int>(-1,-1,-1));
  }
  mand_grid_connected_valid = false;
  mand_grid_contains_half_valid = false;
  mand_grid_trap_valid = false;
}




void IFSGui::recompute_point_data() {
  std::stringstream T;
  int difficulty;
  T.str(""); T << "Location: " << IFS.z;
  W_point_point.update_text(T.str());
  
  if (!point_connected_check) {
    T.str(""); T << "(disabled)";
  } else {
    T.str("");
    T << ( (point_is_connected = IFS.is_connected(point_connected_depth, difficulty)) ? "yes" : "no" );
  }
  W_point_connected_status.update_text(T.str());
  
  if (!point_contains_half_check) {
     T.str(""); T << "(disabled)";
  } else {
    T.str("");
    T << ( (point_is_contains_half = IFS.contains_half(point_contains_half_depth, difficulty)) ? "yes" : "no" );
  }
  W_point_contains_half_status.update_text(T.str());
  
  if (!point_uv_words_check) {
    T.str(""); T << "(disabled)";
  } else {
    IFS.find_closest_uv_words(point_uv_words, point_uv_words_depth);
    T.str("");
    T << point_uv_words[0].first << " " << point_uv_words[0].second;
  }
  W_point_uv_words_status.update_text(T.str());
  
  
}



void IFSGui::find_traps_along_path(int verbose) {
  if (!path.is_valid || currently_drawing_path || path.path.size() == 0) return;
  //find the extents of the path, plus a little extra
  cpx box_ll = path.path[0];
  cpx box_ur = path.path[0];
  for (int i=0; i<(int)path.path.size(); ++i) {
    if (path.path[i].real() < box_ll.real()) {
      box_ll = cpx(path.path[i].real(), box_ll.imag());
    }
    if (path.path[i].real() > box_ur.real()) {
      box_ur = cpx(path.path[i].real(), box_ur.imag());
    }
    if (path.path[i].imag() < box_ll.imag()) {
      box_ll = cpx(box_ll.real(), path.path[i].imag());
    }
    if (path.path[i].imag() > box_ur.real()) {
      box_ur = cpx(box_ur.real(), path.path[i].imag());
    }
  }
  double bw = box_ur.real() - box_ll.real();
  double bh = box_ur.imag() - box_ll.imag();
  box_ur = box_ur + cpx(0.1*bw, 0.1*bh);
  box_ll = box_ll - cpx(0.1*bw, 0.1*bh);
  cpx av = (box_ll + box_ur)/2.0;
  
  //find the TLB for this region
  std::vector<Ball> TLB;
  double TLB_neighborhood;
  ifs temp_IFS;
  temp_IFS.set_params(av, av);
  if (!temp_IFS.TLB_for_region(TLB, TLB_neighborhood, box_ll, box_ur, 15, verbose)) {
    std::cout << "Couldn't find TLB for box " << box_ll << " " << box_ur << " at depth " << mand_trap_depth << "\n";
    return;
  }
  
  //find traps along all the path segments
  path.traps.resize(0);
  path.trap_colors.resize(0);
  path.has_traps = true;
  int upper_index_bound = (path.closed ? path.path.size()-1 : path.path.size()-2);
  for (int i=0; i<=upper_index_bound; ++i) {
    int ip1 = (i == (int)path.path.size()-1 ? 0 : i+1);
    cpx current_z = path.path[i];
    cpx end_z = path.path[ip1];
    if (verbose>0) {
      std::cout << "Finding traps along segment " << current_z << " - " << end_z << "\n";
    }
    do {
      double epsilon = -1;
      int difficulty = -1;
      temp_IFS.set_params(current_z, current_z);
      if ( (difficulty = temp_IFS.check_TLB(TLB, epsilon, TLB_neighborhood, mand_trap_depth)) < 0 ) {
        std::cout << "Failed to find trap at " << current_z << "\n";
        return;
      }
      if (verbose>0) {
        std::cout << "Found trap " << current_z << " - " << epsilon << "\n";
      }
      path.traps.push_back(Ball(current_z, epsilon));
      double gamount = double(difficulty)/double(mand_trap_depth);
      path.trap_colors.push_back( get_rgb_color(0.5,gamount,1) );
      cpx v = end_z-current_z;
      current_z = current_z + epsilon*(v/abs(v));
      
      //draw it to show what's happening
      mand_draw_ball(path.traps.back(), path.trap_colors.back());
      
    } while ( abs(path.traps.back().center - end_z) >= path.traps.back().radius );
  }
}






/****************************************************************************
 * main gui functions
 ****************************************************************************/

void IFSGui::pack_widget_upper_right(const Widget* w1, Widget* w2) {
  //figure out where it can go
  int desired_x,desired_y;
  if (w1 != NULL) {
    desired_x = w1->ul.x + w1->width;
    desired_y = w1->ul.y;
  } else {
    desired_x = 0;
    desired_y = 0;
  }
  
  
  //std::cout << "Packing widget of size " << w2->width << " " << w2->height << "\n";
  //std::cout << "Desired x: " << desired_x << "\n";
  
  //go through and check the other widgets to see how 
  //far down they obstruct this one
  int greatest_y_obstruction = 0;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w1) continue;
    if (widgets[i]->ul.x == desired_x && 
        widgets[i]->ul.y + widgets[i]->height > greatest_y_obstruction) {
      greatest_y_obstruction = widgets[i]->ul.y + widgets[i]->height;
      //std::cout << "Found widget " << i << " obstructs to height " << greatest_y_obstruction << "\n";
    }
  }
  if (greatest_y_obstruction + w2->height > main_window_height) {
    //std::cout << "Cannot pack widget -- too tall\n";
    return;
  }
  int y = (desired_y > greatest_y_obstruction ? desired_y : greatest_y_obstruction);
  
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
  
  //std::cout << "Packed widget at " << w2->ul << "\n";
  
  //record it in the list of widgets
  widgets.push_back(w2);
}


void IFSGui::detach_widget(Widget* w) {
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w) {
      widgets.erase(widgets.begin()+i);
      break;
    }
  }
  w->clear();
}



int IFSGui::get_rgb_color(double r, double g, double b) {
  XColor temp;
  temp.flags = DoRed | DoGreen | DoBlue;
  temp.red = (int)(r*65535);
  temp.green = (int)(g*65535);
  temp.blue = (int)(b*65535);
  if (XAllocColor(display, DefaultColormap(display, screen), &temp) == 0) {
    std::cout << "Color not found?\n";
  }
  return temp.pixel;
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
  int height_rest = display_height - 150;
  int x = (width_rest > height_rest ? height_rest : width_rest);
  
  if (window_mode == MANDLEBROT) {
    main_window_height = x + 100;
    main_window_width = x + mand_sidebar_size;
  } else if (window_mode == LIMIT) {
    main_window_height = x + 100;
    main_window_width = x + limit_sidebar_size;
  } else {
    main_window_height = x + 100;
    main_window_width = 2*x + mand_sidebar_size + limit_sidebar_size;
  }
  
  //compute the widths
  limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(x);
  
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
  
  //used to fill the text
  std::stringstream T; 
  
  //stuff for the IFS computations
  W_switch_to_limit = WidgetButton(this, "Switch to limit", -1, 20, &IFSGui::S_switch_to_limit);
  W_switch_to_mandlebrot = WidgetButton(this, "Switch to mandlebrot", -1,20, &IFSGui::S_switch_to_mandlebrot);
  W_switch_to_combined = WidgetButton(this, "Switch to combined", -1, 20, &IFSGui::S_switch_to_combined);
  
  W_point_title = WidgetText(this, "Current IFS status:", -1, 20);
  W_point_point = WidgetText(this, "Location: initializing", 500, 20);
  W_point_connected_check = WidgetCheck(this, "Connectedness", 105, 20, point_connected_check, &IFSGui::S_point_connected);
  W_point_connected_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_point_connected_decrease_depth);
  T.str(""); T << point_connected_depth;
  W_point_connected_depth_label = WidgetText(this, T.str(), -1, 20);
  W_point_connected_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_point_connected_increase_depth);
  W_point_connected_status = WidgetText(this, "initializing", -1, 20);
  
  W_point_contains_half_check = WidgetCheck(this, "Contains 1/2", 105, 20, point_contains_half_check, &IFSGui::S_point_contains_half);
  W_point_contains_half_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_point_contains_half_decrease_depth);
  T.str(""); T << point_contains_half_depth;
  W_point_contains_half_depth_label = WidgetText(this, T.str(), -1, 20);
  W_point_contains_half_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_point_contains_half_increase_depth);
  W_point_contains_half_status = WidgetText(this, "initializing", -1, 20);
  
  W_point_uv_words_check = WidgetCheck(this, "uv words", 105, 20, point_uv_words_check, &IFSGui::S_point_uv_words);
  W_point_uv_words_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_point_uv_words_decrease_depth);
  T.str(""); T << point_uv_words_depth;
  W_point_uv_words_depth_label = WidgetText(this, T.str(), -1, 20);
  W_point_uv_words_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_point_uv_words_increase_depth);
  W_point_uv_words_status = WidgetText(this, "initializing", x, 20);
  
  
  //if the limit set is shown:
  if (window_mode != MANDLEBROT) {
    W_limit_plot = WidgetDraw(this, x,x, &IFSGui::S_limit_draw);
    W_limit_depth_title = WidgetText(this, "Depth: ", -1, 20);
    W_limit_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_limit_decrease_depth);
    T.str(""); T << limit_depth;
    W_limit_depth_label = WidgetText(this, T.str(), -1, 20);
    W_limit_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_limit_increase_depth);
    W_limit_depth_auto = WidgetCheck(this, "Auto depth", -1, 20, limit_auto_depth, &IFSGui::S_limit_auto_depth);
    W_limit_chunky = WidgetCheck(this, "Chunky", -1, 20, limit_chunky, &IFSGui::S_limit_switch_chunky);
    W_limit_colors = WidgetCheck(this, "Colors", -1, 20, limit_colors, &IFSGui::S_limit_switch_colors);
    W_limit_zoom_title = WidgetText(this, "Zoom: ", -1, 20);
    W_limit_zoom_in = WidgetButton(this, "in", 30, 20, &IFSGui::S_limit_zoom_in);
    W_limit_zoom_out = WidgetButton(this, "out", 30, 20, &IFSGui::S_limit_zoom_out);
    W_limit_center_title = WidgetText(this, "(Click to center)", -1, 20);
    W_limit_uv_graph = WidgetCheck(this, "Plot uv graph", -1, 20, limit_uv_graph, &IFSGui::S_limit_uv_graph);
    W_limit_uv_graph_depth_title = WidgetText(this, "Depth:", -1, 20);
    W_limit_uv_graph_depth_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_limit_uv_graph_decrease_depth);
    T.str(""); T << limit_uv_graph_depth;
    W_limit_uv_graph_depth_label = WidgetText(this, T.str(), -1, 20);
    W_limit_uv_graph_depth_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_limit_uv_graph_increase_depth);
    
    pack_widget_upper_right(NULL, &W_limit_plot);
    if (window_mode == LIMIT) {
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_mandlebrot);
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_combined);
    }
    pack_widget_upper_right(&W_limit_plot, &W_limit_center_title);
    pack_widget_upper_right(&W_limit_plot, &W_limit_depth_title);
    pack_widget_upper_right(&W_limit_depth_title, &W_limit_depth_leftarrow);
    pack_widget_upper_right(&W_limit_depth_leftarrow, &W_limit_depth_label);
    pack_widget_upper_right(&W_limit_depth_label, &W_limit_depth_rightarrow);
    pack_widget_upper_right(&W_limit_plot, &W_limit_depth_auto);
    pack_widget_upper_right(&W_limit_plot, &W_limit_zoom_title);
    pack_widget_upper_right(&W_limit_zoom_title, &W_limit_zoom_in);
    pack_widget_upper_right(&W_limit_zoom_in, &W_limit_zoom_out);
    pack_widget_upper_right(&W_limit_plot, &W_limit_chunky);
    pack_widget_upper_right(&W_limit_plot, &W_limit_colors);
    pack_widget_upper_right(&W_limit_plot, &W_limit_uv_graph);
    pack_widget_upper_right(&W_limit_plot, &W_limit_uv_graph_depth_title);
    pack_widget_upper_right(&W_limit_uv_graph_depth_title, &W_limit_uv_graph_depth_leftarrow);
    pack_widget_upper_right(&W_limit_uv_graph_depth_leftarrow, &W_limit_uv_graph_depth_label);
    pack_widget_upper_right(&W_limit_uv_graph_depth_label, &W_limit_uv_graph_depth_rightarrow);
    
  }
  
  //if the mandlebrot set is shown:
  if (window_mode != LIMIT) {
    W_mand_plot = WidgetDraw(this, x,x, &IFSGui::S_mand_draw);
    W_mand_options_title = WidgetText(this, "Mandlebrot options:", -1, 20);
    W_mand_recenter = WidgetButton(this, "Recenter", -1, 20, &IFSGui::S_mand_recenter);
    W_mand_zoom_title = WidgetText(this, "Zoom: ", -1, 20);
    W_mand_zoom_in = WidgetButton(this, "in", -1, 20, &IFSGui::S_mand_zoom_in);
    W_mand_zoom_out = WidgetButton(this, "out", -1, 20, &IFSGui::S_mand_zoom_out);  
    W_mand_mesh_title = WidgetText(this, "Mesh size:", -1, 20);
    W_mand_mesh_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_decrease_mesh);
    std::stringstream T;  T.str("");  T << mand_pixel_group_size;
    W_mand_mesh_label = WidgetText(this, T.str(), -1, 20);
    W_mand_mesh_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_increase_mesh);
    W_mand_connected_check = WidgetCheck(this, "Connectedness:", 105, 20, (mand_connected ? 1 : 0), &IFSGui::S_mand_connected);
    W_mand_connected_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_connected_decrease_depth);
    T.str("");  T << mand_connected_depth;
    W_mand_connected_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_connected_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_connected_increase_depth);
    W_mand_contains_half_check = WidgetCheck(this, "Contains 1/2:", 105, 20, (mand_contains_half ? 1 : 0), &IFSGui::S_mand_contains_half);
    W_mand_contains_half_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_contains_half_decrease_depth);
    T.str("");  T << mand_contains_half_depth;
    W_mand_contains_half_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_contains_half_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_contains_half_increase_depth);
    W_mand_trap_check = WidgetCheck(this, "Traps:", 105, 20, (mand_trap ? 1 : 0), &IFSGui::S_mand_trap);
    W_mand_trap_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_trap_decrease_depth);
    T.str("");  T << mand_trap_depth;
    W_mand_trap_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_trap_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_trap_increase_depth);
    W_mand_mouse_label = WidgetText(this, "Mouse: initializing", 200, 20);
    
    W_mand_path_create_by_drawing_button = WidgetButton(this, "Draw path", -1, 20, &IFSGui::S_mand_path_create_by_drawing_button);
    W_mand_path_create_by_boundary_button = WidgetButton(this, "Find boundary path", -1, 20, &IFSGui::S_mand_path_create_by_boundary);
    W_mand_path_drawing_title = WidgetText(this, "(Click to draw path)", -1, 20);
    W_mand_path_finish_cancel_button = WidgetButton(this, "Cancel path", -1, 20, &IFSGui::S_mand_path_finish_cancel);
    W_mand_path_finish_path_button = WidgetButton(this, "Finish path", -1, 20, &IFSGui::S_mand_path_finish_path);
    W_mand_path_finish_loop_button = WidgetButton(this, "Finish loop", -1, 20, &IFSGui::S_mand_path_finish_loop);
    W_mand_path_tasks_title = WidgetText(this, "Path options:", -1, 20);
    W_mand_path_delete_button = WidgetButton(this, "Delete path", -2, 20, &IFSGui::S_mand_path_delete);
    W_mand_path_find_traps_button = WidgetButton(this, "Find traps along path", -1, 20, &IFSGui::S_mand_path_find_traps);
    W_mand_path_create_movie_button = WidgetButton(this, "Create movie along path", -1, 20, &IFSGui::S_mand_path_create_movie);
    W_mand_path_movie_length_title = WidgetText(this, "Movie length: ", -1, 20);
    W_mand_path_movie_decrease_length = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_path_movie_decrease_length);
    T.str(""); T << path.movie_length;
    W_mand_path_movie_length_label = WidgetText(this, T.str(), -1, 20);
    W_mand_path_movie_increase_length = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_path_movie_increase_length);
    W_mand_path_movie_with_mandlebrot = WidgetCheck(this, "Movie with mandlebrot", -1, 20, path.movie_with_mandlebrot, &IFSGui::S_mand_path_movie_with_mandlebrot);
    W_mand_path_find_uv_words_button = WidgetButton(this, "Find uv words along path", -1, 20, &IFSGui::S_mand_path_find_uv_words);
    
    if (window_mode == MANDLEBROT) {
      pack_widget_upper_right(NULL, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_combined);
    } else {
      pack_widget_upper_right(&W_limit_center_title, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_mandlebrot);
    }
    pack_widget_upper_right(&W_mand_plot, &W_mand_options_title);
    pack_widget_upper_right(&W_mand_plot, &W_mand_recenter);
    pack_widget_upper_right(&W_mand_plot, &W_mand_zoom_title);
    pack_widget_upper_right(&W_mand_zoom_title, &W_mand_zoom_in);
    pack_widget_upper_right(&W_mand_zoom_in, &W_mand_zoom_out);
    pack_widget_upper_right(&W_mand_plot, &W_mand_mesh_title);
    pack_widget_upper_right(&W_mand_mesh_title, &W_mand_mesh_leftarrow);
    pack_widget_upper_right(&W_mand_mesh_leftarrow, &W_mand_mesh_label);
    pack_widget_upper_right(&W_mand_mesh_label, &W_mand_mesh_rightarrow);
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
    pack_widget_upper_right(&W_mand_plot, &W_mand_mouse_label);
    if (currently_drawing_path) {
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_drawing_title);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_cancel_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_path_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_finish_loop_button);
    } else if (path.is_valid) {
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_tasks_title);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_delete_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_traps_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_movie_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_length_title);
      pack_widget_upper_right(&W_mand_path_movie_length_title, &W_mand_path_movie_decrease_length);
      pack_widget_upper_right(&W_mand_path_movie_decrease_length, &W_mand_path_movie_length_label);
      pack_widget_upper_right(&W_mand_path_movie_length_label, &W_mand_path_movie_increase_length);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_with_mandlebrot);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_uv_words_button);
    } else {
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_by_drawing_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_create_by_boundary_button);
    }
  }
  
  //put the IFS data on the bottom
  pack_widget_upper_right(NULL, &W_point_title);
  pack_widget_upper_right(NULL, &W_point_point);
  pack_widget_upper_right(NULL, &W_point_connected_check);
  pack_widget_upper_right(&W_point_connected_check, &W_point_connected_leftarrow);
  pack_widget_upper_right(&W_point_connected_leftarrow, &W_point_connected_depth_label);
  pack_widget_upper_right(&W_point_connected_depth_label, &W_point_connected_rightarrow);
  pack_widget_upper_right(&W_point_connected_rightarrow, &W_point_connected_status);
  pack_widget_upper_right(NULL, &W_point_contains_half_check);
  pack_widget_upper_right(&W_point_contains_half_check, &W_point_contains_half_leftarrow);
  pack_widget_upper_right(&W_point_contains_half_leftarrow, &W_point_contains_half_depth_label);
  pack_widget_upper_right(&W_point_contains_half_depth_label, &W_point_contains_half_rightarrow);
  pack_widget_upper_right(&W_point_contains_half_rightarrow, &W_point_contains_half_status);
  pack_widget_upper_right(NULL, &W_point_uv_words_check);
  pack_widget_upper_right(&W_point_uv_words_check, &W_point_uv_words_leftarrow);
  pack_widget_upper_right(&W_point_uv_words_leftarrow, &W_point_uv_words_depth_label);
  pack_widget_upper_right(&W_point_uv_words_depth_label, &W_point_uv_words_rightarrow);
  pack_widget_upper_right(&W_point_uv_words_rightarrow, &W_point_uv_words_status);
  
  
  
  //invalidate the grids and stuff
  if (window_mode != LIMIT) {
    mand_reset_mesh();
  }
  
  //draw all the widgets
  for (int i=0; i<(int)widgets.size(); ++i) {
    widgets[i]->initial_draw();
  }
  
  //plot the limit set
  if (window_mode != MANDLEBROT) draw_limit();
  if (window_mode != LIMIT) draw_mand();
  
  //get the point data
  recompute_point_data();
  
}



void IFSGui::main_loop() {
  XEvent e;
  while (true) {
    XNextEvent(display, &e);
    //if it was the keyboard, we deal with it here
    if (e.type == KeyPress) {
      if(XLookupKeysym(&e.xkey, 0) == XK_q){ // left arrow
        break;
      }
     
    //if it involves the mouse, we find the appropriate 
    //widget to send it off to
    } else if (e.type == ButtonPress || e.type == MotionNotify) {
      for (int i=0; i<(int)widgets.size(); ++i) {
        if (widgets[i]->contains_pixel( e.xbutton.x, e.xbutton.y) &&
            widgets[i]->click_signal != NULL) {
          (this->*(widgets[i]->click_signal))(&e);
          break;
        }
      }
    }
  }
}




void IFSGui::launch(IFSWindowMode m, const cpx& c) {
  
  //set the initial window mode
  window_mode = m;
  
  //set the ifs
  IFS.set_params(c,c);
  
  //set the initial settings for limit and mandlebrot
  limit_ll = cpx(-1, -1.5);
  limit_ur = cpx(2, 1.5);
  limit_depth = 12;
  limit_auto_depth = false;
  limit_chunky = true;
  limit_colors = true;
  limit_uv_graph = false;
  limit_uv_graph_depth = 3;
  limit_marked_points.resize(3);
  limit_marked_points[0] = cpx(0,0);
  limit_marked_points[1] = cpx(0.5,0);
  limit_marked_points[2] = cpx(1.0,0);
  
  mand_ll = cpx(-1,-1);
  mand_ur = cpx(1,1);
  mand_pixel_group_size = 4;
  mand_connected = true;
  mand_connected_depth = 13;
  mand_contains_half = false;
  mand_contains_half_depth = 16;
  mand_trap = false;
  mand_trap_depth = 20;
  
  point_connected_check = true;
  point_connected_depth = 18;
  point_is_connected = false;
  point_contains_half_check = true;
  point_contains_half_depth = 18;
  point_is_contains_half = false;
  //point_trap = false;
  point_uv_words_check = false;
  point_uv_words_depth = 18;
  point_uv_words.resize(0);
  
  currently_drawing_path = false;
  
  
  //set up the graphics
  display = XOpenDisplay(NULL);
  screen = DefaultScreen(display);
  if (display == NULL) {
    std::cout << "Failed to open display\n";
    return;
  }
  main_window_initialized = false;
  
  //reset (set) the window
  limit_sidebar_size = 130;
  mand_sidebar_size = 170;
  reset_and_pack_window();
  
  //go for it
  main_loop();
  
  //clean up at the end
  XCloseDisplay(display);
  
}