#include "ifs_gui.h"










/***************************************************************************
 * widget functions
 ***************************************************************************/
bool Widget::contains_pixel(int x, int y) {
  return (ul.x <= x) && (x < ul.x + width) && (ul.y <= y) && (y < ul.y + height);
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
  std::cout << "Drawing the drawing area\n";
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
  std::cout << "Drawing string: " << text << "\n";
  std::cout << "At position: " << ul.x << " " << ul.y << "\n";
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
}

void IFSGui::S_switch_to_mandlebrot(XEvent* e) {
}

void IFSGui::S_switch_to_combined(XEvent* e) {
}
 
//limit set
void IFSGui::S_limit_draw(XEvent* e) {
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
}

void IFSGui::S_limit_zoom_out(XEvent* e) {
}

void IFSGui::S_limit_recenter(XEvent* e) {
}

//mandlebrot
void IFSGui::S_mand_draw(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  int widget_x = e->xbutton.x - W_mand_plot.ul.x;
  int widget_y = e->xbutton.y - W_mand_plot.ul.y;
  cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
  IFS.set_params(c,c);
  draw_mand();
  if (window_mode != MANDLEBROT) draw_limit();
}



void IFSGui::S_mand_connected(XEvent* e) {}
void IFSGui::S_mand_connected_increase_depth(XEvent* e) {}
void IFSGui::S_mand_connected_decrease_depth(XEvent* e) {}
void IFSGui::S_mand_contains_half(XEvent* e) {}
void IFSGui::S_mand_contains_half_increase_depth(XEvent* e) {}
void IFSGui::S_mand_contains_half_decrease_depth(XEvent* e) {}
void IFSGui::S_mand_trap(XEvent* e) {}
void IFSGui::S_mand_trap_increase_depth(XEvent* e) {}
void IFSGui::S_mand_trap_decrease_depth(XEvent* e) {}

//point
void IFSGui::S_point_connected(XEvent* e) {}
void IFSGui::S_point_contains_half(XEvent* e) {}
void IFSGui::S_point_trap(XEvent* e) {}
void IFSGui::S_point_uv_words(XEvent* e) {}








/**************************************************************************
 * back end functions
 * ************************************************************************/
Point2d<int> IFSGui::limit_cpx_to_pixel(const cpx& c) {
  int x = int( (c.real() - limit_ll.real()) / limit_pixel_width );
  int real_y = int( (c.imag() - limit_ll.imag()) / limit_pixel_width );
  int x11_y = W_limit_plot.height - real_y;
  return Point2d<int>(x, x11_y);
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
  
  int blue_color = get_rgb_color(0,0.6,1.0);
  int yellow_color = get_rgb_color(1.0,0.6,0);
  
  Ball initial_ball(0.5,(IFS.z-1.0)/2.0,(1.0-IFS.w)/2.0,min_r);
  std::vector<std::pair<bool,Ball> > stack(0);
  stack.push_back(std::make_pair(false, initial_ball));
  while (stack.size() > 0) {
    std::pair<bool,Ball> b = stack.back();
    stack.pop_back();
    //if the ball is disjoint from the window, we might as well 
    //get rid of it
    if (!b.first && b.second.is_disjoint(limit_ll, limit_ur)) continue;
    if (b.second.word_len >= limit_depth) {
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
  LW.redraw();
}


cpx IFSGui::mand_pixel_group_to_cpx(const Point2d<int>& p) {
  double r = (p.x + 0.5)*mand_pixel_group_width + mand_ll.real();
  double i = mand_ur.imag() - (p.y + 0.5)*mand_pixel_group_width;
  return cpx(r,i);
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
    return get_rgb_color(1.0, double(p.z)/100, 0.0);
  } else if (mand_contains_half && p.y > 0) {
    return get_rgb_color(0.1, double(p.y)/100, 0.1);
  } else if (mand_connected && p.x >= 0) {
    return p.x*0x000001;
  } else {
    return WhitePixel(display, screen);
  }
}


void IFSGui::draw_mand() {
  ifs temp_IFS;
  Widget& MW = W_mand_plot;
  int num_pixel_groups = MW.width / mand_pixel_group_size;
  
  std::vector<std::pair<Bitword,Bitword> > uv_words;
  std::vector<Ball> TLB;
  bool found_TLB = false;
  double TLB_neighborhood;
  if (mand_trap && !mand_grid_trap_valid) {
    temp_IFS.TLB_and_uv_words_for_region(TLB, uv_words, TLB_neighborhood,
                                          mand_ll, mand_ur,
                                          0, mand_trap_depth, 0);
    found_TLB = (TLB.size() != 0);
  }
  
  for (int i=0; i<(int)num_pixel_groups; ++i) {
    for (int j=0; j<(int)num_pixel_groups; ++j) {
      
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
        mand_data_grid[i][j].z = temp_IFS.check_TLB(TLB,trap_radius,TLB_neighborhood,mand_trap_depth);
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
  if (mand_contains_half && !mand_grid_contains_half_valid) mand_grid_connected_valid = true;
  if (mand_trap && !mand_grid_trap_valid) mand_grid_connected_valid = true;
  
  //now draw the highlighted point
  Point2d<int> h = mand_cpx_to_pixel(IFS.z);
  int rcol = get_rgb_color(1.0,0.1,0.0);
  XSetForeground(display, MW.gc, rcol);
  XFillArc(display, MW.p, MW.gc, h.x, h.y, 4, 4, 23040, 23040);
  
  MW.redraw();
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
  
  
  std::cout << "Packing widget of size " << w2->width << " " << w2->height << "\n";
  std::cout << "Desired x: " << desired_x << "\n";
  
  //go through and check the other widgets to see how 
  //far down they obstruct this one
  int greatest_y_obstruction = 0;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w1) continue;
    if (widgets[i]->ul.x == desired_x && 
        widgets[i]->ul.y + widgets[i]->height > greatest_y_obstruction) {
      greatest_y_obstruction = widgets[i]->ul.y + widgets[i]->height;
      std::cout << "Found widget " << i << " obstructs to height " << greatest_y_obstruction << "\n";
    }
  }
  if (greatest_y_obstruction + w2->height > main_window_height) {
    std::cout << "Cannot pack widget -- too tall\n";
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
  
  std::cout << "Packed widget at " << w2->ul << "\n";
  
  //record it in the list of widgets
  widgets.push_back(w2);
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
  int height_rest = display_height - 100;
  int x = (width_rest > height_rest ? height_rest : width_rest);
  
  if (window_mode == MANDLEBROT) {
    main_window_height = x;// + 200;
    main_window_width = x + mand_sidebar_size;
  } else if (window_mode == LIMIT) {
    main_window_height = x;// + 200;
    main_window_width = x + limit_sidebar_size;
  } else {
    main_window_height = x;// + 200;
    main_window_width = 2*x + mand_sidebar_size + limit_sidebar_size;
  }
  
  //compute the widths
  limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(x);
  mand_pixel_group_width = mand_pixel_group_size*( (mand_ur.real() - mand_ll.real()) / double(x) );
  mand_pixel_width = (mand_ur.real() - mand_ll.real()) / double(x) ;
  
  //invalidate the grids and stuff
  if (window_mode != LIMIT) {
    int npg = x / mand_pixel_group_size;
    mand_data_grid.resize(npg);
    for (int i=0; i<npg; ++i) mand_data_grid[i] = std::vector<Point3d<int> >(npg, Point3d<int>(-1,-1,-1));
    mand_grid_connected_valid = false;
    mand_grid_contains_half_valid = false;
    mand_grid_trap_valid = false;
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
  W_switch_to_limit = WidgetButton(this, "Switch to limit", -1, 20, &IFSGui::S_switch_to_limit);
  W_switch_to_mandlebrot = WidgetButton(this, "Switch to mandlebrot", -1,20, &IFSGui::S_switch_to_mandlebrot);
  W_switch_to_combined = WidgetButton(this, "Switch to combined", -1, 20, &IFSGui::S_switch_to_combined);
  
  W_point_title = WidgetText(this, "Highlighted IFS options:", -1, 20);
  W_point_connected_check = WidgetCheck(this, "Connectedness", -1, 20, point_connected, &IFSGui::S_point_connected);
  W_point_contains_half_check = WidgetCheck(this, "Contains 1/2", -1, 20, point_contains_half, &IFSGui::S_point_contains_half);
  W_point_trap_check = WidgetCheck(this, "Find trap", -1, 20, point_trap, &IFSGui::S_point_trap);
  W_point_uv_word_check = WidgetCheck(this, "Find uv words", -1, 20, point_uv_words, &IFSGui::S_point_trap);
  
  
  if (window_mode != MANDLEBROT) {
    W_limit_plot = WidgetDraw(this, x,x, &IFSGui::S_limit_draw);
    W_limit_depth_title = WidgetText(this, "Depth: ", -1, 20);
    W_limit_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_limit_decrease_depth);
    std::stringstream T; T.str(""); T << limit_depth;
    W_limit_depth_label = WidgetText(this, T.str().c_str(), -1, 20);
    W_limit_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_limit_increase_depth);
    W_limit_chunky = WidgetCheck(this, "Chunky", -1, 20, limit_chunky, &IFSGui::S_limit_switch_chunky);
    W_limit_colors = WidgetCheck(this, "Colors", -1, 20, limit_colors, &IFSGui::S_limit_switch_colors);
    W_limit_zoom_title = WidgetText(this, "Zoom: ", -1, 20);
    W_limit_zoom_in = WidgetButton(this, "in", 30, 20, &IFSGui::S_limit_zoom_in);
    W_limit_zoom_out = WidgetButton(this, "out", 30, 20, &IFSGui::S_limit_zoom_out);
    W_limit_center_title = WidgetText(this, "(Click to center)", -1, 20);
    
    pack_widget_upper_right(NULL, &W_limit_plot);
    if (window_mode == LIMIT) {
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_mandlebrot);
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_combined);
    }
    pack_widget_upper_right(&W_limit_plot, &W_limit_depth_title);
    pack_widget_upper_right(&W_limit_depth_title, &W_limit_depth_leftarrow);
    pack_widget_upper_right(&W_limit_depth_leftarrow, &W_limit_depth_label);
    pack_widget_upper_right(&W_limit_depth_label, &W_limit_depth_rightarrow);
    pack_widget_upper_right(&W_limit_plot, &W_limit_zoom_title);
    pack_widget_upper_right(&W_limit_zoom_title, &W_limit_zoom_in);
    pack_widget_upper_right(&W_limit_zoom_in, &W_limit_zoom_out);
    pack_widget_upper_right(&W_limit_plot, &W_limit_chunky);
    pack_widget_upper_right(&W_limit_plot, &W_limit_colors);
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
    W_mand_plot = WidgetDraw(this, x,x, &IFSGui::S_mand_draw);
    W_mand_connected_check = WidgetCheck(this, "Connectedness:", 105, 20, (mand_connected ? 1 : 0), &IFSGui::S_mand_connected);
    W_mand_connected_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_connected_decrease_depth);
    std::stringstream T;  T.str("");  T << mand_connected_depth;
    W_mand_connected_depth_label = WidgetText(this, T.str().c_str(), -1, 20);
    W_mand_connected_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_connected_increase_depth);
    W_mand_contains_half_check = WidgetCheck(this, "Contains 1/2:", 105, 20, (mand_contains_half ? 1 : 0), &IFSGui::S_mand_contains_half);
    W_mand_contains_half_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_contains_half_decrease_depth);
    T.str("");  T << mand_contains_half_depth;
    W_mand_contains_half_depth_label = WidgetText(this, T.str().c_str(), -1, 20);
    W_mand_contains_half_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_contains_half_increase_depth);
    W_mand_trap_check = WidgetCheck(this, "Traps:", 105, 20, (mand_trap ? 1 : 0), &IFSGui::S_mand_trap);
    W_mand_trap_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_trap_decrease_depth);
    T.str("");  T << mand_trap_depth;
    W_mand_trap_depth_label = WidgetText(this, T.str().c_str(), -1, 20);
    W_mand_trap_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_trap_increase_depth);
    
    if (window_mode == MANDLEBROT) {
      pack_widget_upper_right(NULL, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_combined);
    } else {
      pack_widget_upper_right(&W_limit_depth_rightarrow, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_mandlebrot);
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
  
  //plot the limit set
  if (window_mode != MANDLEBROT) draw_limit();
  if (window_mode != LIMIT) draw_mand();
  
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
  limit_chunky = true;
  limit_colors = true;
  
  mand_ll = cpx(-1,-1);
  mand_ur = cpx(1,1);
  mand_pixel_group_size = 2;
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
  limit_sidebar_size = 130;
  mand_sidebar_size = 200;
  reset_and_pack_window();
  
  //go for it
  main_loop();
  
  //clean up at the end
  XCloseDisplay(display);
  
}