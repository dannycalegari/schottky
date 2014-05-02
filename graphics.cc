#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

#include "graphics.h"

XGraphics::XGraphics() {
/* //maybe it's better if this does nothing?
  border_width = 4;
  display=XOpenDisplay(NULL);
  screen_num = DefaultScreen(display);  
  display_width = DisplayWidth(display, screen_num);
  display_height = DisplayHeight(display, screen_num);
  width = 800;
  height = 800;
  win = XCreateSimpleWindow(display, RootWindow(display, screen_num), 0, 0, width, 
                            height, border_width, BlackPixel(display, screen_num), WhitePixel(display, screen_num));
  XSelectInput(display, win, ExposureMask | 
                             KeyPressMask | 
                             ButtonPressMask | 
                             PointerMotionMask |
                             StructureNotifyMask);
  gc = DefaultGC(display, screen_num);
  screen_colormap = DefaultColormap(display, screen_num);
  XSetForeground(display, gc, BlackPixel(display, screen_num));
  XSetBackground(display, gc, WhitePixel(display, screen_num));
  XMapWindow(display, win);
  color_list.clear();
  
  line_thickness = 1;
  
  scale = 1.0;
  translate = Point2d<float>(0.0, 0.0);
  
  while (true) {  //wait until the window is actually mapped
    XEvent e;
    XNextEvent(display, &e);
    if (e.type == Expose) break;
  }
  setup_font();
  */
  initialized = false;
}

XGraphics::XGraphics(int w, int h, float s, const Point2d<float>& t) {
  initialize(w,h,s,t);
}



void XGraphics::initialize(int w, int h, float s, const Point2d<float>& t) {
  border_width = 4;
  display=XOpenDisplay(NULL);
  screen_num = DefaultScreen(display);  
  display_width = DisplayWidth(display, screen_num);
  display_height = DisplayHeight(display, screen_num);
  width = w;
  height = h;
  win = XCreateSimpleWindow(display, RootWindow(display, screen_num), 0, 0, width, 
                            height, border_width, BlackPixel(display, screen_num), WhitePixel(display, screen_num));
  XSelectInput(display, win, ExposureMask | 
                             KeyPressMask | 
                             ButtonPressMask | 
                             PointerMotionMask |
                             StructureNotifyMask);
  gc = DefaultGC(display, screen_num);
  screen_colormap = DefaultColormap(display, screen_num);
  XSetForeground(display, gc, BlackPixel(display, screen_num));
  XSetBackground(display, gc, WhitePixel(display, screen_num));
  XMapWindow(display, win);
  color_list.clear();
  
  line_thickness = 1;
  
  scale = s;
  translate = t;
  
  initialized = true;
  
  while (true) {  //wait until the window is actually mapped
    XEvent e;
    XNextEvent(display, &e);
    if (e.type == Expose) break;
  }
  setup_font();
}






XGraphics::~XGraphics() {
  if (initialized) XCloseDisplay(display);
}


int XGraphics::get_color(std::string c) {
  std::map<std::string, int>::iterator it;
  it = color_list.find(c);
  if (it == color_list.end()) {
    XColor temp;
    Status rc = XAllocNamedColor(display, screen_colormap, c.c_str(), &temp, &temp);
    if (rc == 0) {
      std::cout << "Color error\n";
      return 0;
    }
    color_list[c] = temp.pixel;
  }
  return color_list[c];
}
    

int XGraphics::get_rgb_color(double r, double g, double b) {
  XColor temp;
  temp.flags = DoRed | DoGreen | DoBlue;
  temp.red = (int)(r*65535);
  temp.green = (int)(g*65535);
  temp.blue = (int)(b*65535);
  if (XAllocColor(display, screen_colormap, &temp) == 0) {
    std::cout << "Color not found?\n";
  }
  return temp.pixel;
}



void XGraphics::flush() {
  XFlush(display);
}

void XGraphics::list_fonts() {
  std::cout << "Listing fonts: \n";
  int num_returned;
  char** names = XListFonts(display, "*", 10000, &num_returned);
  for (int i=0; i<num_returned; ++i) {
    std::cout << names[i] << "\n";
  }
}


void XGraphics::setup_font(void){
  //just use the fixed font, because it looks fine
  font = XLoadQueryFont(display, "fixed");
  XSetFont(display, gc, font->fid);
  //const char * fontname = "-*-georgia-*-r-*-*-14-*-*-*-*-*-*-*";
  //const char * fontname = "-*-times-*-r-*-*-16-*-*-*-*-*-*-*";
  //font = XLoadQueryFont (display, fontname);
  //if (! font) {
  //  font = XLoadQueryFont (display, "fixed");
  //  std::cout << "couldn't find font!\n";
  //}
  //XSetFont (display, gc, font->fid);
}

Point2d<int> XGraphics::text_offset_left(const std::string& S) {
  XCharStruct te;
  int fdir, fdescent, fascent;
  XTextExtents(font, S.c_str(), S.size(), &fdir, &fascent, &fdescent, &te);
  int height_offset = (te.ascent-te.descent)/2;
  if ( (te.ascent-te.descent)%2 == 1) --height_offset;
  return Point2d<int>(-te.lbearing, -height_offset);
}

//returns the vector to ADD to your vector such that the center of the 
//text extents gets placed directly over your vector
Point2d<int> XGraphics::text_offset_center(const std::string& S) {
  XCharStruct te;
  int fdir, fdescent, fascent;
  XTextExtents(font, S.c_str(), S.size(), &fdir, &fascent, &fdescent, &te);
  int height_offset = (te.ascent-te.descent)/2;
  //if ( (te.ascent-te.descent)%2 == 1) --height_offset;
  int width_offset = (te.rbearing - te.lbearing)/2;
  return Point2d<int>(-width_offset, -height_offset);
}




void XGraphics::erase_field(void){
  XClearWindow(display, win);
}

Point2d<int> XGraphics::mouse_location(){
//    Bool result;
  Window window_returned;
  int root_x, root_y;
  int win_x, win_y;
  unsigned int mask_return;
  Point2d<int> p;
    
  XQueryPointer(display, win, &window_returned,
          &window_returned, &root_x, &root_y, &win_x, &win_y,
          &mask_return);
  p.x=win_x;
  p.y=height-win_y;
  return p;
}

void XGraphics::draw_point(const Point2d<int>& p, long col){
  XSetForeground(display, gc, col);
  XDrawPoint(display, win, gc, p.x, height-p.y);
}

void XGraphics::draw_dot(const Point2d<int>& p, long col) {
    XSetForeground(display, gc, col);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid);
    XDrawArc(display, win, gc, p.x-3, height-p.y-3, 6, 6, 0, 23040);
    XDrawArc(display, win, gc, p.x-2, height-p.y-2, 4, 4, 0, 23040);
    XDrawArc(display, win, gc, p.x-1, height-p.y-1, 2, 2, 0, 23040);
}

void XGraphics::draw_box(const Point2d<int>& p, int w, long col) {
  if(w==1){
    draw_point(p, col);
  } else {
    XPoint corners[4];
    XSetForeground(display, gc, col);
    corners[0].x = p.x;
    corners[0].y = height-p.y-w;
    corners[1].x = corners[0].x+w;
    corners[1].y = corners[0].y;
    corners[2].x = corners[1].x;
    corners[2].y = corners[1].y+w;
    corners[3].x = corners[0].x;
    corners[3].y = corners[0].y+w;
    XFillPolygon(display, win, gc, corners, 4, Convex, CoordModeOrigin);
  }
}


void XGraphics::draw_point(const Point2d<float>& p, long col){
  XSetForeground(display, gc, col);
  Point2d<int> real_p(scale*p.x + translate.x, scale*p.y +translate.y);
  XDrawPoint(display, win, gc, real_p.x, height-real_p.y);
}

void XGraphics::draw_line(const Point2d<int>& p1, const Point2d<int>& p2, long col) {
  XSetForeground(display, gc, col);
  XSetLineAttributes(display, gc, 1.5, LineSolid, 1, 1);
  XDrawLine(display, win, gc, p1.x, height-p1.y, p2.x, height-p2.y);
}

void XGraphics::draw_line(const Point2d<float>& p1, const Point2d<float>& p2, long col) {
  XSetForeground(display, gc, col);
  XSetLineAttributes(display, gc, line_thickness, LineSolid, 1, 1);
  Point2d<int> p1_real((int)(p1.x*scale + translate.x), (int)(p1.y*scale + translate.y));
  Point2d<int> p2_real((int)(p2.x*scale + translate.x), (int)(p2.y*scale + translate.y));
  XDrawLine(display, win, gc, p1_real.x, height-p1_real.y, p2_real.x, height-p2_real.y);
}


void XGraphics::draw_line(const Point2d<float>& p1, const Point2d<float>& p2, long col, int thickness) {
  XSetForeground(display, gc, col);
  XSetLineAttributes(display, gc, thickness, LineSolid, CapButt, JoinMiter);
  Point2d<int> p1_real((int)(p1.x*scale + translate.x), (int)(p1.y*scale + translate.y));
  Point2d<int> p2_real((int)(p2.x*scale + translate.x), (int)(p2.y*scale + translate.y));
  XDrawLine(display, win, gc, p1_real.x, height-p1_real.y, p2_real.x, height-p2_real.y);
  XSetLineAttributes(display, gc, line_thickness, LineSolid, CapButt, JoinMiter);
}
  
void XGraphics::draw_arrowed_labeled_line(const Point2d<float>& p1, 
                                          const Point2d<float>& p2, 
                                          long col, 
                                          int thickness,
                                          const std::string& label) {
  //draw the main line
  XSetForeground(display, gc, col);
  XSetLineAttributes(display, gc, thickness, LineSolid, CapButt, JoinMiter);
  Point2d<int> p1_real((int)(p1.x*scale + translate.x), (int)(p1.y*scale + translate.y));
  Point2d<int> p2_real((int)(p2.x*scale + translate.x), (int)(p2.y*scale + translate.y));
  XDrawLine(display, win, gc, p1_real.x, height-p1_real.y, p2_real.x, height-p2_real.y);
  XSetLineAttributes(display, gc, line_thickness, LineSolid, CapButt, JoinMiter);
  
  //draw the arrow
  Point2d<float> diff = p2-p1;
  Point2d<float> perp(-diff.y, diff.x);
  float scale_to_1 = 1.0/sqrt(diff.x*diff.x + diff.y*diff.y);
  Point2d<float> scaled_diff = scale_to_1*diff;
  perp = scale_to_1*perp;
  Point2d<float> center = p1 + (float)0.5*diff;
  Point2d<float> side1 = center - (float)0.015*scaled_diff + (float)0.015*perp;
  Point2d<float> side2 = center - (float)0.015*scaled_diff - (float)0.015*perp;
  draw_line(center, side1, col, thickness);
  draw_line(center, side2, col, thickness);
  
  //draw the label
  if (label.size() > 0) {
    Point2d<float>text_center = center +(float)0.02*scaled_diff + (float)0.02*perp;
    draw_text_centered(text_center, label, col);
  }
}


void XGraphics::draw_square(int x, int y, int z, long col){
  Point2d<int> p1(x,y);
  Point2d<int> p2(x+z,y);
  Point2d<int> p3(x+z,y+z);
  Point2d<int> p4(x,y+z);
  this->draw_line(p1, p2, col);
  this->draw_line(p2, p3, col);
  this->draw_line(p3, p4, col);
  this->draw_line(p4, p1, col);
  //XSetForeground(display, gc, col);
  //XSetLineAttributes(display, gc, 2, LineSolid, 1, 1);
  //XDrawLine(display, win, gc, x, y, x+z, y);
  //XDrawLine(display, win, gc, x+z, y, x+z, y+z);
  //XDrawLine(display, win, gc, x+z, y+z, x, y+z);
  //XDrawLine(display, win, gc, x, y+z, x, y);
}

void XGraphics::draw_rectangle(int x, int y, int zx, int zy, long col) {
  Point2d<int> p1(x,y);
  Point2d<int> p2(x+zx,y);
  Point2d<int> p3(x+zx,y+zy);
  Point2d<int> p4(x,y+zy);
  this->draw_line(p1, p2, col);
  this->draw_line(p2, p3, col);
  this->draw_line(p3, p4, col);
  this->draw_line(p4, p1, col);
  //XSetForeground(display, gc, BlackPixel(display, screen_num));
  //XSetLineAttributes(display, gc, 2, LineSolid, 1, 1);
  //XDrawLine(display, win, gc, x, y, x+zx, y);
  //XDrawLine(display, win, gc, x+zx, y, x+zx, y+zy);
  //XDrawLine(display, win, gc, x+zx, y+zy, x, y+zy);
  //XDrawLine(display, win, gc, x, y+zy, x, y);
}

void XGraphics::draw_filled_rectangle(int x, int y, int zx, int zy, long col) {
  XSetForeground(display, gc, col);
  XFillRectangle(display, win, gc, x, height-(y+zy), zx, zy);
}


void XGraphics::draw_filled_rectangle(const Point2d<int>& p, int w, int h, long col) {
  XSetForeground(display, gc, col);
  XFillRectangle(display, win, gc, p.x, height-p.y-h, w, h);
}


void XGraphics::draw_box_radius(Point2d<float>& center, float radius, long col) {
  Point2d<int> LL_real = Point2d<int>( (int)((center.x-radius)*scale + translate.x),
                                       (int)((center.y-radius)*scale + translate.y));
  int z = (int)(2*radius*scale);
  XSetForeground(display, gc, col);
  XFillRectangle(display, win, gc, LL_real.x, height-(LL_real.y+z), z, z);
  draw_point(center, col);
}


void XGraphics::draw_filled_polygon(const std::vector<Point2d<float> >& points, int col) {
  std::vector<XPoint> real_points(points.size());
  for (int i=0; i<(int)points.size(); ++i) {
    real_points[i].x = (int)(points[i].x*scale + translate.x);
    real_points[i].y = height - (int)(points[i].y*scale + translate.y);
  }
  XSetForeground(display, gc, col);
  XFillPolygon(display, win, gc, &real_points[0], points.size(), Convex, CoordModeOrigin);
}


void XGraphics::draw_faint_line(const Point2d<int>& p1, 
                                const Point2d<int>& p2, 
                                long col){
  XSetForeground(display, gc, (long) 0xDDDDDD);
  XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
  XDrawLine(display, win, gc, p1.x, p1.y, p2.x, p2.y);
}

void XGraphics::erase_circle(const Point2d<int>& p, int r){
	XSetForeground(display, gc, 0xFFFFFF);
  XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
	XSetFillStyle(display, gc, FillSolid);
  XFillArc(display, win, gc, p.x-r, p.y-r, 2*r, 2*r, 0, 23040);
}

void XGraphics::draw_circle(const Point2d<int>& p, int r, long col){
    XSetForeground(display, gc, col);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid); 
    XDrawArc(display, win, gc, (p.x-r), (height-p.y)-r, 2*r, 2*r, 0, 23040);
}

void XGraphics::draw_circle(const Point2d<float>& p, int r, long col){
    XSetForeground(display, gc, col);
    XSetLineAttributes(display, gc, line_thickness, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid); 
    Point2d<int> real_p((int)(p.x*scale + translate.x), (int)(p.y*scale + translate.y));
    XDrawArc(display, win, gc, (real_p.x-r), (height-(real_p.y+r)), 2*r, 2*r, 0, 23040);
}

void XGraphics::draw_disk(const Point2d<int>& p, double r, long col) {
    XSetForeground(display, gc, col);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid);
    XFillArc(display, win, gc, p.x-r, height-p.y-r, 2*r, 2*r, 0, 23040);
}


void XGraphics::draw_concentric_circles(const Point2d<int>& p, int r, long col){
  int s;
  for(s=3;s>0;s--){
    XSetForeground(display, gc, col*s);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid);
    XDrawArc(display, win, gc, p.x-r*s+1, p.y-r*s+1, 2*r*s-2, 2*r*s-2, 0, 23040);
  }
}

void XGraphics::draw_path(const std::vector<Point2d<int> >& L, long col){
  int i;
  for(i=1; i<(int)L.size(); i++){
    draw_line(L[i-1],L[i],col);
  }
}

void XGraphics::draw_text(const Point2d<int>& p, std::stringstream &T, long col){
  std::string S;
  XSetForeground(display, gc, col);
  S=T.str();
  XDrawString(display,win,gc,p.x,height-p.y,S.c_str(),strlen(S.c_str()));
}

void XGraphics::draw_text(const Point2d<int>& p, std::string &S, long col){
  XSetForeground(display, gc, col);
  XDrawString(display,win,gc,p.x,height-p.y,S.c_str(),strlen(S.c_str()));
}

void XGraphics::draw_text(const Point2d<float>& p, std::string &S, long col) {
  XSetForeground(display, gc, col);
  Point2d<int> real_p = Point2d<int>((int)(p.x*scale + translate.x), 
                                     (int)(p.y*scale + translate.y));
  Point2d<int> offset = text_offset_left(S);
  XDrawString(display, win, gc, real_p.x+offset.x, height-(real_p.y+offset.y), S.c_str(), S.size());
}
  
void XGraphics::draw_text_centered(const Point2d<float>& p, const std::string &S, long col) {
  XSetForeground(display, gc, col);
  Point2d<int> real_p( (int)(p.x*scale + translate.x), (int)(p.y*scale + translate.y) );
  Point2d<int> offset = text_offset_center(S);
  XDrawString(display, win, gc, real_p.x+offset.x, height-(real_p.y+offset.y), S.c_str(), S.size()); 
  //draw_point(p, col);
}
  
void XGraphics::draw_label(const Point2d<int>& p, int i, long col){
  XSetForeground(display, gc, col);
  std::stringstream T;
  T << i;
  Point2d<int> draw_loc;
  draw_loc.x = p.x+5;
  draw_loc.y = p.y-5;
  this->draw_text(draw_loc, T, col);
}

std::string XGraphics::wait_for_key() {
  XFlush(display);
  bool finished=false;
  while(finished==false){ 
    XNextEvent(display, &report);
    if (report.type != KeyPress) continue; //ignore the mouse
    if(XLookupKeysym(&report.xkey, 0) == XK_Left){ // left arrow
      return "LA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Right){ // right arrow
      return "RA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Up){    // up arrow
      return "UA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Down){  // down arrow
      return "DA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_equal){ 
      return "=";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_minus){ 
      return "-";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_o){ 
      return "o";
    } else {
      return "";
    }
  }
  return "u";
}

void XGraphics::get_next_event(XEvent& xe) {
  XNextEvent(display, &xe);
}


