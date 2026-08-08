#include <cstdlib>
#include <cmath>

#include "ifs_gui.h"


#include "movie.h"
#include "figure_export.h"
#include <fstream>
#include <cstdio>
#include <ctime>
#include <sys/time.h>
#include <sstream>


/***************************************************************************
 * bounds for the spinners
 *
 * Every +/- pair in this file used to step without a bound in at least one
 * direction, and several of the resulting values are not merely useless but
 * fatal: draw_mand computes 100/mand_trap_depth, mand_reset_mesh divides the
 * plot width by the mesh size, and movie.cc divides by the frame count and
 * loops forever at zero.  So each stepper now goes through clamp_int.
 *
 * The upper bound on a word length is 64 because that is the width of
 * Bitword's std::bitset (ifs.h:101); past that the words that name a trap or
 * a uv pair silently lose their leading letters.
 ***************************************************************************/
/* 64 because that is a STRUCTURAL limit, not a matter of taste: a word is a Bitword, which
   holds a std::bitset<64>, so beyond 64 letters the word silently wraps.
   It was briefly 32, on the theory that a ceiling nobody can reach is not a ceiling -- the
   cost of these passes doubles with every level, so a full-view limit set at depth 30 takes
   about a minute and depth 40 is hours.  That was a mistake: how much time to spend is the
   user's call, and deep zooms need deep words.  Checking that the holes of M keep
   accumulating at 1e-5 and below needs a connectedness depth well past 32, and capping it
   there simply removed the ability to look.  The protection against a runaway is Escape,
   which now aborts the parameter-space repaint AND the limit-set redraw, plus the fact that
   the picture visibly fills in as it goes -- not a low ceiling. */
#define MAX_WORD_DEPTH 64
#define MAX_POINT_DEPTH 64
#define MAX_UV_GRAPH_DEPTH 16    //2^depth labelled circles get drawn, one per ball
/* ifs_movie_from_path refuses anything over 2000 frames, so a limit of 100000 only let
   the user pick a value that could never work -- and the refusal went to stdout, where a
   GUI user never sees it.  movie_length is in SECONDS and is multiplied by the frame
   rate, so this has to leave room for that: see S_mand_path_create_movie. */
#define MAX_MOVIE_LENGTH 2000
#define MIN_PICTURE_SIZE 50
/* 10000 was too generous to be usable: 10^8 cells is two 300 MB RGB copies, an EPS of
   about 150 MB, and several minutes with no progress and no way out.  4000 is 16 Mcell,
   ~50 MB per copy, and still four times the finest on-screen mesh. */
#define MAX_PICTURE_SIZE 4000
//the initial parameter-space view is [-1,1]^2, which already contains all of |s| < 1;
//outside that the maps are not contractions, so there is nothing to draw
#define MAND_MAX_RADIUS 1.0

static int clamp_int(int v, int lo, int hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}


/***************************************************************************
 * parsing the parameter entry (W_point_param_entry)
 *
 * Accepted forms:
 *   re im       -- two numbers, whitespace separated
 *   re+imi      -- ordinary complex notation (re-imi, imi, i, -i also work)
 *   deg         -- a bare number: that many degrees on the circle
 *                  |s| = 1/sqrt(2), the convention certify_arc and funddom
 *                  use (certify_arc.cc: static const double R0 = 1/sqrt(2);
 *                  th = deg*M_PI/180; s = R0*(cos th, sin th))
 *   deg@r       -- degrees at radius r, off that circle
 ***************************************************************************/
static std::string trim_ws(const std::string& s) {
  size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b-a+1);
}

//strtod, but rejects anything that leaves unparsed (non-whitespace) junk,
//so a typo reports "cannot parse" instead of silently becoming 0
static bool parse_double_strict(const std::string& tok, double& out) {
  if (tok.empty()) return false;
  const char* cs = tok.c_str();
  char* endp = 0;
  double v = std::strtod(cs, &endp);
  if (endp == cs) return false;
  while (*endp == ' ' || *endp == '\t') ++endp;
  if (*endp != '\0') return false;
  out = v;
  return true;
}

/* Parse a pair of eventually-periodic words -- "fg(fffffggggg) gf(gggggfffff)", or the
   same thing in bits, "01(0000011111) 10(1111100000)" -- into the sign strings A and B of
   the landmark parameterisation.
 *
 * A coincidence u(0) = v(0) with u = P_u (C_u)^inf and v = P_v (C_v)^inf gives, letter by
 * letter, A_j = (eps(u_j) - eps(v_j))/2 over the prefix and B_k the same over the cycle,
 * with eps(f) = -1 and eps(g) = +1; then sigma is a root of
 * Q(z) = A(z)(1 - z^b) + z^a B(z).  So the whole job here is differencing the two words.
 *
 * The prefixes must have the same length (a) and the cycles the same length (b), or the two
 * orbits are not being compared at matching times.  The first letters must differ, since
 * u_0 = v_0 means the coincidence factors through a shorter one.  A_0 is then +-1, and the
 * convention is A_0 = -1, so u and v are swapped if needed -- which negates Q and leaves its
 * roots alone. */
/* Parse the FINITE form -- two words, no brackets, "fgff gfgg" or "0100 1011" -- into the
   f/g strings the coincidence solver wants.
 *
 * These are roots rather than landmarks: u and v are finite, and one asks for the s at which
 * u(0) = v(0) exactly.  That is a statement in the PAPER's normalization, f(z) = sz-1 and
 * g(z) = sz+1, where w(0) = sum_j eps_j s^j with j = 0 the outermost letter; so the condition
 * is a {0,+-1} polynomial in s and its roots are what the search returns.  (The program's own
 * base-1/2 normalization never enters, which is why the same kernel serves both.)
 *
 * Bits are translated to f/g here, because fd_core_from_coin takes f/g.  Equal lengths and
 * differing first letters are required, and the reason is the same as for the landmark form:
 * the words are compared letter by letter, and if they agreed at the first letter the
 * coincidence would factor through a shorter one. */
bool parse_uv_finite(const std::string& raw, std::string& u, std::string& v,
                     std::string& why) {
  std::string s = trim_ws(raw);
  if (s.empty()) { why = "nothing typed"; return false; }
  if (s.find('(') != std::string::npos || s.find(')') != std::string::npos) {
    why = "for roots give two FINITE words with no brackets, e.g. fgff gfgg"
          " (brackets mean an eventually periodic word, which is the landmark search)";
    return false;
  }
  std::string tok[2];
  { std::string cur; int nt = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
      char c = (i < s.size() ? s[i] : ' ');
      if (c == ' ' || c == '\t' || c == ',') {
        if (!cur.empty()) { if (nt < 2) tok[nt] = cur; ++nt; cur.clear(); }
      } else cur += c;
    }
    if (nt != 2) { why = "give exactly two words, u then v, e.g. fgff gfgg"; return false; } }
  for (int k = 0; k < 2; ++k)
    for (size_t i = 0; i < tok[k].size(); ++i)
      if (tok[k][i] != 'f' && tok[k][i] != 'g' && tok[k][i] != '0' && tok[k][i] != '1') {
        why = std::string("the letters must be f/g or bits 0/1, not '") + tok[k][i] + "'";
        return false;
      }
  if (tok[0].size() != tok[1].size()) {
    why = "u and v must be the same length (they are compared letter by letter)";
    return false;
  }
  if (tok[0].size() < 2) { why = "the words must be at least two letters long"; return false; }
  if (tok[0].size() > 60) { why = "the words must be at most 60 letters long"; return false; }
  u.clear(); v.clear();
  for (size_t i = 0; i < tok[0].size(); ++i) {
    u += (tok[0][i] == 'g' || tok[0][i] == '1') ? 'g' : 'f';
    v += (tok[1][i] == 'g' || tok[1][i] == '1') ? 'g' : 'f';
  }
  if (u[0] == v[0]) {
    why = "u and v must start with different letters, or the coincidence factors "
          "through a shorter one";
    return false;
  }
  return true;
}

bool parse_uv_pair(const std::string& raw, std::string& A, std::string& B,
                   std::string& why) {
  std::string s = trim_ws(raw);
  if (s.empty()) { why = "nothing typed"; return false; }
  /* two whitespace- or comma-separated words */
  std::string tok[2];
  { std::string cur;
    int nt = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
      char c = (i < s.size() ? s[i] : ' ');
      if (c == ' ' || c == '\t' || c == ',') {
        if (!cur.empty()) { if (nt < 2) tok[nt] = cur; ++nt; cur.clear(); }
      } else cur += c;
    }
    if (nt != 2) {
      why = "give exactly two words, u then v, e.g. fg(fffffggggg) gf(gggggfffff)";
      return false;
    } }

  std::string pre[2], cyc[2];
  for (int k = 0; k < 2; ++k) {
    size_t o = tok[k].find('(');
    size_t c = tok[k].find(')');
    if (o == std::string::npos || c == std::string::npos || c < o || c + 1 != tok[k].size()) {
      why = "each word needs its periodic part in brackets, as in fg(fffffggggg)";
      return false;
    }
    pre[k] = tok[k].substr(0, o);
    cyc[k] = tok[k].substr(o + 1, c - o - 1);
    for (int part = 0; part < 2; ++part) {
      const std::string& t = (part == 0 ? pre[k] : cyc[k]);
      for (size_t i = 0; i < t.size(); ++i)
        if (t[i] != 'f' && t[i] != 'g' && t[i] != '0' && t[i] != '1') {
          why = std::string("the letters must be f/g or bits 0/1, not '") + t[i] + "'";
          return false;
        }
    }
  }
  if (pre[0].size() != pre[1].size()) {
    why = "the two prefixes must be the same length (they are compared letter by letter)";
    return false;
  }
  if (cyc[0].size() != cyc[1].size()) {
    why = "the two periodic parts must be the same length";
    return false;
  }
  if (pre[0].empty()) {
    why = "the prefix cannot be empty: the parameterisation needs a >= 1";
    return false;
  }
  if (cyc[0].empty()) { why = "the periodic part cannot be empty"; return false; }
  if (pre[0].size() > 12 || cyc[0].size() > 12) {
    why = "a and b are each at most 12 (the kernel's LM_MAXN)";
    return false;
  }
  /* eps(f) = -1, eps(g) = +1; bits are 0 = f, 1 = g */
  int u0 = (pre[0][0] == 'g' || pre[0][0] == '1') ? +1 : -1;
  int v0 = (pre[1][0] == 'g' || pre[1][0] == '1') ? +1 : -1;
  if (u0 == v0) {
    why = "u and v must start with different letters, or the coincidence factors "
          "through a shorter one";
    return false;
  }
  int swap = (u0 > v0) ? 1 : 0;      //want A_0 = -1, i.e. u starting f
  const std::string& Pu = pre[swap], &Pv = pre[1-swap];
  const std::string& Cu = cyc[swap], &Cv = cyc[1-swap];
  A.clear(); B.clear();
  for (size_t j = 0; j < Pu.size(); ++j) {
    int eu = (Pu[j] == 'g' || Pu[j] == '1') ? +1 : -1;
    int ev = (Pv[j] == 'g' || Pv[j] == '1') ? +1 : -1;
    int d = (eu - ev)/2;
    A += (d < 0 ? '-' : (d > 0 ? '+' : '0'));
  }
  for (size_t k = 0; k < Cu.size(); ++k) {
    int eu = (Cu[k] == 'g' || Cu[k] == '1') ? +1 : -1;
    int ev = (Cv[k] == 'g' || Cv[k] == '1') ? +1 : -1;
    int d = (eu - ev)/2;
    B += (d < 0 ? '-' : (d > 0 ? '+' : '0'));
  }
  return true;
}

bool parse_parameter(const std::string& raw, cpx& out) {
  static const double R0 = 0.7071067811865476;  //1/sqrt(2), certify_arc.cc's R0
  std::string s = trim_ws(raw);
  if (s.empty()) return false;

  //"deg@r"
  size_t at = s.find('@');
  if (at != std::string::npos) {
    double deg, r;
    if (!parse_double_strict(trim_ws(s.substr(0,at)), deg)) return false;
    if (!parse_double_strict(trim_ws(s.substr(at+1)), r))   return false;
    if (r <= 0.0) return false;
    double th = deg*M_PI/180.0;
    out = cpx(r*std::cos(th), r*std::sin(th));
    return true;
  }

  //"re im" -- exactly two whitespace-separated tokens
  {
    std::istringstream iss(s);
    std::string a, b, extra;
    if ((iss >> a) && (iss >> b) && !(iss >> extra)) {
      double re, im;
      if (!parse_double_strict(a, re) || !parse_double_strict(b, im)) return false;
      out = cpx(re, im);
      return true;
    }
  }

  //"re+imi" / "re-imi" / "imi" / "i" / "-i" -- ordinary complex notation
  if (s[s.size()-1] == 'i' || s[s.size()-1] == 'I') {
    std::string body = s.substr(0, s.size()-1);
    //the last +/- that is not the leading sign and not an exponent sign
    //(as in "1e-5") separates the real part from the imaginary coefficient
    int split = -1;
    for (int i=(int)body.size()-1; i>0; --i) {
      char c = body[i];
      if ((c == '+' || c == '-') && body[i-1] != 'e' && body[i-1] != 'E') { split = i; break; }
    }
    double re, im;
    if (split < 0) {
      //pure imaginary: "0.3i", "-0.3i", "i", "-i"
      std::string im_s = body;
      if (im_s.empty() || im_s == "+") im_s = "1";
      else if (im_s == "-")            im_s = "-1";
      if (!parse_double_strict(im_s, im)) return false;
      re = 0.0;
    } else {
      std::string re_s = body.substr(0, split);
      std::string im_s = body.substr(split);   //keeps the sign
      if (im_s == "+") im_s = "1";
      else if (im_s == "-") im_s = "-1";
      if (!parse_double_strict(re_s, re)) return false;
      if (!parse_double_strict(im_s, im)) return false;
    }
    out = cpx(re, im);
    return true;
  }

  //bare "deg"
  {
    double deg;
    if (!parse_double_strict(s, deg)) return false;
    double th = deg*M_PI/180.0;
    out = cpx(R0*std::cos(th), R0*std::sin(th));
    return true;
  }
}





/***************************************************************************
 * widget functions
 ***************************************************************************/
bool Widget::contains_pixel(int x, int y) {
  return (ul.x <= x) && (x < ul.x + width) && (ul.y <= y) && (y < ul.y + height);
}

bool Widget::intersects_rectangle(const Point2d<int>& p, int w, int h) {
  return !(ul.x > p.x + w || ul.y > p.y + h || p.x > ul.x + width || p.y > ul.y + height);
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
  
  XFontStruct* font = ifsg->gui_font();
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
  XSetLineAttributes(ifsg->display, gc, 0, LineSolid, CapButt, JoinMiter);
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

void WidgetButton::redraw() {
  initial_draw();
}


/* The interface font, loaded once.
   Every widget used to call XLoadQueryFont(display, "fixed") for itself: 83 widgets, 83
   font loads, none ever freed, and -- worse -- none checking the result before
   dereferencing font->fid, so a machine without "fixed" would have segfaulted rather
   than complained.  Load one font, keep it, and fall back down a list.

   "fixed" is a 6x13 bitmap and looks like 1987.  Helvetica at 12 is a large visible
   improvement for no new dependency, since it ships with XQuartz and every X11 install.
   Every candidate is deliberately iso8859-1: the iso10646-1 versions of the same faces
   are 16-bit fonts, and XDrawString -- which is how every widget here draws -- would
   render nothing at all on one. */
XFontStruct* IFSGui::gui_font() {
  if (the_gui_font != NULL) return the_gui_font;
  /* THE WIDTH BUDGET IS 85 PIXELS, and it decides this list.
     A WidgetCheck is built at a fixed width of 105 with its text starting at x=20, and
     update_text never resizes, so any label wider than 85 is silently clipped.  Bold 12 is
     wanted for legibility -- medium 11 was too thin to read comfortably -- and it only fits
     because the longest labels were shortened at the same time ("Connectedness:" became
     "Connected:", "Plot uv graph" became "uv graph").  Widths of the worst survivor,
     "Contains 0:", measured:
         helvetica medium 11   63
         helvetica bold 11     67
         helvetica bold 12     77   <- used; 8px of margin
         lucida bold sans 11   85   exactly at the limit
         lucida bold sans 12   91   CLIPPED
     Check any new candidate, and any new label, against 85 before adding it. */
  static const char* candidates[] = {
    "-*-helvetica-bold-r-normal--12-*-iso8859-1",
    "-*-helvetica-bold-r-normal--11-*-iso8859-1",
    "-*-helvetica-medium-r-normal--11-*-iso8859-1",
    "fixed",
    NULL
  };
  for (int i=0; candidates[i] != NULL; ++i) {
    the_gui_font = XLoadQueryFont(display, candidates[i]);
    if (the_gui_font != NULL) return the_gui_font;
  }
  std::cerr << "schottky: could not load any font, not even \"fixed\".\n";
  exit(1);
  return NULL;
}


WidgetText::WidgetText(IFSGui* i, const std::string& t, int w, int h) {
  ifsg = i;
  text = t;
  height = h;
  click_signal = NULL;
  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);
  
  XFontStruct* font = ifsg->gui_font();
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
  
  XFontStruct* font = ifsg->gui_font();
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


//A one-line editable field: static "label" prefix, then editable "text".
//Setup is WidgetText's ctor verbatim (GC, "fixed" font, XTextExtents,
//autosize, pixmap); the rendering itself is new and lives in redraw(), which
//this ctor just calls once to fill in the pixmap (harmless to do before ul
//is assigned by pack_widget_upper_right: redraw's XCopyArea targets whatever
//ul happens to be at the time, but nothing flushes to the screen until the
//final initial_draw() pass in reset_and_pack_window overwrites it anyway).
WidgetEntry::WidgetEntry(IFSGui* i, const std::string& lab, const std::string& init,
                          int w, int h, void (IFSGui::*on_enter)(XEvent*)) {
  ifsg = i;
  label = lab;
  text = init;
  committed = init;
  height = h;
  click_signal = &IFSGui::S_entry_click;  //must be non-NULL, or hit-testing skips this widget entirely
  enter_signal = on_enter;
  focused = false;
  caret = (int)text.size();
  max_chars = 80;

  gc = XCreateGC(ifsg->display, RootWindow(ifsg->display, ifsg->screen), 0, NULL);

  XFontStruct* font = ifsg->gui_font();
  XSetFont(ifsg->display, gc, font->fid);
  XCharStruct te;
  int fdir, fdescent, fascent;
  std::string measure = label + text;
  XTextExtents(font, measure.c_str(), measure.size(), &fdir, &fascent, &fdescent, &te);

  int desired_width = te.rbearing - te.lbearing;
  int height_offset = (te.ascent - te.descent)/2;
  if (w > 0) {
    width = w;
  } else {
    width = desired_width + 10;
  }
  text_position = Point2d<int>(5, height/2 + height_offset);
  ascent = te.ascent;
  descent = te.descent;
  char_w = XTextWidth(font, "0", 1);
  label_px = XTextWidth(font, label.c_str(), (int)label.size());

  p = XCreatePixmap(ifsg->display, ifsg->main_window,
                    width, height, DefaultDepth(ifsg->display, ifsg->screen));

  redraw();
}

void WidgetEntry::initial_draw() {
  XCopyArea(ifsg->display, p, ifsg->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}

//Re-renders the whole pixmap from (label, text, caret, focused) every time --
//the WidgetCheck pattern, the only precedent this toolkit has for a widget
//whose appearance depends on mutable state rather than being fixed at
//construction.
void WidgetEntry::redraw() {
  XSetForeground(ifsg->display, gc, WhitePixel(ifsg->display, ifsg->screen));
  XFillRectangle(ifsg->display, p, gc, 0, 0, width, height);
  XSetForeground(ifsg->display, gc, BlackPixel(ifsg->display, ifsg->screen));
  XDrawRectangle(ifsg->display, p, gc, 0, 0, width-1, height-1);
  //focus is shown by a second, inset rectangle -- 2px "thicker" -- rather
  //than a color, since X pixel values here are opaque (see get_rgb_color)
  //and every other widget already relies on shape, not color, to read.
  if (focused) {
    XDrawRectangle(ifsg->display, p, gc, 1, 1, width-3, height-3);
  }
  XDrawString(ifsg->display, p, gc, text_position.x, text_position.y, label.c_str(), label.size());
  XDrawString(ifsg->display, p, gc, text_position.x + label_px, text_position.y, text.c_str(), text.size());
  if (focused) {
    int cx = text_position.x + label_px + caret*char_w;
    XDrawLine(ifsg->display, p, gc, cx, text_position.y - ascent, cx, text_position.y + descent);
  }
  initial_draw();
}

//KeyReleaseMask is not in the XSelectInput mask, so auto-repeat arrives as
//ordinary repeated KeyPress events and needs no special handling here.
void WidgetEntry::handle_key(XEvent* e) {
  char buf[16]; KeySym ks; int n = XLookupString(&e->xkey, buf, sizeof buf, &ks, NULL);
  if      (ks == XK_Return || ks == XK_KP_Enter) { committed = text;
                                                   if (enter_signal) (ifsg->*enter_signal)(e); redraw(); return; }
  else if (ks == XK_Escape)    { text = committed; caret = (int)text.size(); focused = false;
                                 ifsg->focus_widget = NULL; redraw(); return; }
  else if (ks == XK_Tab)       { ifsg->focus_next_entry(); return; }   //walks `widgets` for wants_keys()
  else if (ks == XK_BackSpace) { if (caret > 0) { text.erase(--caret, 1); } }
  else if (ks == XK_Delete)    { if (caret < (int)text.size()) text.erase(caret, 1); }
  else if (ks == XK_Left)      { if (caret > 0) --caret; }
  else if (ks == XK_Right)     { if (caret < (int)text.size()) ++caret; }
  else if (ks == XK_Home)      { caret = 0; }
  else if (ks == XK_End)       { caret = (int)text.size(); }
  else if (n == 1 && buf[0] >= 32 && buf[0] < 127 && (int)text.size() < max_chars)
                               { text.insert(text.begin()+caret, buf[0]); ++caret; }
  redraw();
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

void WidgetLeftArrow::redraw() {
  initial_draw();
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

void WidgetRightArrow::redraw() {
  initial_draw();
}














/****************************************************************************
 * signal handlers
 ****************************************************************************/
 
void IFSGui::S_switch_to_limit(XEvent* e) {
  if (e->type != ButtonPress) return;
  window_mode = LIMIT;
  reset_and_pack_window();
}

void IFSGui::S_switch_to_mandelbrot(XEvent* e) {
  if (e->type != ButtonPress) return;
  window_mode = MANDELBROT;
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
  int d = clamp_int(limit_depth+1, 1, MAX_WORD_DEPTH);
  if (d == limit_depth) return;
  limit_depth = d;
  std::stringstream T; T.str(""); T << limit_depth;
  W_limit_depth_label.update_text(T.str());
  draw_limit();
}

void IFSGui::S_limit_decrease_depth(XEvent* e) {
  if (e->type == KeyPress || e->type == MotionNotify) return;
  int d = clamp_int(limit_depth-1, 1, MAX_WORD_DEPTH);
  if (d == limit_depth) return;
  limit_depth = d;
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

/* The trap search, in one place.
   A trap at s is a pair of words (u,v) -- u beginning f, v beginning g -- whose balls
   u(B) and v(B) meet and whose displacement, renormalized by s^{-|u|}, lands strictly
   inside the trap-like set T_0.  Finding one certifies a whole disk of parameters about
   s inside M, which is the point of the whole exercise.

   Both the point panel and the limit-set overlay need this, and they must not disagree,
   so neither computes it for itself.  The box is deliberately tiny: we are asking about
   this parameter, not about a neighborhood of it. */
bool IFSGui::find_trap_words(int uv_depth,
                             std::vector<std::pair<Bitword,Bitword> >& out,
                             std::string& why) {
  if (trap_cache_valid && trap_cache_depth == uv_depth && trap_cache_z == IFS.z) {
    out = trap_cache_words;
    why = trap_cache_why;
    return trap_cache_ok;
  }
  out.resize(0);
  why.clear();
  std::vector<Ball> TLB;
  /* TLB_BUILD_DEPTH was 15, which was measurably too shallow: at 15 the trap-like set
     has only 2 balls and the first trap at 60 degrees needs words of length 20, while at
     18 it has 10 -- as many as this construction yields -- and the trap appears at length
     9.  Cost doubles per level (0.08s at 18, 1.4s at 22), so 18 is where the ball count
     saturates and the search is still interactive.  Measured over 40..90 degrees in 2.5
     degree steps: a trap at 16 of 21 angles, worst case 0.09s.  The misses are real, not
     an artefact -- 45 degrees is the twindragon coincidence core, where no finite trap
     exists and certify_arc's own sweep fails too. */
  const int TLB_BUILD_DEPTH = 18;
  const double TLB_BOX_HALFWIDTH = 1e-9;
  cpx box_ll = IFS.z - cpx(TLB_BOX_HALFWIDTH, TLB_BOX_HALFWIDTH);
  cpx box_ur = IFS.z + cpx(TLB_BOX_HALFWIDTH, TLB_BOX_HALFWIDTH);
  double TLB_Z, TLB_C;
  trap_cache_valid = true; trap_cache_z = IFS.z; trap_cache_depth = uv_depth;
  if (!IFS.TLB_for_region(TLB, box_ll, box_ur, TLB_BUILD_DEPTH, &TLB_C, &TLB_Z, 0)) {
    why = "no trap-like balls at this parameter";
    trap_cache_ok = false; trap_cache_words.resize(0); trap_cache_why = why;
    return false;
  }
  double trap_radius;
  int diff = IFS.check_TLB(TLB, &TLB_C, &TLB_Z, trap_radius, &out, uv_depth);
  /* check_TLB has one exit (TLB_C == NULL) that reports a depth without recording a
     word, so a non-negative answer does not by itself mean the list is populated. */
  if (diff < 0 || out.empty()) {
    out.resize(0);
    std::stringstream T;
    T << "no trap with words of length <= " << uv_depth;
    why = T.str();
    trap_cache_ok = false; trap_cache_words.resize(0); trap_cache_why = why;
    return false;
  }
  trap_cache_ok = true; trap_cache_words = out; trap_cache_why.clear();
  return true;
}


void IFSGui::S_limit_trap(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_trap = !limit_trap;
  W_limit_trap.checked = limit_trap;
  W_limit_trap.redraw();
  draw_limit();
}

void IFSGui::S_limit_trap_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(limit_trap_depth+1, 1, MAX_WORD_DEPTH);
  if (d == limit_trap_depth) return;
  limit_trap_depth = d;
  std::stringstream T; T << limit_trap_depth;
  W_limit_trap_depth_label.update_text(T.str());
  if (limit_trap) draw_limit();
}

void IFSGui::S_limit_trap_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(limit_trap_depth-1, 1, MAX_WORD_DEPTH);
  if (d == limit_trap_depth) return;
  limit_trap_depth = d;
  std::stringstream T; T << limit_trap_depth;
  W_limit_trap_depth_label.update_text(T.str());
  if (limit_trap) draw_limit();
}


/* Take the limit-set view to the trap.  Not a convenience: the certifying balls have
   radius about min_r*|s|^|u|, which at |u| = 9 and |s| = 0.7 is a couple of hundredths
   against a view three units wide -- a few pixels.  Without this the user has to find a
   four-pixel mark by eye before the drawing tells them anything. */
void IFSGui::S_limit_trap_zoom(XEvent* e) {
  if (e->type != ButtonPress) return;
  if (!limit_trap) {                     /* switch it on rather than do nothing */
    limit_trap = true;
    W_limit_trap.checked = true;
    W_limit_trap.redraw();
    draw_limit();
  }
  if (!limit_trap_located) {
    set_status("no trap to go to at this parameter");
    return;
  }
  /* four ball-radii of margin: enough to see both clusters and the gap they interleave
     across, without losing the surrounding limit set entirely */
  double radius = 4.0*limit_trap_radius;
  limit_ll = limit_trap_center - cpx(radius, radius);
  limit_ur = limit_trap_center + cpx(radius, radius);
  limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(W_limit_plot.width);
  draw_limit();
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
  int d = clamp_int(limit_uv_graph_depth-1, 1, MAX_UV_GRAPH_DEPTH);
  if (d == limit_uv_graph_depth) return;
  limit_uv_graph_depth = d;
  std::stringstream T; T.str(""); T << limit_uv_graph_depth;
  W_limit_uv_graph_depth_label.update_text(T.str());
  draw_limit();
}

void IFSGui::S_limit_uv_graph_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(limit_uv_graph_depth+1, 1, MAX_UV_GRAPH_DEPTH);
  if (d == limit_uv_graph_depth) return;
  limit_uv_graph_depth = d;
  std::stringstream T; T.str(""); T << limit_uv_graph_depth;
  W_limit_uv_graph_depth_label.update_text(T.str());
  draw_limit();
}

void IFSGui::S_limit_nifs(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_nifs = !limit_nifs;
  W_limit_nifs.checked = limit_nifs;
  W_limit_nifs.redraw();
  if (limit_2d || limit_gifs) {
    limit_2d = W_limit_2d.checked = false;
    //limit_gifs = W_limit_gifs.checked = false;
    W_limit_2d.redraw();
    //W_limit_gifs.redraw();
  }
  /* FRAME THE DIFFERENCE SET.  Lambda - Lambda is centred at 0 and is bigger than Lambda,
     which sits around 1/2; drawn in the limit set's own window it came out off centre and
     clipped, and since there is no way to pan, there was no way to see the rest of it.
     The interesting point is 2/s -- s is in M exactly when 2/s lies in the difference set
     -- so the frame is chosen to include it.  Toggling back restores the previous view. */
  if (limit_nifs) {
    nifs_saved_ll = limit_ll;
    nifs_saved_ur = limit_ur;
    nIFS nifs(3, IFS.z);
    nifs.centers[0] = -1; nifs.centers[1] = 0; nifs.centers[2] = 1;
    double r = nifs.minimal_initial_radius();
    if (r > 0) {
      double want = std::abs(1.0/IFS.z);
      if (want > r) r = want;              //keep 1/s in view
      r *= 1.15;
      limit_ll = cpx(-r, -r);
      limit_ur = cpx( r,  r);
      limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(W_limit_plot.width);
    }
  } else {
    limit_ll = nifs_saved_ll;
    limit_ur = nifs_saved_ur;
    limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(W_limit_plot.width);
  }
  draw_limit();
  mand_grid_connected_valid = false;
  if (mand_connected) draw_mand();
}

void IFSGui::S_limit_gifs(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_gifs = !limit_gifs;
  W_limit_gifs.checked = limit_gifs;
  W_limit_gifs.redraw();
  if (limit_2d || limit_nifs) {
    limit_2d = W_limit_2d.checked = false;
    limit_nifs = W_limit_nifs.checked = false;
    W_limit_2d.redraw();
    W_limit_nifs.redraw();
  }
  draw_limit();
  mand_grid_connected_valid = false;
  if (mand_connected) draw_mand();
}


void IFSGui::S_limit_2d(XEvent* e) {
  if (e->type != ButtonPress) return;
  limit_2d = !limit_2d;
  W_limit_2d.checked = limit_2d;
  W_limit_2d.redraw();
  if (limit_nifs || limit_gifs) {
    limit_nifs = W_limit_nifs.checked = false;
    //limit_gifs = W_limit_gifs.checked = false;
    W_limit_nifs.redraw();
    //W_limit_gifs.redraw();
    mand_grid_connected_valid = false;
    if (mand_connected) draw_mand();
  }
  draw_limit();
}

//mandelbrot
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
    /* With landmarks shown, a click within a few pixels of one SNAPS to it, rather
       than to wherever the pointer happened to land.  That matters because a limit
       trap only exists at an exact renormalization point: landing near sigma is not
       good enough, and a hand-placed click never hits it.  Only on a real
       ButtonPress -- a drag must still sweep the parameter freely. */
    /* mand_roots as well as mand_landmarks: the two layers share the marked-point list, so
       whichever of them is showing, a click near a mark should snap to it.  Gating this on
       the landmarks alone meant the roots were drawn but not selectable -- and a mark you
       cannot click is not much use, since selecting is how the annulus button is aimed. */
    if ((mand_landmarks || mand_roots) && e->type == ButtonPress) {
      int k = mand_landmark_near(widget_x, widget_y, 6);
      if (k >= 0) {
        /* A second click on the mark already selected means "zoom to it".  The
           asymptotic self-similarity only becomes visible once you are deep in, and
           clicking twice is a lot less work than repeated Recenter and Zoom. */
        bool again = (mand_landmark_selected == k);
        mand_landmark_selected = k;
        const fd_landmark& L = mand_landmark_list[k];
        cpx sig(L.sigma_re, L.sigma_im);
        std::stringstream T;
        T.precision(15);
        /* a root's spec begins "coin:", a landmark's "lm:"; they are different objects and
           the status line should not call one by the other's name */
        T << (L.spec[0] == 'c' ? "root " : "landmark ") << L.spec
          << "  sigma = " << L.sigma_re << " " << L.sigma_im
          << "  a=" << L.a << " b=" << L.b << " deg=" << L.deg;
        set_status(T.str());
        change_highlighted_ifs(sig);
        if (again) { IFS.set_params(sig, sig); mand_recenter(); mand_zoom(0.25); }
        else draw_mand();
        return;
      }
    }
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
    T.precision(15);
    T << "Re: " << std::real(c);
    W_mand_mouse_X.update_text(T.str());
    T.str("");
    T << "Im: " << std::imag(c);
    W_mand_mouse_Y.update_text(T.str());
  }
  
}


void IFSGui::S_mand_recenter(XEvent* e) {
  if (e->type == ButtonPress) mand_recenter();
}

void IFSGui::S_mand_zoom_in(XEvent* e) {
  if (e->type == ButtonPress) mand_zoom(0.5);
}

void IFSGui::S_mand_zoom_out(XEvent* e) {
  if (e->type != ButtonPress) return;
  /* There is nothing to see outside |s| = 1: the maps are not contractions there, so no
     attractor is defined and the picture is blank by construction.  The initial view is
     already the whole meaningful region, so stop rather than pan out into white space. */
  double radius = (mand_ur.real() - mand_ll.real())/2.0;
  if (radius >= MAND_MAX_RADIUS - 1e-12) {
    set_status("already showing all of |s| < 1; there is nothing outside it");
    return;
  }
  mand_zoom(2.0);
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
    //a group wider than an eighth of the plot leaves fewer than 8 groups across,
    //and doubling without a bound eventually makes mand_num_pixel_groups zero,
    //at which point the grid is empty and "Find boundary path" indexes into it
    int g = clamp_int(mand_pixel_group_size*2, 1, W_mand_plot.width/8);
    if (g == mand_pixel_group_size) return;
    mand_pixel_group_size = g;
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
    int d = clamp_int(mand_connected_depth+1, 1, MAX_WORD_DEPTH);
    if (d == mand_connected_depth) return;
    mand_connected_depth = d;
    std::stringstream T; T.str(""); T << mand_connected_depth;
    W_mand_connected_depth_label.update_text(T.str());
    mand_grid_connected_valid = false;
    if (mand_connected) draw_mand();
  }
}

void IFSGui::S_mand_connected_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    int d = clamp_int(mand_connected_depth-1, 1, MAX_WORD_DEPTH);
    if (d == mand_connected_depth) return;
    mand_connected_depth = d;
    std::stringstream T; T.str(""); T << mand_connected_depth;
    W_mand_connected_depth_label.update_text(T.str());
    mand_grid_connected_valid = false;
    if (mand_connected) draw_mand();
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
    int d = clamp_int(mand_contains_half_depth+1, 1, MAX_WORD_DEPTH);
    if (d == mand_contains_half_depth) return;
    mand_contains_half_depth = d;
    std::stringstream T; T.str(""); T << mand_contains_half_depth;
    W_mand_contains_half_depth_label.update_text(T.str());
    mand_grid_contains_half_valid = false;
    if (mand_contains_half) draw_mand();
  }
}

void IFSGui::S_mand_contains_half_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    int d = clamp_int(mand_contains_half_depth-1, 1, MAX_WORD_DEPTH);
    if (d == mand_contains_half_depth) return;
    mand_contains_half_depth = d;
    std::stringstream T; T.str(""); T << mand_contains_half_depth;
    W_mand_contains_half_depth_label.update_text(T.str());
    mand_grid_contains_half_valid = false;
    if (mand_contains_half) draw_mand();
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
    int d = clamp_int(mand_trap_depth+1, 1, MAX_WORD_DEPTH);
    if (d == mand_trap_depth) return;
    mand_trap_depth = d;
    std::stringstream T; T.str(""); T << mand_trap_depth;
    W_mand_trap_depth_label.update_text(T.str());
    mand_grid_trap_valid = false;
    if (mand_trap) draw_mand();
  }
}
void IFSGui::S_mand_trap_decrease_depth(XEvent* e) {
  if (e->type == ButtonPress) {
    //draw_mand divides 100 by this, so zero is a SIGFPE and not just a bad picture
    int d = clamp_int(mand_trap_depth-1, 1, MAX_WORD_DEPTH);
    if (d == mand_trap_depth) return;
    mand_trap_depth = d;
    std::stringstream T; T.str(""); T << mand_trap_depth;
    W_mand_trap_depth_label.update_text(T.str());
    mand_grid_trap_valid = false;
    if (mand_trap) draw_mand();
  }
}

/* Landmark points -- the renormalization points where CKW Lemma 9.2.5 applies, so
   the parameters at which a limit-trap search is meaningful at all.
 *
 * TWO WAYS TO GET THEM, and the choice is forced by the window.
 *
 * Enumerating every landmark with a+b <= N costs 3^N polynomials, and the count of points
 * grows by about 4.3 per level: 2421 at N = 7, 46201 at N = 9, past 100000 at N = 10.  So
 * the exhaustive list stops being useful at exactly the point where it starts to matter.
 * And it is the wrong thing to want at a deep zoom: a hole spiral accumulates on a
 * landmark of HIGH complexity, and one wants that one point, not a hundred thousand
 * others.  fd_landmarks_near searches instead of enumerating -- branch and bound on the
 * coefficients of Q, pruned by the geometric tail bound (see funddom.c) -- so a small
 * window is cheap right up to the a+b ceiling.
 *
 * So: whenever the window is small enough for the pruning bound to apply, use the targeted
 * search over the window's circumscribed disc and allow the full complexity range; when it
 * is not (a wide view, where the bound degenerates), fall back to the exhaustive list and
 * hold N down to 9, which is the most it can deliver without truncating.
 *
 * The cache therefore has to key on the window as well as on N, and the SELECTION has to
 * survive a rebuild -- it is an index, so it is re-found by its sigma; otherwise zooming
 * in on a chosen landmark would silently unselect it and the annulus button would then
 * claim nothing was selected. */
void IFSGui::mand_rebuild_landmarks() {
  /* a list that came from a typed u,v pair is not a function of the window, so nothing
     about panning or zooming should discard it */
  if (mand_landmark_list_from_uv) return;
  double half_w = 0.5*(mand_ur.real() - mand_ll.real());
  double half_h = 0.5*(mand_ur.imag() - mand_ll.imag());
  cpx ctr = 0.5*(mand_ll + mand_ur);
  double rad = std::sqrt(half_w*half_w + half_h*half_h);   //circumscribes the window
  //the pruning bound needs |sigma| + radius < 1, and there is no point searching a
  //window that covers most of the disc anyway
  bool targeted = (std::abs(ctr) + rad < 0.98) && (rad < 0.25);

  /* THE TWO HALVES ARE CACHED SEPARATELY, each against its own key.
   *
   * They used to share one, and the key recorded whether the ROOTS were included but not
   * whether the LANDMARKS were -- so this sequence left the display empty and believed it was
   * up to date: Landmarks on (list = landmarks), Landmarks off (no rebuild, list kept),
   * Roots on (rebuild -- landmarks excluded, since the layer is off), Roots off (rebuild --
   * now nothing is on, so the list is emptied), Landmarks on (every field of the key matches,
   * so the empty list is served as current, and no landmarks appear).  Keying each half on
   * its own "is this layer on" flag is what fixes that.
   *
   * Separate caches also mean toggling one layer does not re-run the other's search, which
   * matters: the root enumeration is 16 seconds at degree 10. */
  bool lm_stale = (mand_landmark_list_on != mand_landmarks ||
                   mand_landmark_list_N  != mand_landmarks_N ||
                   mand_landmark_list_targeted != targeted ||
                   (targeted && (mand_landmark_list_ll != mand_ll ||
                                 mand_landmark_list_ur != mand_ur)));
  bool rt_stale = (mand_root_list_on != mand_roots ||
                   mand_root_list_deg != mand_roots_deg ||
                   mand_root_list_targeted != targeted ||
                   (targeted && (mand_root_list_ll != mand_ll ||
                                 mand_root_list_ur != mand_ur)));
  if (!lm_stale && !rt_stale) return;

  /* ~3.5 seconds' worth at the measured 28000 leaves a second.  These run on every zoom and
     nothing can interrupt them, so they must be bounded rather than merely usually quick. */
  const long LM_LEAF_BUDGET = 100000;
  int N = mand_landmarks_N;
  if (!targeted && N > 9) N = 9;      //more than the exhaustive list can deliver

  if (lm_stale) {
    mand_lm_cache.clear();
    mand_lm_truncated = false;
    if (mand_landmarks) {
      { std::stringstream T;
        T << "landmarks: " << (targeted ? "searching for" : "enumerating")
          << " renormalization points with a+b <= " << N;
        if (targeted) T << " in this window";
        else if (N >= 9) T << " -- this takes a few seconds";
        set_status(T.str()); }
      //fd_landmarks returning exactly the cap means it truncated, so size generously
      //rather than silently showing a partial set
      std::vector<fd_landmark> buf(60000);
      int trunc = 0;
      int n = targeted
            ? fd_landmarks_near(ctr.real(), ctr.imag(), rad, N, &buf[0], (int)buf.size(),
                                LM_LEAF_BUDGET, &trunc)
            : fd_landmarks(N, &buf[0], (int)buf.size());
      /* the targeted bound can still decline (a window straddling |s| = 1, say); fall back
         rather than showing an empty set as though there were none */
      if (n < 0 && targeted) { n = fd_landmarks(N > 9 ? 9 : N, &buf[0], (int)buf.size());
                               trunc = 0; }
      if (n > 0) { buf.resize(n); mand_lm_cache = buf; }
      mand_lm_truncated = (trunc != 0) || (n == 60000);
    }
    mand_landmark_list_on = mand_landmarks;
    mand_landmark_list_N = mand_landmarks_N;
    mand_landmark_list_targeted = targeted;
    mand_landmark_list_ll = mand_ll;
    mand_landmark_list_ur = mand_ur;
  }

  if (rt_stale) {
    mand_rt_cache.clear();
    mand_rt_truncated = false;
    if (mand_roots) {
      { std::stringstream R;
        R << "roots: " << (targeted ? "searching for" : "enumerating")
          << " finite coincidences of degree <= " << mand_roots_deg;
        if (targeted) R << " in this window";
        else if (mand_roots_deg >= 10)
          R << " over the whole plane -- this takes "
            << (mand_roots_deg >= 12 ? "40" : "16") << " seconds or so; zoom in and it is quick";
        set_status(R.str()); }
      std::vector<fd_landmark> rbuf(200000);
      int trunc = 0;
      /* Targeted in a zoomed window, exhaustive at a wide view -- the same choice, for the
         same reason, as the landmarks: only words with a long shared prefix can have
         u(0) = v(0) near a given s, so a small disc prunes almost everything. */
      int nr = targeted
             ? fd_roots_near(ctr.real(), ctr.imag(), rad, mand_roots_deg,
                             &rbuf[0], (int)rbuf.size(), LM_LEAF_BUDGET, &trunc)
             : fd_roots(mand_roots_deg, &rbuf[0], (int)rbuf.size());
      if (nr > 0) { rbuf.resize(nr); mand_rt_cache = rbuf; }
      mand_rt_truncated = (trunc != 0) || (nr == 200000);
    }
    mand_root_list_on = mand_roots;
    mand_root_list_deg = mand_roots_deg;
    mand_root_list_targeted = targeted;
    mand_root_list_ll = mand_ll;
    mand_root_list_ur = mand_ur;
  }

  /* remember the selected point so it can be re-found after the rebuild: the selection is an
     index, and either half changing length moves it */
  bool had_sel = (mand_landmark_selected >= 0 &&
                  mand_landmark_selected < (int)mand_landmark_list.size());
  double sel_re = 0.0, sel_im = 0.0;
  if (had_sel) { sel_re = mand_landmark_list[mand_landmark_selected].sigma_re;
                 sel_im = mand_landmark_list[mand_landmark_selected].sigma_im; }

  mand_landmark_list = mand_lm_cache;
  mand_landmark_list.insert(mand_landmark_list.end(),
                            mand_rt_cache.begin(), mand_rt_cache.end());
  mand_landmark_selected = -1;
  if (had_sel)
    for (int i = 0; i < (int)mand_landmark_list.size(); ++i)
      if (std::abs(mand_landmark_list[i].sigma_re - sel_re) < 1e-12 &&
          std::abs(mand_landmark_list[i].sigma_im - sel_im) < 1e-12) {
        mand_landmark_selected = i; break;
      }

  int n = (int)mand_lm_cache.size(), nroots = (int)mand_rt_cache.size();
  std::stringstream T;
  if (mand_roots) {
    T << "roots: " << nroots << " finite coincidence" << (nroots == 1 ? "" : "s")
      << " of degree <= " << mand_roots_deg;
    if (targeted) T << " within " << rad << " of this window's centre";
    if (mand_rt_truncated) T << " (list incomplete: zoom in further or lower the degree)";
    T << ".  ";
  }
  if (mand_landmarks) {
    T << "landmarks: " << n << " renormalization point" << (n == 1 ? "" : "s")
      << " with a+b <= " << N;
    if (targeted) {
      T << " within " << rad << " of this window's centre";
      /* the highest complexity present is the useful number here: a hole spiral
         accumulates on a landmark of high a+b, so if the largest is small the search
         has not reached deep enough to see it */
      int best = 0;
      for (int i = 0; i < n; ++i)
        if (mand_lm_cache[i].a + mand_lm_cache[i].b > best)
          best = mand_lm_cache[i].a + mand_lm_cache[i].b;
      if (n > 0) T << "; deepest is a+b = " << best;
      if (mand_lm_truncated)
        T << " -- SEARCH CUT SHORT at this radius, so the list is incomplete: zoom in"
             " further, or lower the complexity bound";
    } else if (mand_lm_truncated) T << " (list truncated)";
    if (!targeted && mand_landmarks_N > 9)
      T << "  [a+b > 9 needs a zoomed-in window: the exhaustive list cannot go further]";
  }
  if (!mand_landmarks && !mand_roots) T << "landmarks and roots both off: nothing marked";
  set_status(T.str());
}

/* |s| = 1/2 and |s| = 1/sqrt2, on demand.  Drawn in the overlay pass so they sit on
   top of the raster and survive a repaint of the cached grid. */
/* mand_get_color applies the layers in a fixed PRIORITY order and returns the first
   that has a value, so with two layers on a given color is ambiguous unless you know
   which won.  Nothing on screen said so; this does. */
void IFSGui::mand_update_legend() {
  /* The column is only ~190px, and update_text does not resize, so this has to stay
     under about 27 characters; the full sentence is in the hover help. */
  std::string t = "col: ";
  bool any = false;
  if (mand_trap)          { t += "trap";                 any = true; }
  if (mand_contains_half) { t += (any ? ">" : ""); t += "half"; any = true; }
  if (mand_connected)     { t += (any ? ">" : ""); t += "conn"; any = true; }
  if (!any) t += "none on";
  W_mand_legend.update_text(t);
}

void IFSGui::mand_draw_guide_circles() {
  if (window_mode == LIMIT) return;
  if (!mand_circle_half && !mand_circle_sqrt2) return;
  Widget& MW = W_mand_plot;
  for (int k = 0; k < 2; ++k) {
    if (k == 0 && !mand_circle_half) continue;
    if (k == 1 && !mand_circle_sqrt2) continue;
    double r = (k == 0 ? 0.5 : 1.0/sqrt(2.0));
    //both red: a blue circle vanishes against the dark blue of the connected region
    XSetForeground(display, MW.gc, get_rgb_color(0.90, 0.10, 0.10));
    Point2d<int> c0 = mand_cpx_to_pixel(cpx(0.0, 0.0));
    int rp = (int)(r/mand_pixel_width);
    XDrawArc(display, MW.p, MW.gc, c0.x-rp, c0.y-rp, 2*rp, 2*rp, 0, 360*64);
  }
}

/* A scale bar, and the coordinate axes when they are in view.
   Parameter space is self-similar under the renormalization at a landmark, so a deep zoom
   looks very much like a shallow one and there is nothing on screen to say which it is.
   The bar is given a round length -- 1, 5, or 2 times a power of ten -- chosen so that it
   comes out between about a sixth and a third of the pane, and labelled with that length,
   so the label is always exact rather than the bar being a fixed size with an awkward
   number attached. */
void IFSGui::mand_draw_scale() {
  if (!mand_scale || window_mode == LIMIT) return;
  Widget& MW = W_mand_plot;

  /* the axes, in gray, drawn first so the bar and the picture sit on top */
  int gray = get_rgb_color(0.55, 0.55, 0.55);
  XSetForeground(display, MW.gc, gray);
  Point2d<int> origin = mand_cpx_to_pixel(cpx(0.0, 0.0));
  if (origin.y >= 0 && origin.y < MW.height)
    XDrawLine(display, MW.p, MW.gc, 0, origin.y, MW.width-1, origin.y);
  if (origin.x >= 0 && origin.x < MW.width)
    XDrawLine(display, MW.p, MW.gc, origin.x, 0, origin.x, MW.height-1);

  /* choose a round length whose bar is a reasonable fraction of the pane */
  double span = mand_pixel_width * MW.width;
  double target = span/4.0;
  double p10 = pow(10.0, floor(log10(target)));
  double mant = target/p10;
  double nice = (mant >= 5.0 ? 5.0 : (mant >= 2.0 ? 2.0 : 1.0)) * p10;
  int barpx = (int)(nice/mand_pixel_width);
  if (barpx < 8 || barpx > MW.width - 40) return;    /* degenerate: say nothing */

  std::stringstream T;
  T.precision(10);
  T << nice;
  std::string lab = T.str();

  int x0 = 12, y0 = MW.height - 14;
  int black = BlackPixel(display, screen);
  int white = WhitePixel(display, screen);
  /* a white pad under the bar and its label, so both stay readable over a dark region */
  XFontStruct* font = gui_font();
  int tw = XTextWidth(font, lab.c_str(), lab.size());
  int padw = (barpx > tw ? barpx : tw) + 8;
  XSetForeground(display, MW.gc, white);
  XFillRectangle(display, MW.p, MW.gc, x0-4, y0-font->ascent-6, padw, font->ascent+14);
  XSetForeground(display, MW.gc, black);
  XDrawLine(display, MW.p, MW.gc, x0, y0, x0+barpx, y0);
  XDrawLine(display, MW.p, MW.gc, x0, y0-3, x0, y0+3);
  XDrawLine(display, MW.p, MW.gc, x0+barpx, y0-3, x0+barpx, y0+3);
  XDrawString(display, MW.p, MW.gc, x0, y0-5, lab.c_str(), lab.size());
}

void IFSGui::S_mand_scale(XEvent* e) {
  if (e->type != ButtonPress) return;
  mand_scale = !mand_scale;
  W_mand_scale_check.checked = mand_scale;
  W_mand_scale_check.redraw();
  draw_mand();
}

void IFSGui::S_mand_circle_half(XEvent* e) {
  if (e->type != ButtonPress) return;
  mand_circle_half = !mand_circle_half;
  W_mand_circle_half_check.checked = mand_circle_half;
  W_mand_circle_half_check.redraw();
  draw_mand();
}

void IFSGui::S_mand_circle_sqrt2(XEvent* e) {
  if (e->type != ButtonPress) return;
  mand_circle_sqrt2 = !mand_circle_sqrt2;
  W_mand_circle_sqrt2_check.checked = mand_circle_sqrt2;
  W_mand_circle_sqrt2_check.redraw();
  draw_mand();
}

void IFSGui::mand_draw_landmarks() {
  mand_drawn_lm = mand_drawn_rt = 0;   //so an early return cannot leave a stale count
  if ((!mand_landmarks && !mand_roots) || window_mode == LIMIT) return;
  Widget& MW = W_mand_plot;
  int col = get_rgb_color(0.05, 0.85, 0.1);       //green: distinct from every layer
  int sel = get_rgb_color(1.0, 0.1, 0.9);
  int rootcol = get_rgb_color(1.0, 0.85, 0.0);    //amber: the finite-coincidence roots

  /* A DARK HALO UNDER EVERY DOT, and the whole lot BATCHED.
   *
   * The halo is what makes a mark legible whatever layer is underneath: green reads well on
   * the connectedness layer, which is dark, but not on "contains 0", which is a light mauve of
   * much the same brightness, and a mark that vanishes when you change layer is no use.  One
   * black pixel of surround lets the dot supply its own contrast.
   *
   * The batching matters because a root set is not small: degree <= 10 is 95776 points, and
   * one XFillArc each -- two, with the halo -- would be 190000 X requests per repaint.
   * XFillArcs takes an array, so the whole field costs three requests: haloes, landmarks,
   * roots.  Same trick as draw_limit. */
  std::vector<XArc> halo, dots[2];     //[0] landmarks, [1] roots
  int selected_x = 0, selected_y = 0;
  bool have_sel = false;
  for (int i = 0; i < (int)mand_landmark_list.size(); ++i) {
    cpx c(mand_landmark_list[i].sigma_re, mand_landmark_list[i].sigma_im);
    if (c.real() < mand_ll.real() || c.real() > mand_ur.real() ||
        c.imag() < mand_ll.imag() || c.imag() > mand_ur.imag()) continue;
    Point2d<int> p = mand_cpx_to_pixel(c);
    bool is_sel = (i == mand_landmark_selected);
    if (is_sel) { selected_x = p.x; selected_y = p.y; have_sel = true; }
    int r = is_sel ? 4 : 1;
    XArc a;
    a.x = (short)(p.x-r-1); a.y = (short)(p.y-r-1);
    a.width = (unsigned short)(2*r+2); a.height = (unsigned short)(2*r+2);
    a.angle1 = 0; a.angle2 = 23040;
    halo.push_back(a);
    a.x = (short)(p.x-r); a.y = (short)(p.y-r);
    a.width = (unsigned short)(2*r); a.height = (unsigned short)(2*r);
    /* a root's spec begins "coin:", a landmark's "lm:" -- the two kinds are different
       objects and worth telling apart at a glance */
    dots[mand_landmark_list[i].spec[0] == 'c' ? 1 : 0].push_back(a);
  }
  mand_drawn_lm = (int)dots[0].size();      //for SCHOTTKY_DUMP_STATE
  mand_drawn_rt = (int)dots[1].size();
  if (!halo.empty()) {
    XSetForeground(display, MW.gc, BlackPixel(display, screen));
    XFillArcs(display, MW.p, MW.gc, &halo[0], (int)halo.size());
  }
  for (int k = 0; k < 2; ++k)
    if (!dots[k].empty()) {
      XSetForeground(display, MW.gc, k ? rootcol : col);
      XFillArcs(display, MW.p, MW.gc, &dots[k][0], (int)dots[k].size());
    }
  if (have_sel) {   //the selection last, and ringed, so it stays visible in a dense field
    XSetForeground(display, MW.gc, BlackPixel(display, screen));
    XFillArc(display, MW.p, MW.gc, selected_x-5, selected_y-5, 10, 10, 0, 23040);
    XDrawArc(display, MW.p, MW.gc, selected_x-8, selected_y-8, 16, 16, 0, 23040);
    XSetForeground(display, MW.gc, sel);
    XFillArc(display, MW.p, MW.gc, selected_x-4, selected_y-4, 8, 8, 0, 23040);
    XDrawArc(display, MW.p, MW.gc, selected_x-7, selected_y-7, 14, 14, 0, 23040);
  }
}

/* index of the landmark nearest to a click in W_mand_plot coordinates, or -1 */
int IFSGui::mand_landmark_near(int wx, int wy, int tol_px) {
  int best = -1;
  double best_d2 = (double)tol_px*tol_px;
  for (int i = 0; i < (int)mand_landmark_list.size(); ++i) {
    cpx c(mand_landmark_list[i].sigma_re, mand_landmark_list[i].sigma_im);
    Point2d<int> p = mand_cpx_to_pixel(c);
    double dx = p.x - wx, dy = p.y - wy;
    double d2 = dx*dx + dy*dy;
    if (d2 <= best_d2) { best_d2 = d2; best = i; }
  }
  return best;
}

void IFSGui::S_mand_roots(XEvent* e) {
  if (e->type != ButtonPress) return;
  mand_roots = !mand_roots;
  W_mand_roots_check.checked = mand_roots;
  W_mand_roots_check.redraw();
  if (mand_landmark_list_from_uv) {   //asking for the enumerated sets back
    mand_landmark_list_from_uv = false;
    mand_landmark_list_N = -1;
  }
  mand_rebuild_landmarks();           //the keys decide what actually gets rebuilt
  draw_mand();
}

void IFSGui::S_mand_roots_increase(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(mand_roots_deg + 1, 1, 20);
  if (d == mand_roots_deg) return;
  mand_roots_deg = d;
  std::stringstream T; T << mand_roots_deg;
  W_mand_roots_label.update_text(T.str());
  if (mand_landmark_list_from_uv) { mand_landmark_list_from_uv = false;
                                    mand_landmark_list_N = -1; }
  if (mand_roots) { mand_rebuild_landmarks(); draw_mand(); }
}

void IFSGui::S_mand_roots_decrease(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(mand_roots_deg - 1, 1, 20);
  if (d == mand_roots_deg) return;
  mand_roots_deg = d;
  std::stringstream T; T << mand_roots_deg;
  W_mand_roots_label.update_text(T.str());
  if (mand_landmark_list_from_uv) { mand_landmark_list_from_uv = false;
                                    mand_landmark_list_N = -1; }
  if (mand_roots) { mand_rebuild_landmarks(); draw_mand(); }
}

void IFSGui::S_mand_landmarks(XEvent* e) {
  if (e->type != ButtonPress) return;
  mand_landmarks = !mand_landmarks;
  W_mand_landmarks_check.checked = mand_landmarks;
  W_mand_landmarks_check.redraw();
  /* toggling is how the user asks for the enumerated landmarks back after a u,v search */
  if (mand_landmark_list_from_uv) {
    mand_landmark_list_from_uv = false;
    mand_landmark_list_N = -1;
  }
  /* unconditionally, not just when turning ON: turning the layer off has to drop its marks
     from the shared list, or they stay on screen while the other layer is showing */
  mand_rebuild_landmarks();
  draw_mand();
}

void IFSGui::S_mand_landmarks_decrease(XEvent* e) {
  if (e->type != ButtonPress) return;
  int N = clamp_int(mand_landmarks_N - 1, 3, 12);
  if (N == mand_landmarks_N) return;
  mand_landmarks_N = N;
  if (mand_landmark_list_from_uv) { mand_landmark_list_from_uv = false;
                                    mand_landmark_list_N = -1; }
  std::stringstream T; T << mand_landmarks_N;
  W_mand_landmarks_label.update_text(T.str());
  if (mand_landmarks) { mand_rebuild_landmarks(); draw_mand(); }
}

void IFSGui::S_mand_landmarks_increase(XEvent* e) {
  if (e->type != ButtonPress) return;
  /* 12 is the kernel's ceiling (LM_MAXN).  Above 9 the exhaustive list cannot deliver
     -- it passes 100000 points and truncates -- so a+b > 9 is only reachable with a
     zoomed-in window, where mand_rebuild_landmarks switches to the targeted search.
     7 (2421 points, 0.12s) is still the default, because at a wide view anything more
     reads as noise rather than as a set of marks. */
  int N = clamp_int(mand_landmarks_N + 1, 3, 12);
  if (N == mand_landmarks_N) return;
  mand_landmarks_N = N;
  if (mand_landmark_list_from_uv) { mand_landmark_list_from_uv = false;
                                    mand_landmark_list_N = -1; }
  std::stringstream T; T << mand_landmarks_N;
  W_mand_landmarks_label.update_text(T.str());
  if (mand_landmarks) { mand_rebuild_landmarks(); draw_mand(); }
}

/* The limit traps of a whole fundamental annulus at the selected landmark point.
 *
 * At a renormalization point sigma the set of limit-trap parameters C is invariant
 * under C -> sigma^b C, so the annulus rho|sigma|^b <= |C| <= rho is a fundamental
 * domain for that action on E_sigma = the quotient of C^* by <sigma^b>.  Covering
 * one such annulus therefore settles EVERY C, which is what puts sigma in the
 * interior of M.  This is the picture of Figure 25 of the paper: the disc of C
 * colored by the least tail length c at which each is certified, with the two
 * bounding circles of one annulus drawn in red.
 *
 * The mathematics is funddom's, reached through funddom_core.h, so there is only one
 * implementation of it.  The two things this has to get right locally are the
 * trap-like balls -- which must be computed at sigma, and there must be MANY of
 * them, since the ten that trap_like_balls_from_balls yields cover essentially
 * nothing -- and rho, which has to sit inside the range where Y(C) can still land
 * in the difference set at all. */
void IFSGui::S_mand_annulus(XEvent* e) {
  if (e->type != ButtonPress) return;
  if (mand_landmark_selected < 0 ||
      mand_landmark_selected >= (int)mand_landmark_list.size()) {
    set_status("no landmark selected: switch Landmarks on and click one first");
    return;
  }
  const fd_landmark& L = mand_landmark_list[mand_landmark_selected];
  cpx sigma(L.sigma_re, L.sigma_im);

  /* A ROOT IS NOT A LANDMARK, and this figure is only about landmarks.
   *
   * The selection can now be either kind: "lm:<A>:<B>:<k>" is a landmark, "coin:<u>:<v>:<k>"
   * a root of a finite coincidence.  This picture is the set of C admitting a LIMIT trap for
   * the asymptotic family sigma + C sigma^{bn}, which is what a landmark's infinite
   * coincidence provides.  At a finite coincidence u(0) = v(0) the agreement is already exact
   * and there is no asymptotic renormalization to take a limit of, so there is no annulus of
   * E_sigma to cover -- the formulae would still evaluate (Delta = 0 there, and rho_max comes
   * out positive), which is exactly why this has to refuse explicitly rather than quietly
   * produce a picture of nothing. */
  std::string spec(L.spec);
  if (spec.compare(0, 5, "coin:") == 0) {
    set_status("the annulus figure is for LANDMARK points, not roots: " + spec + " is a finite"
               " coincidence u(0) = v(0), which is exact rather than asymptotic, so there is no"
               " fundamental annulus of E_sigma to cover.  Select a landmark (green) instead");
    return;
  }
  std::string A, B; int pick = 0;
  { size_t p1 = spec.find(':'), p2 = spec.find(':', p1+1), p3 = spec.find(':', p2+1);
    A = spec.substr(p1+1, p2-p1-1);
    B = spec.substr(p2+1, (p3==std::string::npos ? spec.size() : p3) - p2 - 1);
    if (p3 != std::string::npos) pick = atoi(spec.c_str()+p3+1); }
  fd_core core;
  if (fd_core_from_lm(A.c_str(), B.c_str(), pick, &core) <= 0) {
    set_status("could not rebuild the renormalization data for " + spec);
    return;
  }

  set_status("annulus: finding trap-like balls at sigma ...");
  ifs F; F.set_params(sigma, sigma); F.depth = mand_annulus_ball_depth;
  std::vector<Ball> TLB;
  double d = 1e-9;
  /* ngaps, ntrials, nradial = 40, 14, 8 -- the same sweep the command line uses
     (`certify_arc dumptlbmanyz <re> <im> 20 1e-9 40 14 8`) and the README quotes.  This
     used to be 40, 8, 6, which is a materially weaker set: at the landmark
     lm:-:++----++:0 it yielded 502 balls where the CLI settings give 1087, and the
     coverage of one fundamental annulus went 90.69% with the weak set against 97.65%
     with the strong one.  Beyond about a thousand balls it saturates, so there is no
     point going further. */
  if (!F.trap_like_balls_many(TLB, sigma-cpx(d,d), sigma+cpx(d,d),
                              mand_annulus_ball_depth, 40, 14, 8, true, 0)
      || TLB.size() == 0) {
    set_status("annulus: no trap-like balls at this landmark");
    return;
  }
  std::vector<double> balls3;
  balls3.reserve(3*TLB.size());
  for (int i=0; i<(int)TLB.size(); ++i) {
    balls3.push_back(TLB[i].center.real());
    balls3.push_back(TLB[i].center.imag());
    balls3.push_back(TLB[i].radius);
  }
  if (fd_solver_init(&core, &balls3[0], (int)TLB.size()) != 0) {
    set_status("annulus: could not initialise the limit-trap solver");
    return;
  }
  fd_info info;
  fd_solver_info(&info);
  if (!(info.rho_max > 0.0)) {
    fd_solver_free(); set_status("annulus: no usable rho at this landmark"); return;
  }
  /* THE ASYMPTOTIC REGIME, and why rho cannot simply be a fraction of rho_max.
   *
   * The true statement is about sigma + C sigma^{bn} as n -> infinity, and the set
   * of C admitting a limit trap is exactly invariant under C -> sigma^b C.  Near
   * |C| = rho_max, though, the picture still remembers the finite scale: what one
   * sees there is the limit-trap structure in a small but nontrivial NEIGHBOURHOOD
   * of sigma, not the scale-invariant asymptotic object the red annulus is meant to
   * be a fundamental domain OF.  To get the asymptotic picture one has to go deep,
   * which is what make_asymp_fig.py does for the paper's figure: at s0 it uses
   * rho = 5e-3 against rho_max = 0.354, i.e. about rho_max/70.
   *
   * That is not free.  Certified at tail length c implies sigma^b C certified at
   * c + b, so each factor |sigma|^b inward costs b more levels; descending from
   * rho_max to rho costs about log(rho_max/rho)/log(1/|sigma|) extra levels.  So
   * cmax is derived from the descent rather than guessed -- a small rho with a
   * cmax chosen for rho_max would show a spuriously uncovered picture. */
  const double ASYMP_FRACTION = 1.0/70.0;
  double rho_outer = ASYMP_FRACTION*info.rho_max;
  double abs_sigma = std::abs(sigma);
  /* WHICH PERIOD TO RENDER, and why not the outermost one.
   *
   * The certification is ONE-SIDED: certified at tail length c implies sigma^b C is
   * certified at c + b, but NOT conversely.  So a period further in is strictly EASIER
   * for the test than the one outside it -- and since the set of C admitting a limit
   * trap is exactly invariant under C -> sigma^b C, covering ANY ONE period settles
   * every C.  Rendering the outermost period therefore reports the hardest case and
   * understates the answer, sometimes badly: at the landmark lm:-:++----++:0 the
   * outermost period plateaus at 97.68% however much depth it is given, while the next
   * one in reaches 100.0000% -- checked to 1080x1080, i.e. 1.17M pixels with none
   * uncovered.  The outermost period is also the one that still remembers the finite
   * scale, since its outer edge sits at rho_max/70 rather than at zero.
   *
   * So the picture is the disc |C| <= rho_outer*|sigma|^(PERIOD*b), and the annulus
   * marked in red inside it is that period -- a fundamental domain of E_sigma either
   * way, just a scale-invariant copy further in. */
  const int PERIOD = 1;
  double rho = rho_outer*std::pow(abs_sigma, (double)(PERIOD*core.b));
  int extra = 0;
  if (abs_sigma > 0.0 && abs_sigma < 1.0)
    extra = (int)std::ceil(std::log(info.rho_max/rho_outer)/std::log(1.0/abs_sigma));
  /* ...and (PERIOD+1)*b more.  The descent term above only pays for getting from
     rho_max down to rho_outer.  One fundamental annulus spans a factor |sigma|^b in
     radius and each such factor inward costs b more levels, so reaching the INNER edge
     of period PERIOD costs (PERIOD+1)*b beyond that descent.  At b = 1 this is nearly
     invisible, which is why it went unnoticed at s0; at b = 8 it is sixteen levels, and
     without them the deep part of the picture comes out spuriously uncovered -- which is
     exactly how a picture that ought to be invariant under C -> sigma^b C ends up
     looking as though it is not. */
  const int RES = mand_annulus_res;
  int cmax = mand_annulus_cmax + extra + (PERIOD+1)*core.b;
  /* 120, not 45.  cmax has to reach base + extra + (PERIOD+1)*b, and b is not small at
     every landmark: at lm:-:+-----++++:0, where b = 10, the formula asks for 60, so a cap
     of 45 silently truncated the depth exactly where a large b makes depth matter most.
     funddom's own limit is 253 (the raster byte reserves 254 and 255); the reason not to go
     there is cost, which grows with cmax -- though mildly, since the search prunes. */
  if (cmax > 120) cmax = 120;
  figexp::Figure fig;
  fig.title = "limit traps at " + spec;
  fig.set_window(-rho, -rho, rho, rho);
  fig.raster.set_size(RES, RES);
  /* Coverage is quoted over the FUNDAMENTAL ANNULUS, not the whole disc.  The disc
     runs in to C = 0, where certification needs unboundedly deep c (each factor
     |sigma|^b inward costs b more levels), so a disc figure can never reach 100% at
     finite cmax and quoting it understates the result.  One annulus is a fundamental
     domain of E_sigma, so covering it settles every C -- that is the number the
     mathematics is about, and the one comparable with the paper. */
  double inner_r = rho*std::pow(abs_sigma, core.b);
  long inside = 0, covered = 0;          /* over the whole disc, for reference */
  long ann_inside = 0, ann_covered = 0;  /* over one fundamental annulus */
  /* For the verdict below: how much of the annulus lies in T_sigma at all (the ceiling
     any covering could reach), and how much was already certified b levels shallower.
     Both come free from the per-pixel level that is computed anyway. */
  long ann_in_T = 0, ann_cov_shallow = 0;
  /* Depth of the survival test that decides T_sigma.  This was 8, which is far too weak
     to be worth reporting: at the CKW hexahole -- a PROVED boundary point, and the one
     control where the answer is known -- depth 8 says 100% of the annulus is in T_sigma,
     while depth 26 finds 0.107% of it that provably is not.  That 0.107% is the whole
     signal (see the verdict below), so a test that cannot see it makes the ceiling
     useless.  It is also nearly free: 26 costs about 2% more than 8, because surv()
     prunes as hard as the certification does. */
  const int SURV_DEPTH = 26;
  bool aborted = false;
  for (int iy = 0; iy < RES && !aborted; ++iy) {
    if ((iy % 16) == 0) {
      std::stringstream T;
      T << "annulus at " << spec << ": row " << iy << " of " << RES
        << "  (" << TLB.size() << " balls, rho=" << rho << ", Escape to stop)";
      set_progress(double(iy)/double(RES), T.str());
      XEvent ev;
      while (XCheckMaskEvent(display, KeyPressMask, &ev))
        if (XLookupKeysym(&ev.xkey, 0) == XK_Escape) aborted = true;
    }
    for (int ix = 0; ix < RES; ++ix) {
      double x = (2.0*(ix+0.5)/RES - 1.0)*rho;
      double y = (2.0*(iy+0.5)/RES - 1.0)*rho;
      unsigned char r = 255, g = 255, b = 255;      //white: not certified
      double q2 = x*x + y*y;
      if (q2 <= rho*rho) {
        ++inside;
        bool in_annulus = (q2 >= inner_r*inner_r);
        if (in_annulus) ++ann_inside;
        int lev = fd_level(x, y, cmax);
        if (lev < 0) {
          if (!fd_survives(x, y, SURV_DEPTH)) { r = g = b = 205; } //gray: outside T_sigma
          else if (in_annulus) ++ann_in_T;                  //in T_sigma, not certified
        } else {
          ++covered;
          if (in_annulus) {
            ++ann_covered;
            ++ann_in_T;                                     //certified => in T_sigma
            if (lev <= cmax - core.b) ++ann_cov_shallow;
          }
          /* viridis by tail length: dark for shallow, bright yellow for deep.  Same
             ramp as render_funddom.py, which is why it lives in figure_export. */
          double t = (cmax > 1 ? double(lev)/double(cmax) : 0.0);
          unsigned char rgb[3];
          figexp::viridis(t, rgb);
          r = rgb[0]; g = rgb[1]; b = rgb[2];
        }
      }
      fig.raster.set_pixel(ix, iy, r, g, b);
    }
  }
  fd_solver_free();
  if (aborted) { set_status("annulus: stopped, nothing written"); return; }

  /* the two bounding circles of one fundamental annulus, as true vector paths */
  fig.add_circle(0.0, 0.0, rho,     figexp::Style::stroke(0.85, 0.05, 0.05, 1.2));
  fig.add_circle(0.0, 0.0, inner_r, figexp::Style::stroke(0.85, 0.05, 0.05, 1.2));

  static int annulus_counter = 0;
  char base[64];
  std::snprintf(base, sizeof base, "schottky_annulus_%03d", annulus_counter++);
  figexp::Options opt;
  opt.raster_px = RES;
  opt.draw_frame = true;
  std::string wrote, err, failed;
  const char* ext[3] = {".png", ".eps", ".pdf"};
  for (int k=0; k<3; ++k) {
    std::string name = std::string(base) + ext[k];
    if (figexp::write_auto(fig, opt, name, &err)) wrote += (wrote.empty()?"":" ") + name;
    else failed += (failed.empty()?"":"; ") + err;
  }
  /* THE VERDICT.
   *
   * A coverage percentage on its own says nothing: at the CKW hexahole, a proved BOUNDARY
   * point, coverage reaches 99.87% and stops, and at s0, a proved INTERIOR point, it
   * reaches 99.93% and keeps climbing.  Thresholding a number that close is hopeless, and
   * an earlier attempt to do it -- with constants fitted to those two points -- reported
   * s0 itself as "inconclusive".
   *
   * What separates them is not the coverage but the CEILING.  Two facts do the work:
   *
   *   (1) The survival test is a NECESSARY condition.  If C fails it, C admits no limit
   *       trap at all, however deep one searches.  So the fraction of the annulus in
   *       T_sigma is a hard ceiling on any covering, and a ceiling BELOW 100% means a
   *       positive-measure set of C provably has none -- an obstruction, not a shortage
   *       of depth.  Being scale-invariant, that hole recurs at every scale: exactly the
   *       boundary-point signature.  It is also independent of cmax, which is what makes
   *       it a stable thing to test: at the hexahole the ceiling reads 99.8927% at cmax
   *       26, 34 and 42 alike.
   *
   *   (2) Conversely the certified set is invariant under C -> sigma^b C, so covering ONE
   *       fundamental annulus covers every C, and sigma is interior.  Reaching 100% here
   *       is therefore a real (floating-point) certificate, not a suggestive number.
   *
   * Hence three outcomes, in this order, and each is a statement rather than a guess:
   * ceiling short of 100% -> boundary; else coverage at 100% -> interior; else the run is
   * too shallow and says so.  Measured, with PERIOD = 1 and SURV_DEPTH = 26:
   *     CKW hexahole (BOUNDARY, proved)  ceiling 99.8927%, coverage 99.8742% -> boundary?
   *     s0           (INTERIOR, proved)  ceiling 100%,     coverage 100.0000% -> interior?
   *     lm:-:++----++:0                  ceiling 100%,     coverage 100.0000% -> interior?
   * All three come out right, and the first two are the only points where the answer is
   * known independently.
   *
   * The tolerances are in PIXELS, not in percent, because that is what the raster can
   * actually resolve: a percentage threshold silently means a different number of pixels
   * at a different RES. */
  double ann_pct   = ann_inside ? 100.0*ann_covered/ann_inside : 0.0;
  double ceil_pct  = ann_inside ? 100.0*ann_in_T/ann_inside    : 0.0;
  double gain_pct  = ann_inside ? 100.0*(ann_covered - ann_cov_shallow)/ann_inside : 0.0;
  double short_pct = ceil_pct - ann_pct;
  long ann_missing = ann_inside - ann_covered;     /* uncovered pixels of the annulus */
  long ann_out_T   = ann_inside - ann_in_T;        /* pixels PROVABLY not in T_sigma */
  /* A handful of pixels is noise -- a curve clipping a pixel corner, or one parameter
     landing a rounding error away from a ball.  A hundred is a region. */
  const long PIX_NOISE = 5;      /* this many uncovered still counts as "covered" */
  const long PIX_REGION = 20;    /* this many outside T_sigma counts as a real hole */
  const char* verdict_short;
  std::string verdict;
  if (ann_out_T >= PIX_REGION) {
    verdict_short = "boundary?";
    verdict = "a positive-measure set of C (that gray region) FAILS the survival test, so"
              " it admits no limit trap at any depth; the hole is scale-invariant, so"
              " sigma most likely lies on the BOUNDARY of M -- as at the CKW hexahole,"
              " whose ceiling is 99.89%";
  } else if (ann_missing <= PIX_NOISE) {
    verdict_short = "interior?";
    verdict = "one whole fundamental annulus of E_sigma is certified, and the limit-trap"
              " set is invariant under C -> sigma^b C, so every C is covered: sigma lies"
              " in the INTERIOR of M, up to floating point";
  } else {
    verdict_short = "too shallow";
    verdict = "the annulus is not yet closed, but nothing in it is ruled out either (the"
              " ceiling is 100%), so this is a shortage of depth rather than an"
              " obstruction -- raise cmax and it should close";
  }
  { std::stringstream V;
    V.precision(3);
    V << ".  " << ann_missing << " of " << ann_inside << " pixels uncovered, of which "
      << ann_out_T << " provably admit no limit trap; the last " << core.b
      << " levels of depth bought " << gain_pct << "%";
    verdict += V.str(); }

  { std::ofstream sc((std::string(base) + ".txt").c_str());
    if (sc) {
      sc.precision(17);
      sc << "# limit traps over a fundamental annulus of E_sigma\n";
      sc << "spec        " << spec << "\n";
      sc << "sigma       " << L.sigma_re << " " << L.sigma_im << "\n";
      sc << "a b         " << core.a << " " << core.b << "\n";
      sc << "Delta       " << core.Delta_re << " " << core.Delta_im << "\n";
      sc << "Pprime      " << core.Pp_re << " " << core.Pp_im << "\n";
      sc << "period      " << PERIOD
         << "   # which fundamental annulus; deeper is strictly easier for the test,\n";
      sc << "            "
         << "    # and covering any ONE settles every C (see S_mand_annulus)\n";
      sc << "rho_outer   " << rho_outer << "   # outer edge of period 0\n";
      sc << "rho         " << rho << "   # outer edge of the period rendered\n";
      sc << "rho_max     " << info.rho_max << "\n";
      sc << "inner_rho   " << inner_r << "\n";
      sc << "cmax        " << cmax << "  (" << mand_annulus_cmax << " + " << extra
         << " for the descent to rho_outer + " << (PERIOD+1)*core.b << " = ("
         << PERIOD << "+1)*b to reach the inner edge of period " << PERIOD << ")\n";
      sc << "regime      asymptotic: rho_outer = rho_max/" << (int)(1.0/ASYMP_FRACTION)
         << ", so the picture should be invariant under C -> sigma^b C\n";
      sc << "resolution  " << RES << "\n";
      sc << "ball_depth  " << mand_annulus_ball_depth << "\n";
      sc << "balls       " << TLB.size() << "\n";
      sc << "surv_depth  " << SURV_DEPTH << "   # depth of the T_sigma survival test\n";
      sc << "R_Gamma     " << info.R_Gamma << "\n";
      sc << "annulus_certified_percent " << ann_pct << "\n";
      sc << "disc_certified_percent    " << (inside ? 100.0*covered/inside : 0.0) << "\n";
      sc << "annulus_in_T_sigma_percent " << ceil_pct
         << "   # the CEILING: nothing outside T_sigma can ever be certified\n";
      sc << "shortfall_percent          " << short_pct << "\n";
      sc << "annulus_pixels             " << ann_inside << "\n";
      sc << "uncovered_pixels           " << ann_missing << "\n";
      sc << "no_limit_trap_pixels       " << ann_out_T
         << "   # PROVABLY admit no limit trap: >= " << PIX_REGION << " means boundary\n";
      sc << "gain_from_last_" << core.b << "_levels  " << gain_pct
         << "   # is more depth still buying coverage?\n";
      sc << "reference_hexahole_ceiling 99.8927   # a PROVED boundary point: ceiling < 100%\n";
      sc << "reference_s0_ceiling       100.0     # a PROVED interior point: ceiling = 100%\n";
      sc << "reference_lm--0_ceiling    97.4612   # lm:--0:0-+-0+-:2, a vivid boundary case\n";
      sc << "verdict     " << verdict << "\n";
    } }

  std::stringstream T;
  T.precision(4);
  T << "annulus at " << spec << " period " << PERIOD << ": " << ann_pct
    << "% certified of a " << ceil_pct << "% ceiling -- " << verdict
    << "  (" << TLB.size() << " balls, cmax " << cmax << ") -- wrote "
    << base << ".{png,eps,pdf,txt}";
  if (!failed.empty()) T << "  [" << failed << "]";
  set_status(T.str());
  /* and keep a short form on screen: the status line is transient, and this is the
     number the whole exercise exists to produce */
  { std::stringstream R;
    R.precision(4);
    R << ann_pct << "% " << verdict_short;   /* stays under the 190px column */
    mand_annulus_last = R.str();
    W_mand_annulus_result.update_text(mand_annulus_last); }
}


void IFSGui::S_mand_output_window(XEvent* e) {
  if (e->type != ButtonPress) return;
  int old_prec = std::cout.precision();
  std::cout.precision(15);
  std::cout << "Lower left: " << mand_ll << ", upper right: " << mand_ur << "\n";
  std::cout.precision(old_prec);
}


void IFSGui::S_mand_output_picture(XEvent* e) {
  if (e->type != ButtonPress) return;
  /* Say so before starting.  This recomputes the whole view at the chosen resolution with
     no progress and no abort -- at the top of the range that is 16 Mcell and tens of
     seconds -- so the one thing it can do is not leave the user wondering whether the
     click registered. */
  { std::stringstream T;
    T << "writing a " << mand_output_picture_size << "x" << mand_output_picture_size
      << " picture -- recomputing the view, which cannot be stopped once started";
    set_status(T.str()); }
  std::vector<std::vector<Point3d<unsigned char> > > bmp(mand_output_picture_size);
  for (int i=0; i<mand_output_picture_size; ++i) {
    bmp[i].resize(mand_output_picture_size);
  }
  IFS.draw_mand_to_array(bmp, mand_ll, mand_ur, 
                         (mand_connected ? mand_connected_depth : 0), 
                         (mand_contains_half ? mand_contains_half_depth : 0));
  //double PI = 3.14159265358979;
  //IFS.draw_mand_to_array_radial(bmp, 
  //                              mand_output_picture_size, 
  //                              0.89*1/sqrt(2.0), 1.01*(1/sqrt(2.0)),
  //                              PI/2.0, 0.04,
  //                              (mand_connected ? mand_connected_depth : 0), 
  //                              (mand_contains_half ? mand_contains_half_depth : 0));
  /* Export through figure_export rather than write_bitmap.  Three things this buys:
     the vector formats a paper wants, overlays that stay vector (the circle |s|=1/sqrt2
     and the current parameter, which no bitmap could carry sharply), and a sidecar .txt
     recording the window and every enabled layer, which is what makes a published figure
     regenerable.  All three formats are written from one click -- there is no format
     widget yet, and offering the choice this way costs nothing and never overwrites,
     since the basename auto-increments. */
  static int export_counter = 0;
  char base[64];
  std::snprintf(base, sizeof base, "schottky_mand_%03d", export_counter++);

  figexp::Figure F;
  F.title = "schottky parameter space";
  F.set_window(mand_ll.real(), mand_ll.imag(), mand_ur.real(), mand_ur.imag());
  int np = mand_output_picture_size;
  F.raster.set_size(np, np);
  //draw_mand_to_array indexes bmp[column][row] with the row increasing UPWARD
  //(ifs_picture.cc: y = region_ll.imag() + (j+0.5)*pixel_width), which is already
  //figexp::Raster's convention, so there is no flip here
  for (int i=0; i<np; ++i)
    for (int j=0; j<np; ++j)
      F.raster.set_pixel(i, j, bmp[i][j].x, bmp[i][j].y, bmp[i][j].z);

  //the circle |s| = 1/sqrt2, on which all the certified-arc work lives
  F.add_circle(0.0, 0.0, 1.0/sqrt(2.0), figexp::Style::stroke(0.85, 0.1, 0.1, 0.7));
  //and where the highlighted parameter currently sits
  F.add_dot(IFS.z.real(), IFS.z.imag(), figexp::Style::fill(0.0, 0.35, 0.9), 2.5);

  figexp::Options opt;
  opt.raster_px = np;
  opt.draw_frame = true;
  std::string wrote, failed, err;
  const char* ext[3] = {".png", ".eps", ".pdf"};
  for (int k=0; k<3; ++k) {
    std::string name = std::string(base) + ext[k];
    if (figexp::write_auto(F, opt, name, &err)) wrote += (wrote.empty() ? "" : " ") + name;
    else { failed += (failed.empty() ? "" : "; ") + err; }
  }

  //the sidecar: everything needed to reproduce this picture
  {
    std::string sname = std::string(base) + ".txt";
    std::ofstream sc(sname.c_str());
    if (sc) {
      sc.precision(17);
      sc << "# schottky parameter-space figure\n";
      sc << "lower_left  " << mand_ll.real() << " " << mand_ll.imag() << "\n";
      sc << "upper_right " << mand_ur.real() << " " << mand_ur.imag() << "\n";
      sc << "pixels      " << np << "\n";
      sc << "parameter   " << IFS.z.real() << " " << IFS.z.imag() << "\n";
      sc << "mesh_size   " << mand_pixel_group_size << "\n";
      sc << "connected   " << (mand_connected ? mand_connected_depth : 0) << "\n";
      sc << "contains_half " << (mand_contains_half ? mand_contains_half_depth : 0) << "\n";
      sc << "trap        " << (mand_trap ? mand_trap_depth : 0) << "\n";
      wrote += " " + sname;
    }
  }

  if (failed.empty()) set_status("wrote " + wrote);
  else                set_status("export failed: " + failed);
}
  
  
void IFSGui::S_mand_output_picture_increase_size(XEvent* e) {
  if (e->type != ButtonPress) return;
  int s = clamp_int(mand_output_picture_size + 50, MIN_PICTURE_SIZE, MAX_PICTURE_SIZE);
  if (s == mand_output_picture_size) return;
  mand_output_picture_size = s;
  std::stringstream T;
  T.str(""); T << mand_output_picture_size;
  W_mand_output_picture_size_label.update_text(T.str());
}

void IFSGui::S_mand_output_picture_decrease_size(XEvent* e) {
  if (e->type != ButtonPress) return;
  //the old guard was >= 50, which let the size land on exactly 0 and then
  //asked draw_mand_to_array to fill an empty array
  int s = clamp_int(mand_output_picture_size - 50, MIN_PICTURE_SIZE, MAX_PICTURE_SIZE);
  if (s == mand_output_picture_size) return;
  mand_output_picture_size = s;
  std::stringstream T;
  T.str(""); T << mand_output_picture_size;
  W_mand_output_picture_size_label.update_text(T.str());
}
  







//point

//enter_signal for W_point_param_entry.  The widget owns only the text typed
//into it; parsing and acting on it is entirely this handler's job (see the
//"Reading the value" note on WidgetEntry).
void IFSGui::S_point_param_entered(XEvent* e) {
  cpx c;
  if (!parse_parameter(W_point_param_entry.text, c)) {
    set_status("cannot parse: " + W_point_param_entry.text);
    return;
  }
  if (window_mode == LIMIT) {
    //no mandelbrot plot exists in this mode to repaint a patch of or
    //recenter on -- just move the highlighted parameter and redraw the
    //limit set, which is the only picture on screen.
    IFS.set_params(c,c);
    draw_limit();
    recompute_point_data();
    return;
  }
  if (c.real() < mand_ll.real() || c.real() > mand_ur.real() ||
      c.imag() < mand_ll.imag() || c.imag() > mand_ur.imag()) {
    //outside the current mandelbrot window: recenter on it rather than
    //moving the highlighted parameter somewhere off screen with no feedback
    IFS.set_params(c,c);
    mand_recenter();
  } else {
    change_highlighted_ifs(c);   //repaints only a small patch (ifs_gui.cc)
  }
  set_status("moved to " + W_point_param_entry.text);
}

/* The landmark points belonging to a u,v pair the user typed.
 *
 * Given the pair, sigma is not chosen but SOLVED FOR: it is a root of
 * Q(z) = A(z)(1 - z^b) + z^a B(z).  Q has integer coefficients, so its roots come in a
 * Galois orbit and there is generally more than one -- but only those with 1/2 < |sigma| < 1
 * are of any use (outside the unit circle the maps are not contractions, and inside |s| = 1/2
 * the attractor is a Cantor set, so nothing can meet), and of each conjugate pair only the
 * one with Im sigma > 0 is new, since M is symmetric under conjugation.  Those are exactly
 * the roots fd_core_from_lm calls admissible, and it returns their count, so the whole orbit
 * comes out of the same code path funddom uses -- which means every spec reported here can be
 * handed straight back to the command line.
 *
 * The results are loaded into the landmark list, so they are drawn, clicked and selected by
 * the machinery that already exists: a click within 6 pixels snaps to one, a second click
 * zooms to it, and the annulus button then works on it unchanged.  The list stays until the
 * user asks for the enumerated landmarks again (by toggling Landmarks or moving N), which
 * is what mand_landmark_list_from_uv guards. */
void IFSGui::S_point_uv_entered(XEvent* e) {
  /* Brackets mean an eventually periodic word, i.e. the landmark search; two bare words mean
     a finite coincidence, i.e. the roots search.  Dispatching on the syntax means Enter does
     the right thing without the user having to remember which control to press, and the
     `roots` button is there for when one wants to be explicit. */
  if (W_point_uv_entry.text.find('(') == std::string::npos) { do_uv_roots_search(); return; }
  std::string A, B, why;
  if (!parse_uv_pair(W_point_uv_entry.text, A, B, why)) {
    set_status("uv pair: " + why);
    return;
  }
  fd_core core;
  int n = fd_core_from_lm(A.c_str(), B.c_str(), 0, &core);
  if (n < 0) {
    set_status("uv pair: A = " + A + ", B = " + B +
               " is not an admissible sign pair for the landmark parameterisation");
    return;
  }
  if (n == 0) {
    set_status("uv pair: A = " + A + ", B = " + B + " gives no root with 1/2 < |s| < 1 and"
               " Im s > 0, so this pair has no landmark point (the roots are all outside the"
               " annulus where an attractor can be connected)");
    return;
  }
  std::vector<fd_landmark> found;
  for (int pick = 0; pick < n; ++pick) {
    fd_core c2;
    if (fd_core_from_lm(A.c_str(), B.c_str(), pick, &c2) <= 0) continue;
    fd_landmark L;
    L.sigma_re = c2.sigma_re; L.sigma_im = c2.sigma_im;
    L.a = c2.a; L.b = c2.b; L.deg = c2.a + c2.b;
    std::snprintf(L.spec, sizeof L.spec, "lm:%s:%s:%d", A.c_str(), B.c_str(), pick);
    found.push_back(L);
  }
  if (found.empty()) { set_status("uv pair: could not rebuild the roots"); return; }

  uv_show_results(found, "uv pair -> A = " + A + ", B = " + B, "landmark point");
}

/* The tail shared by the landmark and the root searches: show the points, select one, take
   the view to it, and say what was found.  Both put their results into the landmark list, so
   the drawing, the 6px click-snap, the second-click zoom and the annulus button all work
   without knowing which search produced them. */
void IFSGui::uv_show_results(const std::vector<fd_landmark>& found,
                             const std::string& lead, const char* what) {
  mand_landmark_list = found;
  mand_landmark_list_from_uv = true;
  mand_landmark_selected = 0;
  if (!mand_landmarks) {                 //no point finding them and not showing them
    mand_landmarks = true;
    W_mand_landmarks_check.checked = true;
    W_mand_landmarks_check.redraw();
  }
  /* take the view to the first root if it is off screen, exactly as typing a parameter
     into the s = box does -- otherwise the answer is invisible */
  cpx sig(found[0].sigma_re, found[0].sigma_im);
  IFS.set_params(sig, sig);
  if (sig.real() < mand_ll.real() || sig.real() > mand_ur.real() ||
      sig.imag() < mand_ll.imag() || sig.imag() > mand_ur.imag()) {
    mand_recenter();
    recompute_point_data();        //mand_recenter does not, change_highlighted_ifs does
  } else {
    change_highlighted_ifs(sig);   //this recomputes the point data itself
    draw_mand();
  }

  std::stringstream T;
  T.precision(15);
  T << lead << " (a=" << found[0].a << " b=" << found[0].b << "): " << found.size()
    << " " << what << (found.size() == 1 ? "" : "s") << " -- ";
  for (int i = 0; i < (int)found.size(); ++i) {
    if (i) T << " ; ";
    T << "[" << i << "] " << found[i].sigma_re << (found[i].sigma_im < 0 ? "" : "+")
      << found[i].sigma_im << "i";
  }
  T << ".  Click one to select it, twice to zoom in; " << found[0].spec
    << " names the first on the command line";
  set_status(T.str());
}

/* The `roots` search: the s at which two FINITE words agree, u(0) = v(0).
 *
 * Unlike a landmark, where the coincidence is asymptotic, here it is exact and finite -- so
 * sigma is a root of a {0,+-1} polynomial and there are generally several, a Galois orbit,
 * of which the useful ones are those with 1/2 < |sigma| < 1 and Im sigma > 0.  Since u and v
 * then denote the same affine map, sigma is a renormalization point with b = 1, so everything
 * downstream (including the annulus run) works at it exactly as at a landmark.
 *
 * The solver is funddom's own coincidence_data, reached through fd_core_from_coin, so the
 * specs reported here (coin:u:v:k) can be handed straight to the command line. */
void IFSGui::S_point_roots(XEvent* e) {
  if (e->type != ButtonPress) return;
  do_uv_roots_search();
}

void IFSGui::do_uv_roots_search() {
  std::string u, v, why;
  if (!parse_uv_finite(W_point_uv_entry.text, u, v, why)) {
    set_status("roots: " + why);
    return;
  }
  fd_core c0;
  int n = fd_core_from_coin(u.c_str(), v.c_str(), 0, &c0);
  if (n == -2) { set_status("roots: those words are too long"); return; }
  if (n < 0) {
    set_status("roots: u = " + u + ", v = " + v + " is not a legal coincidence pair");
    return;
  }
  if (n == 0) {
    set_status("roots: u = " + u + ", v = " + v + " has no root with 1/2 < |s| < 1 and"
               " Im s > 0 -- the coincidence polynomial has roots, but none in the annulus"
               " where the attractor can be connected with interior");
    return;
  }
  std::vector<fd_landmark> found;
  for (int pick = 0; pick < n; ++pick) {
    fd_core c2;
    if (fd_core_from_coin(u.c_str(), v.c_str(), pick, &c2) <= 0) continue;
    fd_landmark L;
    L.sigma_re = c2.sigma_re; L.sigma_im = c2.sigma_im;
    L.a = c2.a; L.b = c2.b; L.deg = c2.a;
    std::snprintf(L.spec, sizeof L.spec, "coin:%s:%s:%d", u.c_str(), v.c_str(), pick);
    found.push_back(L);
  }
  if (found.empty()) { set_status("roots: could not rebuild the roots"); return; }
  uv_show_results(found, "roots of u = " + u + ", v = " + v, "root");
}

void IFSGui::S_point_connected(XEvent* e) {
  if (e->type != ButtonPress) return;
  point_connected_check = !point_connected_check;
  W_point_connected_check.checked = point_connected_check;
  W_point_connected_check.redraw();
  recompute_point_data();
}

void IFSGui::S_point_connected_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_connected_depth+1, 1, MAX_POINT_DEPTH);
  if (d == point_connected_depth) return;
  point_connected_depth = d;
  std::stringstream T; T.str(""); T << point_connected_depth;
  W_point_connected_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_connected_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_connected_depth-1, 1, MAX_POINT_DEPTH);
  if (d == point_connected_depth) return;
  point_connected_depth = d;
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
  int d = clamp_int(point_contains_half_depth+1, 1, MAX_POINT_DEPTH);
  if (d == point_contains_half_depth) return;
  point_contains_half_depth = d;
  std::stringstream T; T.str(""); T << point_contains_half_depth;
  W_point_contains_half_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_contains_half_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_contains_half_depth-1, 1, MAX_POINT_DEPTH);
  if (d == point_contains_half_depth) return;
  point_contains_half_depth = d;
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
  int d = clamp_int(point_uv_words_depth+1, 1, MAX_POINT_DEPTH);
  if (d == point_uv_words_depth) return;
  point_uv_words_depth = d;
  std::stringstream T; T.str(""); T << point_uv_words_depth;
  W_point_uv_words_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_uv_words_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_uv_words_depth-1, 1, MAX_POINT_DEPTH);
  if (d == point_uv_words_depth) return;
  point_uv_words_depth = d;
  std::stringstream T; T.str(""); T << point_uv_words_depth;
  W_point_uv_words_depth_label.update_text(T.str());
  recompute_point_data();
}


void IFSGui::S_point_trap(XEvent* e) {
  if (e->type != ButtonPress) return;
  point_trap_check = !point_trap_check;
  W_point_trap_check.checked = point_trap_check;
  W_point_trap_check.redraw();
  recompute_point_data();
}

void IFSGui::S_point_trap_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_trap_depth+1, 1, MAX_POINT_DEPTH);
  if (d == point_trap_depth) return;
  point_trap_depth = d;
  std::stringstream T; T.str(""); T << point_trap_depth;
  W_point_trap_depth_label.update_text(T.str());
  recompute_point_data();
}

void IFSGui::S_point_trap_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(point_trap_depth-1, 1, MAX_POINT_DEPTH);
  if (d == point_trap_depth) return;
  point_trap_depth = d;
  std::stringstream T; T.str(""); T << point_trap_depth;
  W_point_trap_depth_label.update_text(T.str());
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
      /* >= 0, not > 0.  A pixel is in M exactly when this entry is non-negative -- that is
         what mand_get_color tests, and what the connectedness pass writes (-1 for "not
         connected", otherwise a shade).  Testing > 0 threw away every pixel whose shade is
         0, which is most of M: those are the ones that render black.  The trace then
         started inside what it took to be a hole and returned a path bounding nothing,
         while the movie path just below already used >= 0.  Same grid, two conventions. */
      grid[i][j] = (mand_data_grid[i][mand_num_pixel_groups-j-1].x >= 0);
    }
  }
  Point2d<int> p = mand_cpx_to_pixel_group(IFS.z);
  p.y = mand_num_pixel_groups - p.y -1;
  path = IFSPath();
  IFS.hole_boundary_containing_point_from_grid(path.path, path.closed, grid, p, mand_ll, mand_ur, 0);
  /* An empty trace is the normal outcome when the selected parameter is not inside a hole
     of M at all, and it used to be indistinguishable from success: path.is_valid went true
     regardless and the whole "Path options" panel appeared, its buttons all operating on
     nothing.  Say what happened and leave the panel alone. */
  if (path.path.size() == 0) {
    path = IFSPath();
    set_status(mand_connected
      ? "no hole boundary here: the selected parameter is not inside a component of the"
        " complement of M in this window -- pick one in a white region and try again"
      : "no hole boundary here: the Connected layer is off, so there is no computed grid"
        " of M to trace the boundary of -- switch it on first");
    return;
  }
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
  //create the mandelbrot connectedness grid
  std::vector<std::vector<bool> > mand_connected_grid;
  if (path.movie_with_mandelbrot) {
    mand_connected_grid = std::vector<std::vector<bool> >(mand_data_grid.size(), std::vector<bool>(mand_data_grid.size(), false));
    for (int i=0; i<(int)mand_data_grid.size(); ++i) {
      for (int j=0; j<(int)mand_data_grid.size(); ++j) {
        mand_connected_grid[i][mand_data_grid.size()-j-1] = (mand_data_grid[i][j].x >= 0);
      }
    }
  }
  /* ifs_movie_from_path refuses more than 2000 frames, and said so only on stdout -- so
     from the GUI the button simply did nothing.  Check here, where the status line is. */
  long frames = (long)path.movie_fps*(long)path.movie_length;
  if (frames > 2000) {
    std::stringstream T;
    T << "movie: " << path.movie_fps << " fps x " << path.movie_length << " s = " << frames
      << " frames, and at most 2000 can be written -- lower the movie length to "
      << (2000/(path.movie_fps > 0 ? path.movie_fps : 1)) << " s or less";
    set_status(T.str());
    return;
  }
  { std::stringstream T;
    T << "movie: writing " << frames << " frames to ifs_movie* -- this cannot be stopped"
         " once started";
    set_status(T.str()); }
  bool movie_ok = ifs_movie_from_path(IFS, path.path, path.closed, "ifs_movie",
                            limit_ll, limit_ur, limit_depth,
                            path.movie_with_mandelbrot, &mand_ll, &mand_ur, &mand_connected_grid,
                            W_limit_plot.width, W_limit_plot.height,
                            path.movie_fps, path.movie_length, 1);
  set_status(movie_ok ? "movie: wrote ifs_movie*" : "movie: failed (see the terminal)");
}

void IFSGui::S_mand_path_movie_decrease_length(XEvent* e) {
  if (e->type != ButtonPress) return;
  //the old guard was >= 1, which let the length land on 0; movie.cc then
  //divides by it and its frame loop never terminates
  int n = clamp_int(path.movie_length-1, 1, MAX_MOVIE_LENGTH);
  if (n == path.movie_length) return;
  path.movie_length = n;
  std::stringstream T; T.str(""); T << path.movie_length;
  W_mand_path_movie_length_label.update_text(T.str());
}
void IFSGui::S_mand_path_movie_increase_length(XEvent* e) {
  if (e->type != ButtonPress) return;
  int n = clamp_int(path.movie_length+1, 1, MAX_MOVIE_LENGTH);
  if (n == path.movie_length) return;
  path.movie_length = n;
  std::stringstream T; T.str(""); T << path.movie_length;
  W_mand_path_movie_length_label.update_text(T.str());
}

void IFSGui::S_mand_path_movie_with_mandelbrot(XEvent* e) {
  if (e->type != ButtonPress) return;
  path.movie_with_mandelbrot = !path.movie_with_mandelbrot;
  W_mand_path_movie_with_mandelbrot.checked = path.movie_with_mandelbrot;
  W_mand_path_movie_with_mandelbrot.redraw();
}

void IFSGui::S_mand_path_find_uv_words(XEvent* e) {
  if (e->type != ButtonPress || !path.is_valid) return;
  ifs temp_IFS;
  temp_IFS.find_closest_uv_words_along_path(path.path, path.closed, point_uv_words_depth); 
}


void IFSGui::S_mand_path_find_half_words(XEvent* e) {
  if (e->type != ButtonPress || !path.is_valid) return;
  ifs temp_IFS;
  path.half_words = temp_IFS.get_certified_half_balls_along_path(path.path, limit_depth, 1);
  path.has_half_words = (path.half_words.size() > 0);
  if (path.has_half_words) {
    path.half_start = 0;
    std::stringstream T; T.str(""); T << path.half_start;
    W_mand_path_half_start_label.update_text(T.str());
    path.half_end = path.half_words.size()-1;
    T.str(""); T << path.half_end;
    W_mand_path_half_end_label.update_text(T.str());
    draw_mand();
  }
}


void IFSGui::S_mand_path_half_increase_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(path.half_depth+1, 1, MAX_WORD_DEPTH);
  if (d == path.half_depth) return;
  path.half_depth = d;
  std::stringstream T; T.str(""); T << path.half_depth;
  W_mand_path_half_depth_label.update_text(T.str());
  draw_mand();
}

void IFSGui::S_mand_path_half_decrease_depth(XEvent* e) {
  if (e->type != ButtonPress) return;
  int d = clamp_int(path.half_depth-1, 1, MAX_WORD_DEPTH);
  if (d == path.half_depth) return;
  path.half_depth = d;
  std::stringstream T; T.str(""); T << path.half_depth;
  W_mand_path_half_depth_label.update_text(T.str());
  draw_mand();
}

//the four start/stop steppers below only wrapped when has_half_words was set,
//so before "Find half words along path" has run they walked away from 0 with
//nothing to bound them -- and each step paid for a full draw_mand
void IFSGui::S_mand_path_half_increase_start(XEvent* e) {
  if (e->type != ButtonPress || !path.has_half_words) return;
  ++path.half_start;
  if (path.half_start >= (int)path.half_words.size()) path.half_start = 0;
  std::stringstream T; T.str(""); T << path.half_start;
  W_mand_path_half_start_label.update_text(T.str());
  draw_mand();
}

void IFSGui::S_mand_path_half_decrease_start(XEvent* e) {
  if (e->type != ButtonPress || !path.has_half_words) return;
  --path.half_start;
  if (path.half_start < 0) path.half_start = path.half_words.size() - 1;
  std::stringstream T; T.str(""); T << path.half_start;
  W_mand_path_half_start_label.update_text(T.str());
  draw_mand();
}

void IFSGui::S_mand_path_half_increase_end(XEvent* e) {
  if (e->type != ButtonPress || !path.has_half_words) return;
  ++path.half_end;
  if (path.half_end >= (int)path.half_words.size()) path.half_end = 0;
  std::stringstream T; T.str(""); T << path.half_end;
  W_mand_path_half_end_label.update_text(T.str());
  draw_mand();
}

void IFSGui::S_mand_path_half_decrease_end(XEvent* e) {
  if (e->type != ButtonPress || !path.has_half_words) return;
  --path.half_end;
  if (path.half_end < 0) path.half_end = path.half_words.size()-1;
  std::stringstream T; T.str(""); T << path.half_end;
  W_mand_path_half_end_label.update_text(T.str());
  draw_mand();
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
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_with_mandelbrot);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_uv_words_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_half_words_button);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_half_depth_title);
  pack_widget_upper_right(&W_mand_path_half_depth_title, &W_mand_path_half_depth_leftarrow);
  pack_widget_upper_right(&W_mand_path_half_depth_leftarrow, &W_mand_path_half_depth_label);
  pack_widget_upper_right(&W_mand_path_half_depth_label, &W_mand_path_half_depth_rightarrow);
  pack_widget_upper_right(&W_mand_plot, &W_mand_path_half_start_stop_title);
  pack_widget_upper_right(&W_mand_path_half_start_stop_title, &W_mand_path_half_start_leftarrow);
  pack_widget_upper_right(&W_mand_path_half_start_leftarrow, &W_mand_path_half_start_label);
  pack_widget_upper_right(&W_mand_path_half_start_label, &W_mand_path_half_start_rightarrow);
  pack_widget_upper_right(&W_mand_path_half_start_rightarrow, &W_mand_path_half_end_leftarrow);
  pack_widget_upper_right(&W_mand_path_half_end_leftarrow, &W_mand_path_half_end_label);
  pack_widget_upper_right(&W_mand_path_half_end_label, &W_mand_path_half_end_rightarrow);
  
  W_mand_path_tasks_title.initial_draw();
  W_mand_path_delete_button.initial_draw();
  W_mand_path_find_traps_button.initial_draw();
  W_mand_path_create_movie_button.initial_draw();
  W_mand_path_movie_length_title.initial_draw();
  W_mand_path_movie_decrease_length.initial_draw();
  std::stringstream T; T.str(""); T << path.movie_length;
  W_mand_path_movie_length_label.update_text(T.str());
  W_mand_path_movie_with_mandelbrot.checked = path.movie_with_mandelbrot;
  W_mand_path_movie_with_mandelbrot.initial_draw();
  W_mand_path_movie_with_mandelbrot.redraw();
  W_mand_path_movie_increase_length.initial_draw();
  W_mand_path_find_uv_words_button.initial_draw();
  W_mand_path_find_half_words_button.initial_draw();
  W_mand_path_half_depth_title.initial_draw();
  W_mand_path_half_depth_leftarrow.initial_draw();
  T.str(""); T << path.half_depth;
  W_mand_path_half_depth_label.update_text(T.str());
  W_mand_path_half_depth_rightarrow.initial_draw();
  W_mand_path_half_start_stop_title.initial_draw();
  W_mand_path_half_start_leftarrow.initial_draw();
  T.str(""); T << path.half_start;
  W_mand_path_half_start_label.update_text(T.str());
  W_mand_path_half_start_rightarrow.initial_draw();
  W_mand_path_half_end_leftarrow.initial_draw();
  T.str(""); T << path.half_end;
  W_mand_path_half_end_label.update_text(T.str());
  W_mand_path_half_end_rightarrow.initial_draw();
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
    detach_widget(&W_mand_path_movie_with_mandelbrot);
    detach_widget(&W_mand_path_find_uv_words_button);
    detach_widget(&W_mand_path_find_half_words_button);
    detach_widget(&W_mand_path_half_depth_title);
    detach_widget(&W_mand_path_half_depth_leftarrow);
    detach_widget(&W_mand_path_half_depth_label);
    detach_widget(&W_mand_path_half_depth_rightarrow);
    detach_widget(&W_mand_path_half_start_stop_title);
    detach_widget(&W_mand_path_half_start_leftarrow);
    detach_widget(&W_mand_path_half_start_label);
    detach_widget(&W_mand_path_half_start_rightarrow);
    detach_widget(&W_mand_path_half_end_leftarrow);
    detach_widget(&W_mand_path_half_end_label);
    detach_widget(&W_mand_path_half_end_rightarrow);
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

void IFSGui::draw_nifs_limit() {

  //create the nifs object
  //nIFS nifs(4, IFS.z);
  nIFS nifs(3,IFS.z);
  //nifs.centers[3] = -2; nifs.centers[4] = 2;
  nifs.centers[0] = -1; nifs.centers[1] = 0; nifs.centers[2] = 1;
  double min_r = nifs.minimal_initial_radius();
  
  //std::cout << "Got min initial radius of " << min_r << "\n";
  
  if (min_r < 0) return;

  Widget& LW = W_limit_plot;
  
  //clear the limit widget
  XSetForeground(display, LW.gc, WhitePixel(display, screen));
  XFillRectangle(display, LW.p, LW.gc, 0, 0, LW.width, LW.height);
  XSetForeground(display, LW.gc, BlackPixel(display, screen));
  XDrawRectangle(display, LW.p, LW.gc, 0, 0, LW.width-1, LW.height-1);
  XSetFillStyle(display, LW.gc, FillSolid);
  
  std::vector<int> colors(5);
  colors[0] = get_rgb_color(1,0,0);
  colors[1] = get_rgb_color(0,1,0);
  colors[2] = get_rgb_color(0,0,1);
  colors[3] = get_rgb_color(1,0,1);
  colors[4] = get_rgb_color(0,1,1);
  
  nBall initial_ball(0,min_r,1);

  /* The point the picture exists to locate -- and it is 1/s, not 2/s.
     Lambda is connected iff f(Lambda) meets g(Lambda), i.e. iff s a - 1 = s b + 1 for some
     a, b in Lambda, i.e. iff 2/s lies in Lambda - Lambda.  But the set drawn here is the
     attractor E of z -> s z + c with c in {-1, 0, +1}, whereas the differences of the two
     translations are {-2, 0, +2}: so Lambda - Lambda = 2E, and the test on the plotted set
     is 1/s in E.  (Checked numerically as well as algebraically: on |s| = 1/sqrt2, where s
     is certainly in M, 1/s lands in E at 45, 60 and 84 degrees and 2/s lands outside.) */
  cpx crit_point = 1.0/IFS.z;

  std::vector< nBall_stuff > stack(0);
  stack.push_back( nBall_stuff(false, -1, 0, initial_ball) );
  while (stack.size() > 0) {
    nBall_stuff bs = stack.back();
    stack.pop_back();
    //if the ball is disjoint from the window, we might as well get rid of it
    if (!bs.contained && bs.ball.is_disjoint(limit_ll, limit_ur)) continue;
    if ( bs.depth >= limit_depth || bs.ball.radius < limit_pixel_width/2.0 ) {
      Point2d<int> p = limit_cpx_to_pixel(bs.ball.center);
      double r = bs.ball.radius / limit_pixel_width;
      if (r <= 1.0) r = 1.0;
      if (limit_colors) {
        XSetForeground(display, LW.gc, colors[bs.last_gen]);
      }
      if (limit_chunky) {
        XFillArc(display, LW.p, LW.gc, p.x-r, p.y-r, int(2.0*r), int(2.0*r), 23040, 23040);
      } else {
        XDrawPoint(display, LW.p, LW.gc, p.x, p.y);
      }
      continue;
    }
    //if the ball isn't disjoint from the window, maybe it is contained in it?
    for (int i=0; i<3; ++i) {
      nBall new_ball = nifs.act_on_right(i, bs.ball);
      if (bs.contained) {
        stack.push_back(nBall_stuff(true, (bs.last_gen==-1 ? i : bs.last_gen), bs.depth+1, new_ball));
      } else {
        if (bs.ball.is_contained(limit_ll, limit_ur)) {
          stack.push_back(nBall_stuff(true, (bs.last_gen==-1 ? i : bs.last_gen), bs.depth+1, new_ball));
        } else {
          stack.push_back(nBall_stuff(false, (bs.last_gen==-1 ? i: bs.last_gen), bs.depth+1, new_ball));
        }
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

  /* the cross at 1/s, drawn last so it is never buried under the set */
  if (limit_ll.real() < crit_point.real() && crit_point.real() < limit_ur.real() &&
      limit_ll.imag() < crit_point.imag() && crit_point.imag() < limit_ur.imag()) {
    Point2d<int> p = limit_cpx_to_pixel(crit_point);
    XSetForeground(display, LW.gc, BlackPixel(display, screen));
    XDrawLine(display, LW.p, LW.gc, p.x-6, p.y, p.x+6, p.y);
    XDrawLine(display, LW.p, LW.gc, p.x, p.y-6, p.x, p.y+6);
  }

  LW.redraw();
}


void IFSGui::draw_gifs_limit() {

  double PI = 3.14159265358979323846;

  //create the gifs object
  std::vector<cpx> c(4);
  std::vector<cpx> f(4);
  c[0] = cpx(0,1.0); 
  c[1] = cpx(-sin(PI/3.0), -cos(PI/3.0)); 
  c[2] = cpx(sin(PI/3.0), -cos(PI/3.0));
  c[3] = cpx(0,0);
  f[0] = IFS.z.real(); f[1] = IFS.z.real(); f[2] = IFS.z.real();
  f[3] = IFS.z.imag();
  
  gIFS gifs(f,c);
  
  gBall initial_ball;
  
  if (!gifs.minimal_ball(initial_ball)) return;

  Widget& LW = W_limit_plot;
  
  //clear the limit widget
  XSetForeground(display, LW.gc, WhitePixel(display, screen));
  XFillRectangle(display, LW.p, LW.gc, 0, 0, LW.width, LW.height);
  XSetForeground(display, LW.gc, BlackPixel(display, screen));
  XDrawRectangle(display, LW.p, LW.gc, 0, 0, LW.width-1, LW.height-1);
  XSetFillStyle(display, LW.gc, FillSolid);
  
  std::vector<int> colors(5);
  colors[0] = get_rgb_color(1,0,0);
  colors[1] = get_rgb_color(0,1,0);
  colors[2] = get_rgb_color(0,0,1);
  colors[3] = get_rgb_color(1,0,1);
  
  std::vector< gBall_stuff > stack(0);
  stack.push_back( gBall_stuff(false, -1, 0, initial_ball) );
  while (stack.size() > 0) {
    gBall_stuff bs = stack.back();
    stack.pop_back();
    //if the ball is disjoint from the window, we might as well get rid of it
    if (!bs.contained && bs.ball.is_disjoint(limit_ll, limit_ur)) continue;
    if ( bs.depth >= limit_depth || bs.ball.radius < limit_pixel_width/2.0 ) {
      Point2d<int> p = limit_cpx_to_pixel(bs.ball.center);
      double r = bs.ball.radius / limit_pixel_width;
      if (r <= 1.0) r = 1.0;
      if (limit_colors) {
        XSetForeground(display, LW.gc, colors[bs.last_gen]);
      }
      if (limit_chunky) {
        XFillArc(display, LW.p, LW.gc, p.x-r, p.y-r, int(2.0*r), int(2.0*r), 23040, 23040);
      } else {
        XDrawPoint(display, LW.p, LW.gc, p.x, p.y);
      }
      continue;
    }
    //if the ball isn't disjoint from the window, maybe it is contained in it?
    for (int i=0; i<gifs.centers.size(); ++i) {
      gBall new_ball = gifs.act_on_right(i, bs.ball);
      if (bs.contained) {
        stack.push_back(gBall_stuff(true, (bs.last_gen==-1 ? i : bs.last_gen), bs.depth+1, new_ball));
      } else {
        if (bs.ball.is_contained(limit_ll, limit_ur)) {
          stack.push_back(gBall_stuff(true, (bs.last_gen==-1 ? i : bs.last_gen), bs.depth+1, new_ball));
        } else {
          stack.push_back(gBall_stuff(false, (bs.last_gen==-1 ? i: bs.last_gen), bs.depth+1, new_ball));
        }
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
  
  LW.redraw();
}






void IFSGui::draw_2d_limit() {
  
  Widget& LW = W_limit_plot;
  //clear the limit widget
  XSetForeground(display, LW.gc, WhitePixel(display, screen));
  XFillRectangle(display, LW.p, LW.gc, 0, 0, LW.width, LW.height);
  XSetForeground(display, LW.gc, BlackPixel(display, screen));
  XDrawRectangle(display, LW.p, LW.gc, 0, 0, LW.width-1, LW.height-1);
  XSetFillStyle(display, LW.gc, FillSolid);
  
  //set up the IFS and get the initial rectangle 
  double rz = std::real(IFS.z);
  if (rz < 0.6 || rz > 0.68) {
    /* Say so.  This used to return silently, leaving a blank white pane and no clue that
       the layer had simply declined to draw -- which reads as a broken feature. */
    std::stringstream T;
    T.precision(6);
    T << "2d IFS: only defined for 0.6 <= Re(s) <= 0.68 (the box sizes are tabulated over"
         " that range); Re(s) is " << rz << " here";
    set_status(T.str());
    LW.redraw();
    return;
  }

  ifs2d IFS2d;
  IFS2d.gens.resize(2);
  IFS2d.gens[0] = AffineMap( rz, 0, 1, rz, 0, 0 );
  IFS2d.gens[1] = AffineMap( rz, 0, 1, rz, 1-rz, -1 );
  double box_sizes[17] = {1.98593, 2.08085, 2.1883, 2.31103, 2.45263, 2.61793, 2.81353,
                          3.04875, 3.33714, 3.69922, 4.16762, 4.79755, 5.69039, 7.05483,
                          9.46713, 14.6384, 33.0532};
  /* Index directly and clamp.  Accumulating 0.005 at a time and counting the steps lands
     on 16 at rz = 0.68 only because 0.6 + 16*0.005 happens to round up to >= 0.68; one ulp
     the other way reads box_sizes[17], one past the end. */
  int ub_ind = (int)std::ceil((rz - 0.6)/0.005);
  ub_ind = clamp_int(ub_ind, 0, 16);
  double bs = box_sizes[ub_ind];
  
  std::vector<Point2d<double> > initial_box(4);
  initial_box[0] = Point2d<double>(0.5-bs, 0-bs);
  initial_box[1] = Point2d<double>(0.5+bs, 0-bs);
  initial_box[2] = Point2d<double>(0.5+bs, 0+bs);
  initial_box[3] = Point2d<double>(0.5-bs, 0+bs);
    
  //figure out where the box goes, as averages of its own boundary points,
  //under all the words of length 6
  std::vector<std::vector<Point4d<double> > > targets_under_words(64);
  //std::cout << "Original box: " << initial_box[0] << " " << initial_box[1] << " " << initial_box[2] << " " << initial_box[3] << "\n";
  for (int i=0; i<64; ++i) {
    AffineMap A = IFS2d.semigroup_element(i, 6);
    targets_under_words[i].resize(4);
    //std::cout << i << ": box points: ";
    for (int j=0; j<4; ++j) {
      Point2d<double> point_target = A(initial_box[j]);
      targets_under_words[i][j] = point_as_weighted_average_in_box(point_target, Point2d<double>(0.5,0), bs);
      //std::cout << point_target << " weights: " << targets_under_words[i][j] << "; ";
    }
    //std::cout << "\n";
  }
  Point4d<double> zero_weights = point_as_weighted_average_in_box(Point2d<double>(0,0), Point2d<double>(0.5,0), bs);
  
  //set the colors
  int bcolor = get_rgb_color(0,0,1);
  int rcolor = get_rgb_color(1,0,0);
  
  std::vector< Box_Stuff > stack(1);
  stack[0] = Box_Stuff(false, -1, 0, 0, initial_box);
  
  while ((int)stack.size() > 0) {
    Box_Stuff bs = stack.back();
    stack.pop_back();
    
    //if the box is disjoint from the window, get rid of it
    if (!bs.contained && bs.is_disjoint(cpx(limit_ll.real(), 5*limit_ll.imag()), 
                                        cpx(limit_ur.real(), 5*limit_ur.imag())) ) continue;
    
    //if we are at our target depth, draw it
    if (bs.depth >= limit_depth) {
      //XPoint v[4];
      //for (int i=0; i<(int)4; ++i) {
      //  Point2d<int> pv = limit_cpx_to_pixel(cpx(bs.box[i].x, bs.box[i].y/5.0));
      //  v[i].x = pv.x; v[i].y = pv.y;
      //}
      //XSetForeground(display, LW.gc, (bs.last_gen == 1 ? rcolor : bcolor));
      //XFillPolygon(display, LW.p, LW.gc, v, 4, Convex, CoordModeOrigin);
      XSetForeground(display, LW.gc, (bs.last_gen == 1 ? rcolor : bcolor));
      Point2d<double> zero_image = zero_weights.x*bs.box[0] + zero_weights.y*bs.box[1] + 
                                   zero_weights.z*bs.box[2] + zero_weights.w*bs.box[3];
      Point2d<int> pv = limit_cpx_to_pixel(cpx(zero_image.x, zero_image.y/5.0));
      XDrawPoint(display, LW.p, LW.gc, pv.x, pv.y);
      //std::cout << "Drew: " << std::bitset<64>(bs.word).to_string() << "\n";
      continue;
    }
    
    //if the box isn't disjoint, we must subdivide
    //we check if it is contained in order to speed this up
    for (int i=0; i<64; ++i) {
      std::vector<Point2d<double> > new_box(4);
      for (int j=0; j<4; ++j) {
        Point4d<double> target_weights = targets_under_words[i][j];
        new_box[j] = target_weights.x*bs.box[0] + target_weights.y*bs.box[1] + 
                     target_weights.z*bs.box[2] + target_weights.w*bs.box[3];
      }
      if (bs.contained) {
        stack.push_back( Box_Stuff(true, 
                                   (bs.last_gen==-1 ? (i>>5)&1 : bs.last_gen), 
                                   (bs.last_gen==-1 ? i : (bs.word<<6)|i),
                                   bs.depth+1, 
                                   new_box) );
      } else {
        stack.push_back( Box_Stuff( bs.is_contained(cpx(limit_ll.real(), 5*limit_ll.imag()), 
                                                    cpx(limit_ur.real(), 5*limit_ur.imag())), 
                                    (bs.last_gen==-1 ? (i>>5)&1 : bs.last_gen), 
                                    (bs.last_gen==-1 ? i : (bs.word<<6)|i),
                                    bs.depth+1, 
                                    new_box) );
      }
    }
  }
  
  
  
  LW.redraw();
  
}




void IFSGui::draw_limit() {
  double min_r;
  if (!IFS.minimal_enclosing_radius(min_r)) {
    return;
  }
  
  if (limit_nifs) {
    draw_nifs_limit();
    return;
  } else if (limit_2d) {
    draw_2d_limit();
    return;
  } else if (limit_gifs) {
    draw_gifs_limit();
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
  
  /* BATCH THE LEAVES.  One leaf ball used to be one X request -- and there are 2^n of
     them, so at depth 22 that is four million requests for a picture of at most 319,225
     pixels.  XDrawPoints and XFillArcs take arrays, so the leaves are collected per
     colour and flushed in blocks; the picture is identical and the request count drops
     by three orders of magnitude.  Two buckets because the colour can only be one of
     two (blue = last generator f, yellow = g), and with `Colors` off only one is used. */
  const int LIMIT_BATCH = 4096;
  std::vector<XPoint> pts[2];
  std::vector<XArc>   arcs[2];
  pts[0].reserve(LIMIT_BATCH); pts[1].reserve(LIMIT_BATCH);
  arcs[0].reserve(LIMIT_BATCH); arcs[1].reserve(LIMIT_BATCH);
  /* With `Colors` off the old code set no foreground at all, so the leaves inherited the
     BlackPixel left by the border above.  Batching has to set a foreground explicitly, so
     that has to be black here or the picture would silently turn blue. */
  int leaf_colors[2];
  if (limit_colors) { leaf_colors[0] = blue_color; leaf_colors[1] = yellow_color; }
  else { leaf_colors[0] = leaf_colors[1] = BlackPixel(display, screen); }

  Ball initial_ball(0.5,(IFS.z-1.0)/2.0,(1.0-IFS.w)/2.0,min_r);
  std::vector<std::pair<bool,Ball> > stack(0);
  stack.push_back(std::make_pair(false, initial_ball));
  /* Escape has to work here too.  The word length goes up to MAX_WORD_DEPTH and the cost
     doubles with every level, so this loop can run for many minutes, and until now the
     only way out was to kill the process.  Polling on a counter rather than per ball
     because XCheckMaskEvent is a good deal dearer than one iteration of this loop. */
  bool limit_aborted = false;
  long balls_done = 0;
  while (stack.size() > 0) {
    if ((++balls_done & 0xFFFF) == 0) {
      XEvent ev;
      while (XCheckMaskEvent(display, KeyPressMask, &ev))
        if (XLookupKeysym(&ev.xkey, 0) == XK_Escape) limit_aborted = true;
      if (limit_aborted) break;
    }
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
      int k = (limit_colors && b.second.last_gen_index() != 0) ? 1 : 0;
      if (limit_chunky) {
        /* (int)(p.x - r), NOT p.x - (int)r.  The XFillArc call this replaces passed the
           double expression p.x - r straight into an int parameter, so it truncated the
           whole thing; rounding r first moves edge pixels by one and the picture is
           visibly not the same. */
        XArc a;
        a.x = (short)(int)(p.x - r); a.y = (short)(int)(p.y - r);
        a.width = (unsigned short)(int)(2.0*r);
        a.height = (unsigned short)(int)(2.0*r);
        a.angle1 = 23040; a.angle2 = 23040;
        arcs[k].push_back(a);
        if ((int)arcs[k].size() >= LIMIT_BATCH) {
          XSetForeground(display, LW.gc, leaf_colors[k]);
          XFillArcs(display, LW.p, LW.gc, &arcs[k][0], (int)arcs[k].size());
          arcs[k].clear();
        }
      } else {
        XPoint q; q.x = (short)p.x; q.y = (short)p.y;
        pts[k].push_back(q);
        if ((int)pts[k].size() >= LIMIT_BATCH) {
          XSetForeground(display, LW.gc, leaf_colors[k]);
          XDrawPoints(display, LW.p, LW.gc, &pts[k][0], (int)pts[k].size(), CoordModeOrigin);
          pts[k].clear();
        }
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
  //whatever is left in the buckets
  for (int k=0; k<2; ++k) {
    if (!arcs[k].empty()) {
      XSetForeground(display, LW.gc, leaf_colors[k]);
      XFillArcs(display, LW.p, LW.gc, &arcs[k][0], (int)arcs[k].size());
    }
    if (!pts[k].empty()) {
      XSetForeground(display, LW.gc, leaf_colors[k]);
      XDrawPoints(display, LW.p, LW.gc, &pts[k][0], (int)pts[k].size(), CoordModeOrigin);
    }
  }
  if (limit_aborted)
    set_status("limit set: stopped -- the picture is incomplete; lower the depth");

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
    XFontStruct* font = gui_font();
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
      cpx& c = uv_graph_balls[i].center;
      if (c.real() < limit_ll.real() || 
          c.real() > limit_ur.real() || 
          c.imag() < limit_ll.imag() || 
          c.imag() > limit_ur.imag()) {
        continue;
      }
      Point2d<int> p = limit_cpx_to_pixel(c);
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
    XSetLineAttributes(display, LW.gc, 1, LineSolid, 1, 1);
    for (int i=0; i<(int)uv_graph_edges.size(); ++i){
      Point3d<int>& e = uv_graph_edges[i];
      cpx& cc1 = uv_graph_balls[e.x].center;
      cpx& cc2 = uv_graph_balls[e.y].center;
      if ((cc1.real() < limit_ll.real() || 
           cc1.real() > limit_ur.real() || 
           cc1.imag() < limit_ll.imag() || 
           cc1.imag() > limit_ur.imag()) &&
          (cc2.real() < limit_ll.real() || 
           cc2.real() > limit_ur.real() || 
           cc2.imag() < limit_ll.imag() || 
           cc2.imag() > limit_ur.imag())) {
        continue;
      }
      Point2d<int> c1 = limit_cpx_to_pixel(cc1);
      Point2d<int> c2 = limit_cpx_to_pixel(cc2);
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
  
  if (limit_trap) draw_limit_trap();

  LW.redraw();
}


/* Draw the trap certifying the current parameter, on top of the limit set.
   What makes the picture a proof rather than a decoration is the INTERLEAVING: it is
   not enough that the two balls u(B) and v(B) overlap, the two Cantor clusters inside
   them have to link, so that no perturbation can pull them apart.  So we draw both the
   enclosing balls (outlined, so the limit set stays visible through them) and enough of
   each cluster's descendants to see them alternate.

   The balls are rebuilt from the words rather than carried out of the search, because
   check_TLB reports Bitwords.  act_on_right shifts the word left and ORs the new letter
   into bit 0, and Bitword::reverse_get(n) returns the n-th letter reading outermost
   first, so replaying reverse_get(0..len-1) through act_on_right reproduces exactly the
   ball the search found.  Getting that order backwards would draw a different, wrong
   pair of balls and look entirely plausible. */
void IFSGui::draw_limit_trap() {
  Widget& LW = W_limit_plot;
  std::vector<std::pair<Bitword,Bitword> > tw;
  std::string why;
  limit_trap_located = false;
  if (!find_trap_words(limit_trap_depth, tw, why)) {
    set_status("trap: " + why);
    return;
  }
  double min_r;
  if (!IFS.minimal_enclosing_radius(min_r)) return;
  /* the same seed ball as check_TLB uses, or the words would not name these balls */
  Ball base(0.5, (IFS.z-1.0)/2.0, (1.0-IFS.w)/2.0, 1.01*min_r);
  const Bitword& u = tw[0].first;
  const Bitword& v = tw[0].second;
  Ball bu = base, bv = base;
  for (int i=0; i<u.len; ++i) bu = IFS.act_on_right(u.reverse_get(i), bu);
  for (int i=0; i<v.len; ++i) bv = IFS.act_on_right(v.reverse_get(i), bv);

  /* Descendants, to show the clusters rather than just their hulls.  How many levels is
     not a free choice: each level shrinks the balls by |s|, so drawing past the point
     where they fall below a pixel costs 2^k arcs and adds nothing, while stopping too
     early at a deep zoom leaves a blank disk.  So take the depth from the geometry --
     enough levels for the smallest ball to be about a pixel -- and cap it, since this
     overlay is redrawn on every pan. */
  int cluster_depth = 0;
  if (IFS.az > 0.0 && IFS.az < 1.0) {
    double want = std::log(bu.radius / limit_pixel_width) / std::log(1.0/IFS.az);
    if (want > 0.0) cluster_depth = (int)std::floor(want);
  }
  if (cluster_depth > 10) cluster_depth = 10;      /* 1024 arcs a side */
  std::vector<Ball> CU, CV;
  IFS.compute_balls_right(CU, bu, cluster_depth);
  IFS.compute_balls_right(CV, bv, cluster_depth);

  int u_col = get_rgb_color(0.85, 0.0, 0.85);      /* magenta: the u side  */
  int v_col = get_rgb_color(0.0, 0.55, 0.15);      /* green:   the v side  */
  int line_col = BlackPixel(display, screen);

  /* Interleave the two clusters rather than drawing one and then the other: the two
     overlap, X11 has no alpha, and whichever went second would hide the other
     everywhere they meet -- which is exactly the region the picture is about. */
  int nmax = (int)(CU.size() > CV.size() ? CU.size() : CV.size());
  for (int i=0; i<nmax; ++i) {
    for (int side=0; side<2; ++side) {
      const std::vector<Ball>& C = (side==0 ? CU : CV);
      if (i >= (int)C.size()) continue;
      XSetForeground(display, LW.gc, (side==0 ? u_col : v_col));
      Point2d<int> p = limit_cpx_to_pixel(C[i].center);
      double r = C[i].radius / limit_pixel_width;
      if (r < 1.0) r = 1.0;
      XFillArc(display, LW.p, LW.gc, p.x-int(r), p.y-int(r),
               int(2.0*r), int(2.0*r), 0, 23040);
    }
  }

  /* the two enclosing balls, outlined so what is underneath still reads */
  for (int pass=0; pass<2; ++pass) {
    const Ball& B = (pass==0 ? bu : bv);
    XSetForeground(display, LW.gc, (pass==0 ? u_col : v_col));
    Point2d<int> p = limit_cpx_to_pixel(B.center);
    double r = B.radius / limit_pixel_width;
    if (r < 2.0) r = 2.0;
    XDrawArc(display, LW.p, LW.gc, p.x-int(r), p.y-int(r),
             int(2.0*r), int(2.0*r), 0, 23040);
    XDrawArc(display, LW.p, LW.gc, p.x-int(r)-1, p.y-int(r)-1,
             int(2.0*r)+2, int(2.0*r)+2, 0, 23040);
    XSetForeground(display, LW.gc, line_col);
    XDrawString(display, LW.p, LW.gc, p.x+int(r)+3, p.y,
                (pass==0 ? "u" : "v"), 1);
  }

  /* the displacement that the trap-like set has to contain */
  XSetForeground(display, LW.gc, line_col);
  Point2d<int> pu = limit_cpx_to_pixel(bu.center);
  Point2d<int> pv = limit_cpx_to_pixel(bv.center);
  XDrawLine(display, LW.p, LW.gc, pu.x, pu.y, pv.x, pv.y);

  limit_trap_located = true;
  limit_trap_center = (bu.center + bv.center)/2.0;
  limit_trap_radius = (bu.radius > bv.radius ? bu.radius : bv.radius);

  std::stringstream T;
  T << "trap: u=" << u.str() << " v=" << v.str() << "  |u|=" << u.len
    << " (" << tw.size() << " pair" << (tw.size()==1 ? "" : "s") << " found)";
  /* the balls are usually a few pixels across in a default view, so say how to see them */
  double on_screen = 2.0*limit_trap_radius/limit_pixel_width;
  if (on_screen < 40.0) T << " -- press 'to trap' to zoom in";
  set_status(T.str());
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

int IFSGui::mand_get_color(PointNd<4,int>& p) {
  if (mand_trap && p.z > 0) { //use the trap color
    return p.z;
  } else if (mand_contains_half && p.y > 0) {
    return p.y;
  } else if (mand_connected && p.x >= 0) {
    return p.x*0x000001;
  } else {
    return WhitePixel(display, screen);
  }
}

/* Paint one mesh cell into the client-side raster of draw_mand.  A NULL image means
   XCreateImage failed, in which case fall back to drawing the cell straight into the
   pixmap -- slow, but correct, and it keeps the failure from blanking the picture. */
void IFSGui::mand_put_cell(XImage* img, int i, int j, unsigned long col) {
  Widget& MW = W_mand_plot;
  int x0 = i*mand_pixel_group_size, y0 = j*mand_pixel_group_size;
  if (img == NULL) {
    XSetForeground(display, MW.gc, col);
    XFillRectangle(display, MW.p, MW.gc, x0, y0,
                   mand_pixel_group_size, mand_pixel_group_size);
    return;
  }
  for (int dy = 0; dy < mand_pixel_group_size; ++dy) {
    int y = y0 + dy;
    if (y >= MW.height) break;
    for (int dx = 0; dx < mand_pixel_group_size; ++dx) {
      int x = x0 + dx;
      if (x >= MW.width) break;
      XPutPixel(img, x, y, col);
    }
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
  
  

//draw the mandelbrot set
void IFSGui::draw_mand() {
  ifs temp_IFS;
  Widget& MW = W_mand_plot;
  
  //set up the TLB
  std::vector<Ball> TLB;
  double TLB_C,TLB_Z;
  bool found_TLB = false;
  if (mand_trap && !mand_grid_trap_valid) {
    /* NOTE the set_params below is vestigial, and that matters for reading the cache logic:
       TLB_for_region overwrites z with the centre of the region it is given and restores it
       on the way out, so this grid depends only on (mand_ll, mand_ur) and mand_trap_depth,
       NOT on the selected parameter.  That is why mand_grid_trap_valid needs no invalidation
       in change_highlighted_ifs -- a question worth an hour of audit until the assignment
       above is read.  Left in place rather than deleted because it is harmless and removing
       it would silently change which overload of the region search is exercised. */
    temp_IFS.set_params(IFS.z, IFS.z);
    temp_IFS.TLB_for_region(TLB, mand_ll, mand_ur, 16, &TLB_C, &TLB_Z, 0);
    found_TLB = (TLB.size() != 0);
  }
  
  //is any layer actually being recomputed?  if not this is a cheap repaint from
  //the cache and there is no point in advertising progress through it
  bool computing = ((mand_connected && !mand_grid_connected_valid) ||
                    (mand_contains_half && !mand_grid_contains_half_valid) ||
                    (mand_trap && !mand_grid_trap_valid));
  bool aborted = false;
  clock_t tmr_t0 = clock();
  struct timeval tmr_w0; gettimeofday(&tmr_w0, NULL);
  struct timeval tmr_last_blit = tmr_w0;
  int pending_from = 0;          //first column not yet shown on screen

  /* NO PROGRESS REPORT HERE, deliberately.  The picture fills in column by column as it
     is computed, so the user can watch the progress directly -- a bar saying the same
     thing in words is redundant, and it was not free: redrawing the status text and the
     bar cost about eleven X requests and two XFlushes each time.  The status line still
     says once, up front, that Escape stops it, and the annulus figure (which produces no
     visual feedback at all until it finishes) keeps its bar. */
  if (computing)
    set_status("mandelbrot: computing (Escape to stop)");

  /* ONE XImage FOR THE WHOLE PANE, blitted a column at a time.
   *
   * The cell-by-cell drawing this replaces issued XSetForeground + XFillRectangle per
   * cell -- 638,450 X requests at the finest mesh, for a picture of 319,225 pixels.  That
   * made the repaint X-bound rather than compute-bound: measured at mesh 1 with three
   * layers, 1.12 s of CPU inside this function but 6.0 s of wall clock.  Writing the
   * pixels into a client-side image costs no X traffic at all, and one XPutImage per
   * column puts the finished strip into the pixmap.  XPutPixel rather than raw stores
   * because it is the only way to be correct for every visual and byte order. */
  XImage* raster = XCreateImage(display, DefaultVisual(display, screen),
                                DefaultDepth(display, screen), ZPixmap, 0, NULL,
                                MW.width, MW.height, 32, 0);
  if (raster) {
    raster->data = (char*)malloc((size_t)raster->bytes_per_line * MW.height);
    if (!raster->data) { XDestroyImage(raster); raster = NULL; }
  }
  const unsigned long white = WhitePixel(display, screen);

  for (int i=0; i<(int)mand_num_pixel_groups; ++i) {

    if (computing) {
      //Poll for an abort.  This pane is the only thing in the program that can
      //take minutes, and until now there was no way out of it but to kill the
      //process.  Draining KeyPressMask here does swallow other keystrokes typed
      //during the repaint, which is a fair trade for being able to stop it.
      XEvent ev;
      while (XCheckMaskEvent(display, KeyPressMask, &ev)) {
        if (XLookupKeysym(&ev.xkey, 0) == XK_Escape) aborted = true;
      }
      if (aborted) break;
    }

    for (int j=0; j<(int)mand_num_pixel_groups; ++j) {

      //do the necessary computations
      cpx c = mand_pixel_group_to_cpx(Point2d<int>(i,j));
      /* Nothing outside the closed unit disc is meaningful: for |s| >= 1 the maps
         are not contractions, so there is no attractor and every layer below would
         be reporting on an object that does not exist.  Leaving those pixels the
         background color makes the picture stop at a circle instead of fading into
         black, and saves computing them at all. */
      if (std::abs(c) > 1.0) {
        mand_put_cell(raster, i, j, white);
        continue;
      }
      temp_IFS.set_params(c,c);
      
      if (mand_connected && !mand_grid_connected_valid) {
        //temp_IFS.set_params(c*c,c*c);
        if (true) { //!limit_nifs && !limit_gifs) {
          if (!temp_IFS.is_connected(mand_connected_depth, mand_data_grid[i][j].x) ) {
            mand_data_grid[i][j].x = -1;
          }
        } else if (limit_nifs) {
          nIFS nifs(3, c);
          nifs.centers[0] = -1; nifs.centers[1] = 0; nifs.centers[2] = 1;
          if (!nifs.is_connected(mand_connected_depth, mand_data_grid[i][j].x)) {
            mand_data_grid[i][j].x = -1;
          }
        } else if (limit_gifs) {
         //create the gifs object
          double PI = 3.14159265358979323846;
          std::vector<cpx> cs(4);
          std::vector<cpx> fs(4);
          cs[0] = cpx(0,1.0); 
          cs[1] = cpx(-sin(PI/3.0), -cos(PI/3.0)); 
          cs[2] = cpx(sin(PI/3.0), -cos(PI/3.0));
          cs[3] = cpx(0,0);
          fs[0] = c.real(); fs[1] = c.real(); fs[2] = c.real();
          fs[3] = c.imag();
          gIFS gifs(fs,cs);
          if (!gifs.is_connected(mand_connected_depth, mand_data_grid[i][j].x)) {
            mand_data_grid[i][j].x = -1;
          }
        }
      }
      if (mand_contains_half && !mand_grid_contains_half_valid) {
        //temp_IFS.set_params(sqrt(c), sqrt(c));
        if (temp_IFS.contains_half(mand_contains_half_depth, mand_data_grid[i][j].y)) {
          mand_data_grid[i][j].y = get_rgb_color(0.5, double(mand_data_grid[i][j].y)/100, 0.5);
        } else {
          mand_data_grid[i][j].y = -1;
        }
        //temp_IFS.set_params(c,c);
      }
      if (mand_trap && !mand_grid_trap_valid && found_TLB) {
        double trap_radius;
        int multiplier = 100/mand_trap_depth;
        int diff = multiplier*temp_IFS.check_TLB(TLB,NULL,NULL,trap_radius,NULL,mand_trap_depth);
        mand_data_grid[i][j].z = (diff < 0 ? -1 : get_rgb_color(0, double(diff)/100, 1.0));
      }

      //record the cell in the client-side image; no X traffic
      mand_put_cell(raster, i, j, mand_get_color(mand_data_grid[i][j]));
    }

    /* SHOW THE FINISHED COLUMNS -- but not one blit per column.
     *
     * A blit to the on-screen window is expensive out of all proportion to its size: on
     * XQuartz each one makes the server update the screen, and 565 of them cost about
     * 4 seconds of wall clock against 1.1 seconds of actual computation.  (That is what
     * was left after the earlier fix, which had already cut one blit per CELL down to one
     * per column.)  So accumulate finished columns and blit the whole strip at once, at
     * most every BLIT_SEC, plus whatever remains at the end.  Twenty updates a second is
     * already smoother than the eye needs, and the picture still fills in left to right. */
    const double BLIT_SEC = 0.05;
    struct timeval now; gettimeofday(&now, NULL);
    double since = (now.tv_sec - tmr_last_blit.tv_sec)
                 + 1e-6*(now.tv_usec - tmr_last_blit.tv_usec);
    bool last_col = (i + 1 == (int)mand_num_pixel_groups);
    if (since >= BLIT_SEC || last_col) {
      int gx = pending_from*mand_pixel_group_size;
      int gw = (i - pending_from + 1)*mand_pixel_group_size;
      if (gx + gw > MW.width) gw = MW.width - gx;
      if (gw > 0) {
        if (raster)
          XPutImage(display, MW.p, MW.gc, raster, gx, 0, gx, 0, gw, MW.height);
        XCopyArea(display, MW.p, main_window, MW.gc, gx, 0, gw, MW.height,
                  MW.ul.x + gx, MW.ul.y);
      }
      pending_from = i + 1;
      tmr_last_blit = now;
    }
  }
  /* an abort leaves the last few columns unblitted; the pixmap is what everything else
     reads, so it must still get them even though the picture is being discarded */
  if (raster && pending_from < (int)mand_num_pixel_groups) {
    int gx = pending_from*mand_pixel_group_size;
    int gw = MW.width - gx;
    if (gw > 0) XPutImage(display, MW.p, MW.gc, raster, gx, 0, gx, 0, gw, MW.height);
  }

  if (raster) XDestroyImage(raster);   //frees raster->data too

  if (getenv("SCHOTTKY_TIMING")) {
    struct timeval tv; gettimeofday(&tv, NULL);
    fprintf(stderr, "[TIMING] draw_mand cpu %.3f s  wall %.3f s, %ld columns\n",
            double(clock()-tmr_t0)/CLOCKS_PER_SEC,
            (tv.tv_sec - tmr_w0.tv_sec) + 1e-6*(tv.tv_usec - tmr_w0.tv_usec),
            (long)mand_num_pixel_groups);
  }

  //an aborted pass left most of the grid at its -1 fill, so the layers must stay
  //invalid or the half-finished picture becomes the cached answer forever
  if (!aborted) {
    if (mand_connected && !mand_grid_connected_valid) mand_grid_connected_valid = true;
    if (mand_contains_half && !mand_grid_contains_half_valid) mand_grid_contains_half_valid = true;
    if (mand_trap && !mand_grid_trap_valid) mand_grid_trap_valid = true;
  }

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
    //if the path has half balls, draw them
    //if (path.has_half_words) {
    //  IFS.certify_set_B_path(path.path, limit_depth, 1);
    //}
    if (path.has_half_words) {
      std::vector<std::vector<Ball> > subdivided_balls(0);
      int i=path.half_start;
      int num_done = 0;
      int total_num = path.half_words.size();
      int num_to_do = (path.half_end >= path.half_start ? 
                            path.half_end - path.half_start+1 : 
                            (path.half_end+1) + (path.half_words.size()-path.half_start));
      std::cout << "Num to do: " << num_to_do << "\n";
      do {
        //find the balls
        subdivided_balls.push_back( IFS.subdivide_half_prefix(path.half_words[i], 
                                                              path.path[0],
                                                              path.half_depth,
                                                              mand_ll, mand_ur) );
        std::cout.flush();
        //draw the balls
        int sbs = subdivided_balls.back().size();
        std::cout << "Got the " << sbs << " balls -- about to draw\n";
        int col = get_rgb_color(0, (double)i/(double)total_num,
                                   (double)(total_num-i)/(double)total_num);
        std::cout << "Got color " << col << "\n";
        XSetForeground(display, MW.gc, col);
        for (int j=0; j<(int)subdivided_balls.back().size(); ++j) {
          if (j==0) {
            std::cout << "Drawing ball at " << subdivided_balls.back()[j].center << " of radius " << subdivided_balls.back()[j].radius << 
                                " with ll " << mand_ll << " ur " << mand_ur <<  " pixel width: " << mand_pixel_width << " so pixel radius = " << int(subdivided_balls.back()[j].radius / mand_pixel_width) << "\n";
            std::cout.flush();
          }
          Point2d<int> p = mand_cpx_to_pixel(subdivided_balls.back()[j].center);
          int r = int(subdivided_balls.back()[j].radius / mand_pixel_width);
          if (r < 2) r = 2;
          XFillArc(display, MW.p, MW.gc, p.x-r, p.y-r, 2*r, 2*r, 23040, 23040);
        }  
        MW.redraw();                               
        
        ++i;
        ++num_done;
        if (i == (int)path.half_words.size()) i = 0;
      } while (num_done < num_to_do);
    }


  }

  mand_update_legend();
  mand_draw_scale();
  mand_draw_guide_circles();
  mand_draw_landmarks();
  MW.redraw();

  /* SCHOTTKY_DUMP_STATE: one line describing the marked-point layers, their caches and the
     keys those caches were built against.  This exists so a test can check COHERENCE exactly
     -- "is this layer on?" against "was it on when the cache was built?" -- instead of trying
     to infer it from pixels.  The bug it was written for had a layer on, an empty cache, and a
     key that matched in every field it recorded. */
  if (getenv("SCHOTTKY_DUMP_STATE")) {
    int inwin_lm = 0, inwin_rt = 0;
    for (int i = 0; i < (int)mand_landmark_list.size(); ++i) {
      cpx c(mand_landmark_list[i].sigma_re, mand_landmark_list[i].sigma_im);
      if (c.real() < mand_ll.real() || c.real() > mand_ur.real() ||
          c.imag() < mand_ll.imag() || c.imag() > mand_ur.imag()) continue;
      if (mand_landmark_list[i].spec[0] == 'c') ++inwin_rt; else ++inwin_lm;
    }
    fprintf(stderr,
      "[STATE] lm_on=%d lm_key_on=%d lm_N=%d lm_key_N=%d lm_cache=%d lm_inwin=%d lm_drawn=%d"
      " rt_on=%d rt_key_on=%d rt_deg=%d rt_key_deg=%d rt_cache=%d rt_inwin=%d rt_drawn=%d"
      " tgt=%d lm_key_tgt=%d rt_key_tgt=%d list=%d sel=%d uv=%d"
      " win=%.17g,%.17g,%.17g,%.17g\n",
      mand_landmarks?1:0, mand_landmark_list_on?1:0, mand_landmarks_N, mand_landmark_list_N,
      (int)mand_lm_cache.size(), inwin_lm, mand_drawn_lm,
      mand_roots?1:0, mand_root_list_on?1:0, mand_roots_deg, mand_root_list_deg,
      (int)mand_rt_cache.size(), inwin_rt, mand_drawn_rt,
      0, mand_landmark_list_targeted?1:0, mand_root_list_targeted?1:0,
      (int)mand_landmark_list.size(), mand_landmark_selected,
      mand_landmark_list_from_uv?1:0,
      mand_ll.real(), mand_ll.imag(), mand_ur.real(), mand_ur.imag());
  }

  if (aborted) {
    set_status("mandelbrot: stopped -- the picture is incomplete, so nothing was cached");
  } else if (computing) {
    if (mand_trap && !found_TLB) {
      /* Say HOW MUCH more zoom is wanted.  "No trap-like balls in this window" is true but
         unhelpful: the balls have to be common to the whole window, so the window has to be
         very small -- around 1e-8 across, which is about twenty-eight halvings from the
         initial view.  Told only "well zoomed-in", a user reasonably concludes after a few
         zooms that the feature is broken.  It is not: at 4e-9 across the whole pane comes
         out certified. */
      std::stringstream T;
      T.precision(2);
      T << "mandelbrot: done -- no trap-like balls common to this whole window, which is "
        << (mand_ur.real() - mand_ll.real()) << " across; traps want it below about 1e-8"
           " (right-click to zoom, or use Zoom in repeatedly)";
      set_status(T.str());
    }
    else                         set_status("mandelbrot: done");
  }
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
      int col = (std::abs(mand_pixel_group_to_cpx(Point2d<int>(i,j))) > 1.0
                 ? WhitePixel(display, screen)
                 : mand_get_color(mand_data_grid[i][j]));
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
  
  if (window_mode != MANDELBROT) {
    draw_limit();
  }
  recompute_point_data();
  
}


//this zooms on whatever the highlighted ifs is
void IFSGui::mand_zoom(double radius_multiplier) {
  double radius = (mand_ur.real() - mand_ll.real())/2.0;
  cpx c = IFS.z;
  radius *= radius_multiplier;
  /* clamp, so no path into this function can wander outside the unit disc */
  if (radius > MAND_MAX_RADIUS) radius = MAND_MAX_RADIUS;
  mand_ll = c - cpx(radius, radius);
  mand_ur = c + cpx(radius, radius);
  mand_reset_mesh();
  //the targeted landmark list is built for one window, so it goes stale on a zoom
  if (mand_landmarks) mand_rebuild_landmarks();
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
  if (mand_landmarks) mand_rebuild_landmarks();
  draw_mand();
}

void IFSGui::mand_reset_mesh() {
  mand_pixel_width = (mand_ur.real() - mand_ll.real()) / double(W_mand_plot.width);
  mand_pixel_group_width = mand_pixel_group_size * mand_pixel_width;
  mand_num_pixel_groups = W_mand_plot.width / mand_pixel_group_size;
  mand_data_grid.resize(mand_num_pixel_groups);
  for (int i=0; i<mand_num_pixel_groups; ++i) {
    mand_data_grid[i] = std::vector<PointNd<4,int> >(mand_num_pixel_groups, PointNd<4,int>(-1));
  }
  mand_grid_connected_valid = false;
  mand_grid_contains_half_valid = false;
  mand_grid_trap_valid = false;
}




void IFSGui::recompute_point_data() {
  std::stringstream T;
  T.precision(15);
  int difficulty;
  //full precision and plain space-separated, so the line can be selected and
  //pasted elsewhere (e.g. back into W_point_param_entry as "re im") -- a
  //labelled "Location: (a,b)" string could not be reused that way
  T.str("");
  T << IFS.z.real() << " " << IFS.z.imag() << "  " << std::abs(IFS.z)
    << "  " << std::arg(IFS.z)*180.0/M_PI;
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
    std::vector<Bitword> half_words;
    IFS.half_balls(half_words, point_contains_half_depth, 4);
    for (int i=0; i<(int)half_words.size(); ++i) {
      std::cout << half_words[i];
      double epsilon;
      if (IFS.certify_set_B_point(half_words[i], true, epsilon)) {
        std::cout << " (certified all within " << epsilon << ")\n";
      } else {
        std::cout << " (not certified)\n";
      }
    }
    if (half_words.size() > 0) std::cout << "\n";
  }
  W_point_contains_half_status.update_text(T.str());
  
  if (!point_uv_words_check) {
    T.str(""); T << "(disabled)";
  } else {
    /* CAP THE BEAM.  Unbounded (the default) this search quadruples its working set at
       every level, and at a parameter where the coincidence is exact -- a landmark, which
       is exactly where one asks -- it runs the machine out of memory instead of answering.
       Nothing here can interrupt it either: recompute_point_data has no abort.  20000 is
       far more beam than the answer needs at any depth reachable from the stepper, and the
       truncation keeps the closest candidates, so the minimiser is not thrown away. */
    IFS.find_closest_uv_words(point_uv_words, point_uv_words_depth, 0.0, 20000);
    T.str("");
    //find_closest_uv_words gives up and leaves the list empty when the IFS has
    //no minimal enclosing ball, which happens in the corners of the default view
    if (point_uv_words.empty()) {
      T << "invalid (|s| >= 1)";
    } else {
      T << point_uv_words[0].first << " " << point_uv_words[0].second;
    }
  }
  W_point_uv_words_status.update_text(T.str());
  
  if (!point_trap_check) {
    T.str(""); T << "(disabled)";
  } else {
    std::string why;
    T.str("");
    if (!find_trap_words(point_trap_depth, point_trap_words, why)) {
      T << why;
    } else {
      T << point_trap_words[0].first << " " << point_trap_words[0].second;
    }
  }
  W_point_trap_status.update_text(T.str());
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
    /* box_ur.imag(), not box_ur.real().  Comparing the imaginary part of a path point
       against the REAL part of the running upper bound is a copy-paste slip that has been
       here since 2014: whenever a path's imaginary coordinates never exceed its real ones
       -- the common case -- this test never fires, the top edge of the box stays pinned to
       the first point's height, and the trap search below runs on a box that does not
       contain the path. */
    if (path.path[i].imag() > box_ur.imag()) {
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
  ifs temp_IFS;
  double TLB_C, TLB_Z;
  temp_IFS.set_params(av, av);
  if (!temp_IFS.TLB_for_region(TLB, box_ll, box_ur, 15, &TLB_C, &TLB_Z, verbose)) {
    std::cout << "Couldn't find TLB for box " << box_ll << " " << box_ur << " at depth " << mand_trap_depth << "\n";
    return;
  }
  
  if (verbose>0) {
    std::cout << "Got the TLB constants: " << TLB_C << " " << TLB_Z << "\n";
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
      difficulty = temp_IFS.check_TLB(TLB, &TLB_C, &TLB_Z, epsilon, NULL, mand_trap_depth);
      if ( difficulty < 0 ) {
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

  //Try the desired column, and if it overflows the bottom of the window,
  //retry exactly once more in a new column just to the right of everything
  //already stacked in the old one.  "Column" means what it always meant here:
  //the set of widgets sharing this exact ul.x -- at ul.x==0 that is the plot,
  //the whole point panel, and the status row all at once, so the widest of
  //those (not just w2's own desired_x) is where the new column has to start.
  for (int attempt = 0; attempt < 2; ++attempt) {
    //go through and check the other widgets to see how
    //far down they obstruct this one, and how far right this column extends
    int greatest_y_obstruction = 0;
    int greatest_column_x = desired_x;
    for (int i=0; i<(int)widgets.size(); ++i) {
      if (widgets[i] == w1) continue;
      if (widgets[i]->ul.x == desired_x) {
        if (widgets[i]->ul.y + widgets[i]->height > greatest_y_obstruction) {
          greatest_y_obstruction = widgets[i]->ul.y + widgets[i]->height;
          //std::cout << "Found widget " << i << " obstructs to height " << greatest_y_obstruction << "\n";
        }
        if (widgets[i]->ul.x + widgets[i]->width > greatest_column_x) {
          greatest_column_x = widgets[i]->ul.x + widgets[i]->width;
        }
      }
    }
    if (greatest_y_obstruction + w2->height <= main_window_height) {
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
      return;
    }
    //it didn't fit -- on the first attempt, try again one column over;
    //if that column is no wider than this one (nothing to gain by moving),
    //or this was already the second attempt, fall through to the warning.
    if (attempt == 0 && greatest_column_x > desired_x) {
      desired_x = greatest_column_x;
    } else {
      break;
    }
  }

  //It still doesn't fit.  This used to return here silently, with only a
  //commented-out std::cout diagnostic, so a widget could vanish with no clue
  //why and no visible symptom beyond "the control does nothing".  Now: warn
  //on stderr, and give the widget a well-defined off-screen ul.  Without
  //this, ul is left however it was before this call -- zero on the very
  //first pack, or, on a later mode switch, a stale on-screen position from
  //the last time this member packed successfully -- and anything packed
  //relative to it (pack_widget_upper_right(w2, w3)) would silently inherit
  //that stale position instead of also landing off-screen.
  std::cerr << "pack_widget_upper_right: widget (w=" << w2->width << ", h=" << w2->height
            << ") does not fit in column x=" << desired_x << " (main_window_height="
            << main_window_height << ") even in a new column; it will not be shown.\n";
  w2->ul = Point2d<int>(1000000, 1000000);
}


void IFSGui::detach_widget(Widget* w) {
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w) {
      widgets.erase(widgets.begin()+i);
      break;
    }
  }
  //if the widget being pulled out from under the layout was the one
  //receiving keys, drop focus -- otherwise it would keep receiving KeyPress
  //events after its pixmap/ul may no longer mean anything.
  if (w == focus_widget) focus_widget = NULL;
  w->clear();
}

//Shared click_signal for every WidgetEntry: a click just moves the caret to
//the nearest character boundary under the pointer.  By the time this runs,
//the ButtonPress branch in main_loop has already made this widget
//focus_widget (see the hit-test loop), so focus_widget IS the widget that
//was clicked.
void IFSGui::S_entry_click(XEvent* e) {
  if (e->type != ButtonPress || focus_widget == NULL) return;
  WidgetEntry* W = (WidgetEntry*)focus_widget;
  int wx = e->xbutton.x - W->ul.x - W->text_position.x - W->label_px;
  W->caret = clamp_int((wx + W->char_w/2)/W->char_w, 0, (int)W->text.size());
  W->redraw();
}

//Tab from an entry: walk the packed widgets (in packing order, which is a
//reasonable proxy for visual order) for the next one that wants keys, and
//move focus there.  Wraps around; a no-op if focus_widget isn't packed or no
//other entry exists.
void IFSGui::focus_next_entry() {
  int n = (int)widgets.size();
  if (n == 0) return;
  int start = 0;
  for (int i=0; i<n; ++i) {
    if (widgets[i] == focus_widget) { start = i; break; }
  }
  for (int k=1; k<=n; ++k) {
    Widget* w = widgets[(start+k)%n];
    if (w->wants_keys()) {
      if (focus_widget != NULL) {
        ((WidgetEntry*)focus_widget)->focused = false;
        focus_widget->redraw();
      }
      focus_widget = w;
      ((WidgetEntry*)focus_widget)->focused = true;
      focus_widget->redraw();
      return;
    }
  }
}



//Say what we are doing, and make it appear NOW.  Nothing else in this file
//flushes, so without the XFlush a progress message drawn inside a long loop
//would only reach the server when the loop ended and main_loop called
//XNextEvent again -- i.e. exactly when it is no longer of any use.
void IFSGui::set_status(const std::string& s) {
  if (!status_widget_ready) return;
  W_status.update_text(s);
  XFlush(display);
}

/* set_status plus a thin bar along the bottom of the status row.  A long recompute
   reporting "row 403 of 565" is legible only if you read it; a bar is legible at a
   glance, which is the point of a progress indicator. */
void IFSGui::set_progress(double frac, const std::string& t) {
  if (!status_widget_ready) return;
  /* update_text redraws and blits the widget but does not flush; the single XFlush at
     the end of this function covers both it and the bar.  Calling set_status here
     instead would flush twice for one visible change. */
  W_status.update_text(t);
  if (frac < 0.0) frac = 0.0;
  if (frac > 1.0) frac = 1.0;
  Widget& S = W_status;
  int h = 3, y = S.height - h - 1;
  int w = (int)(frac*(S.width - 2));
  XSetForeground(display, S.gc, get_rgb_color(0.80, 0.80, 0.80));
  XFillRectangle(display, S.p, S.gc, 1, y, S.width - 2, h);
  XSetForeground(display, S.gc, get_rgb_color(0.10, 0.45, 0.85));
  if (w > 0) XFillRectangle(display, S.p, S.gc, 1, y, w, h);
  XCopyArea(display, S.p, main_window, S.gc, 0, 0, S.width, S.height, S.ul.x, S.ul.y);
  XFlush(display);
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




/* Where the user actually sees this window, which is the position of the window MANAGER's
   frame and not of our own window: under a reparenting WM ours sits at a small offset
   inside the frame, so recreating it at our own absolute coordinates would walk it down
   and to the right by the border width on every mode switch.  Walk up to the child of the
   root -- that is the frame -- and take its origin.  False if anything is unavailable, in
   which case the caller keeps whatever it had. */
bool IFSGui::window_frame_origin(int& ox, int& oy) {
  if (!main_window_initialized) return false;
  Window the_root = RootWindow(display, screen);
  Window w = main_window;
  for (int guard = 0; guard < 32; ++guard) {
    Window r = 0, parent = 0, *kids = NULL;
    unsigned int nkids = 0;
    if (!XQueryTree(display, w, &r, &parent, &kids, &nkids)) return false;
    if (kids) XFree(kids);
    if (parent == the_root || parent == 0) break;
    w = parent;
  }
  XWindowAttributes a;
  if (!XGetWindowAttributes(display, w, &a)) return false;
  ox = a.x; oy = a.y;
  return true;
}

void IFSGui::reset_and_pack_window() {
  //destroy the main window, if it exists
  if (main_window_initialized) {
    int fx = 0, fy = 0;
    if (window_frame_origin(fx, fy)) { main_window_x = fx; main_window_y = fy; }
    XDestroyWindow(display, main_window);
  }
  //figure out how big to make the window (here ss = sidebar size)
  //for each of limit and mandelbrot, we need height x+ss and width x+ss,
  //if if it's just one of them, we can set the height to be the display height
  //minus 200 (say), as long as x+ss is smaller than the width
  //otherwise, we can make it x+ss, as long as 2x+ss < width
  int display_width = XDisplayWidth(display, screen);
  int display_height = XDisplayHeight(display, screen);
  int ss = (limit_sidebar_size > mand_sidebar_size ? limit_sidebar_size : mand_sidebar_size);
  int width_rest = (window_mode == BOTH ? (display_width-2*ss)/2 :
                                                 display_width - ss );
  //the bottom panel is 9 rows of 20 (7 of point data, a status line, and one
  //spare); it used to be exactly the 7 rows, which meant any new row silently
  //vanished in pack_widget_upper_right.  The height allowance below went up by
  //the same 40, so a window that is limited by the screen height rather than
  //its width still fits on the screen.
  int height_rest = display_height - 210;
  int x = (width_rest > height_rest ? height_rest : width_rest);

  if (window_mode == MANDELBROT) {
    main_window_height = x + 180;
    main_window_width = x + mand_sidebar_size;
  } else if (window_mode == LIMIT) {
    main_window_height = x + 180;
    main_window_width = x + limit_sidebar_size;
  } else {
    main_window_height = x + 180;
    main_window_width = 2*x + mand_sidebar_size + limit_sidebar_size;
  }
  
  //compute the widths
  limit_pixel_width = (limit_ur.real() - limit_ll.real())/double(x);
  
  //create the window
  /* Reuse wherever the user last had it.  Switching window mode destroys and recreates
     this window, and hard-coding (20,20) snapped it back to the corner every time, undoing
     whatever placement they had chosen.  main_window_x/y are updated from ConfigureNotify. */
  main_window = XCreateSimpleWindow(display,
                                    RootWindow(display, screen),
                                    main_window_x, main_window_y,
                                    main_window_width, main_window_height, 4,
                                    BlackPixel(display, screen), WhitePixel(display, screen));
  //the window used to have no name at all, so window lists and screenshot tools
  //had nothing to call it
  XStoreName(display, main_window, "schottky");
  XSetIconName(display, main_window, "schottky");
  //and ask to be told about the close button rather than being killed by it
  wm_delete_window = XInternAtom(display, "WM_DELETE_WINDOW", False);
  XSetWMProtocols(display, main_window, &wm_delete_window, 1);
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
  hover_widget = NULL;   //the old pointer would dangle across a relayout
  //every widget member is about to be reassigned (new ul, new pixmap in some
  //cases), so a focus left pointing at one of them would mean typing goes to
  //a widget whose position or pixmap no longer belongs to this window mode.
  focus_widget = NULL;

  //used to fill the text
  std::stringstream T; 
  
  //the status line spans the whole window, so a long message or filename fits
  W_status = WidgetText(this, "ready", main_window_width, 20);

  //stuff for the IFS computations
  W_switch_to_limit = WidgetButton(this, "Switch to limit", -1, 20, &IFSGui::S_switch_to_limit);
  W_switch_to_mandelbrot = WidgetButton(this, "Switch to mandelbrot", -1,20, &IFSGui::S_switch_to_mandelbrot);
  W_switch_to_combined = WidgetButton(this, "Switch to combined", -1, 20, &IFSGui::S_switch_to_combined);
  
  W_point_title = WidgetText(this, "Current IFS status:", -1, 20);
  //wide enough for "re im |s| arg" at 15 digits of precision each, with room
  //to spare, so nothing is clipped and the line stays pasteable in full
  W_point_point = WidgetText(this, "initializing", 700, 20);
  W_point_param_entry = WidgetEntry(this, "s = ", "", 340, 20, &IFSGui::S_point_param_entered);
    /* Beside it on the same row: a u,v pair, which is the other way of naming a parameter --
       not by its coordinates but by the coincidence it realises.  It shares the row because
       the bottom panel has the width to spare and no rows to spare. */
    W_point_uv_entry = WidgetEntry(this, "uv = ", "", 360, 20, &IFSGui::S_point_uv_entered);
  W_point_connected_check = WidgetCheck(this, "Connected", 105, 20, point_connected_check, &IFSGui::S_point_connected);
  W_point_connected_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_point_connected_decrease_depth);
  T.str(""); T << point_connected_depth;
  W_point_connected_depth_label = WidgetText(this, T.str(), -1, 20);
  W_point_connected_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_point_connected_increase_depth);
  W_point_connected_status = WidgetText(this, "initializing", -1, 20);
  
  W_point_contains_half_check = WidgetCheck(this, "Contains 0", 105, 20, point_contains_half_check, &IFSGui::S_point_contains_half);
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
  
  W_point_trap_check = WidgetCheck(this, "Trap", 105, 20, point_trap_check, &IFSGui::S_point_trap);
  W_point_trap_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_point_trap_decrease_depth);
  T.str(""); T << point_trap_depth;
  W_point_trap_depth_label = WidgetText(this, T.str(), -1, 20);
  W_point_trap_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_point_trap_increase_depth);
  W_point_trap_status = WidgetText(this, "initializing", x, 20);

  //if the limit set is shown:
  if (window_mode != MANDELBROT) {
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
    W_limit_uv_graph = WidgetCheck(this, "uv graph", -1, 20, limit_uv_graph, &IFSGui::S_limit_uv_graph);
    W_limit_uv_graph_depth_title = WidgetText(this, "Depth:", -1, 20);
    W_limit_uv_graph_depth_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_limit_uv_graph_decrease_depth);
    T.str(""); T << limit_uv_graph_depth;
    W_limit_uv_graph_depth_label = WidgetText(this, T.str(), -1, 20);
    W_limit_uv_graph_depth_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_limit_uv_graph_increase_depth);
    W_limit_trap = WidgetCheck(this, "Trap:", -1, 20, limit_trap, &IFSGui::S_limit_trap);
    W_limit_trap_depth_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_limit_trap_decrease_depth);
    T.str(""); T << limit_trap_depth;
    W_limit_trap_depth_label = WidgetText(this, T.str(), 20, 20);
    W_limit_trap_depth_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_limit_trap_increase_depth);
    W_limit_trap_zoom = WidgetButton(this, "to trap", -1, 20, &IFSGui::S_limit_trap_zoom);
    W_limit_nifs = WidgetCheck(this, "nIFS", -1, 20, limit_nifs, &IFSGui::S_limit_nifs);
    //W_limit_gifs = WidgetCheck(this, "gIFS", -1, 20, limit_gifs, &IFSGui::S_limit_gifs);
    W_limit_2d = WidgetCheck(this, "2d IFS", -1, 20, limit_2d, &IFSGui::S_limit_2d);
    
    pack_widget_upper_right(NULL, &W_limit_plot);
    if (window_mode == LIMIT) {
      pack_widget_upper_right(&W_limit_plot, &W_switch_to_mandelbrot);
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
    pack_widget_upper_right(&W_limit_plot, &W_limit_trap);
    pack_widget_upper_right(&W_limit_trap, &W_limit_trap_depth_leftarrow);
    pack_widget_upper_right(&W_limit_trap_depth_leftarrow, &W_limit_trap_depth_label);
    pack_widget_upper_right(&W_limit_trap_depth_label, &W_limit_trap_depth_rightarrow);
    pack_widget_upper_right(&W_limit_plot, &W_limit_trap_zoom);
    pack_widget_upper_right(&W_limit_plot, &W_limit_nifs);
    //pack_widget_upper_right(&W_limit_plot, &W_limit_gifs);
    pack_widget_upper_right(&W_limit_plot, &W_limit_2d);
    
  }
  
  //if the mandelbrot set is shown:
  if (window_mode != LIMIT) {
    W_mand_plot = WidgetDraw(this, x,x, &IFSGui::S_mand_draw);
    W_mand_options_title = WidgetText(this, "mandelbrot options:", -1, 20);
    W_mand_recenter = WidgetButton(this, "Recenter", -1, 20, &IFSGui::S_mand_recenter);
    W_mand_zoom_title = WidgetText(this, "Zoom: ", -1, 20);
    W_mand_zoom_in = WidgetButton(this, "in", -1, 20, &IFSGui::S_mand_zoom_in);
    W_mand_zoom_out = WidgetButton(this, "out", -1, 20, &IFSGui::S_mand_zoom_out);  
    W_mand_mesh_title = WidgetText(this, "Mesh size:", -1, 20);
    W_mand_mesh_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_decrease_mesh);
    std::stringstream T;  T.str("");  T << mand_pixel_group_size;
    W_mand_mesh_label = WidgetText(this, T.str(), -1, 20);
    W_mand_mesh_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_increase_mesh);
    W_mand_connected_check = WidgetCheck(this, "Connected:", 105, 20, (mand_connected ? 1 : 0), &IFSGui::S_mand_connected);
    W_mand_connected_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_connected_decrease_depth);
    T.str("");  T << mand_connected_depth;
    W_mand_connected_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_connected_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_connected_increase_depth);
    W_mand_contains_half_check = WidgetCheck(this, "Contains 0:", 105, 20, (mand_contains_half ? 1 : 0), &IFSGui::S_mand_contains_half);
    W_mand_contains_half_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_contains_half_decrease_depth);
    T.str("");  T << mand_contains_half_depth;
    W_mand_contains_half_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_contains_half_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_contains_half_increase_depth);
    W_mand_trap_check = WidgetCheck(this, "Traps:", 105, 20, (mand_trap ? 1 : 0), &IFSGui::S_mand_trap);
    W_mand_trap_depth_leftarrow = WidgetLeftArrow(this, 20,20, &IFSGui::S_mand_trap_decrease_depth);
    T.str("");  T << mand_trap_depth;
    W_mand_trap_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_trap_depth_rightarrow = WidgetRightArrow(this, 20,20, &IFSGui::S_mand_trap_increase_depth);
    W_mand_scale_check = WidgetCheck(this, "Scale + axes", 105, 20, (mand_scale ? 1 : 0), &IFSGui::S_mand_scale);
    W_mand_circle_half_check = WidgetCheck(this, "|s|=1/2", 105, 20, (mand_circle_half ? 1 : 0), &IFSGui::S_mand_circle_half);
    W_mand_circle_sqrt2_check = WidgetCheck(this, "|s|=1/sqrt2", 105, 20, (mand_circle_sqrt2 ? 1 : 0), &IFSGui::S_mand_circle_sqrt2);
    /* The other reading of the uv box: two FINITE words, whose exact coincidence
       u(0) = v(0) gives roots rather than landmark points.  Laid out exactly like Landmarks
       below it, because it is the same kind of thing: a set of marked parameters with a
       complexity bound. */
    W_mand_roots_check = WidgetCheck(this, "Roots:", 105, 20, (mand_roots ? 1 : 0), &IFSGui::S_mand_roots);
    W_mand_roots_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_roots_decrease);
    { std::stringstream RT; RT << mand_roots_deg;
      W_mand_roots_label = WidgetText(this, RT.str(), 20, 20); }
    W_mand_roots_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_roots_increase);
    W_mand_landmarks_check = WidgetCheck(this, "Landmarks:", 105, 20, (mand_landmarks ? 1 : 0), &IFSGui::S_mand_landmarks);
    W_mand_landmarks_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_landmarks_decrease);
    T.str("");  T << mand_landmarks_N;
    W_mand_landmarks_label = WidgetText(this, T.str(), 20, 20);
    W_mand_landmarks_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_landmarks_increase);
    W_mand_annulus_button = WidgetButton(this, "Limit traps in annulus", -1, 20, &IFSGui::S_mand_annulus);
    /* 165, not 190: mand_sidebar_size is 170, so a 190px readout has its right edge (and
       its border) clipped off the end of the panel in every window mode. */
    W_mand_annulus_result = WidgetText(this, mand_annulus_last, 165, 20);
    W_mand_legend = WidgetText(this, "color:", 190, 20);

    W_mand_mouse_label = WidgetText(this, "Mouse:", 50, 20);
    W_mand_mouse_X = WidgetText(this, "Re: initializing", 200, 20);
    W_mand_mouse_Y = WidgetText(this, "Im: initializing", 200, 20);
    W_mand_output_window = WidgetButton(this, "Write window coords", -1, 20, &IFSGui::S_mand_output_window);
    W_mand_output_picture = WidgetButton(this, "Write picture", -1, 20, &IFSGui::S_mand_output_picture);
    W_mand_output_picture_size_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_output_picture_decrease_size);
    T.str(""); T << mand_output_picture_size;
    W_mand_output_picture_size_label = WidgetText(this, T.str(), -1, 20);
    W_mand_output_picture_size_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_output_picture_increase_size);
    
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
    W_mand_path_movie_with_mandelbrot = WidgetCheck(this, "Movie with mandelbrot", -1, 20, path.movie_with_mandelbrot, &IFSGui::S_mand_path_movie_with_mandelbrot);
    W_mand_path_find_uv_words_button = WidgetButton(this, "Find uv words along path", -1, 20, &IFSGui::S_mand_path_find_uv_words);
    W_mand_path_find_half_words_button = WidgetButton(this, "Find half words along path", -1, 20, &IFSGui::S_mand_path_find_half_words);
    W_mand_path_half_depth_title = WidgetText(this, "Half ball depth", -1, 20);
    W_mand_path_half_depth_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_path_half_decrease_depth);
    T.str(""); T << path.half_depth;
    W_mand_path_half_depth_label = WidgetText(this, T.str(), -1, 20);
    W_mand_path_half_depth_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_path_half_increase_depth);
    W_mand_path_half_start_stop_title = WidgetText(this, "Start/stop:", -1, 20);
    W_mand_path_half_start_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_path_half_decrease_start);
    T.str(""); T << path.half_start;
    W_mand_path_half_start_label = WidgetText(this, T.str(), -1, 20);
    W_mand_path_half_start_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_path_half_increase_start);
    W_mand_path_half_end_leftarrow = WidgetLeftArrow(this, 20, 20, &IFSGui::S_mand_path_half_decrease_end);
    T.str(""); T << path.half_end;
    W_mand_path_half_end_label = WidgetText(this, T.str(), -1, 20);
    W_mand_path_half_end_rightarrow = WidgetRightArrow(this, 20, 20, &IFSGui::S_mand_path_half_increase_end);




    if (window_mode == MANDELBROT) {
      pack_widget_upper_right(NULL, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_combined);
    } else {
      pack_widget_upper_right(&W_limit_center_title, &W_mand_plot);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_limit);
      pack_widget_upper_right(&W_mand_plot, &W_switch_to_mandelbrot);
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
    pack_widget_upper_right(&W_mand_plot, &W_mand_scale_check);
    pack_widget_upper_right(&W_mand_plot, &W_mand_circle_half_check);
    pack_widget_upper_right(&W_mand_plot, &W_mand_circle_sqrt2_check);
    /* `roots` sits directly above Landmarks: both put marks on the parameter plane, and the
       two searches the uv box drives belong next to the control that shows their results.
       Everything below moves down a row, which the packer does by itself. */
    pack_widget_upper_right(&W_mand_plot, &W_mand_roots_check);
    pack_widget_upper_right(&W_mand_roots_check, &W_mand_roots_leftarrow);
    pack_widget_upper_right(&W_mand_roots_leftarrow, &W_mand_roots_label);
    pack_widget_upper_right(&W_mand_roots_label, &W_mand_roots_rightarrow);
    pack_widget_upper_right(&W_mand_plot, &W_mand_landmarks_check);
    pack_widget_upper_right(&W_mand_landmarks_check, &W_mand_landmarks_leftarrow);
    pack_widget_upper_right(&W_mand_landmarks_leftarrow, &W_mand_landmarks_label);
    pack_widget_upper_right(&W_mand_landmarks_label, &W_mand_landmarks_rightarrow);
    pack_widget_upper_right(&W_mand_plot, &W_mand_annulus_button);
    pack_widget_upper_right(&W_mand_plot, &W_mand_annulus_result);
    pack_widget_upper_right(&W_mand_plot, &W_mand_legend);
    pack_widget_upper_right(&W_mand_plot, &W_mand_mouse_label);
    pack_widget_upper_right(&W_mand_plot, &W_mand_mouse_X);
    pack_widget_upper_right(&W_mand_plot, &W_mand_mouse_Y);
    pack_widget_upper_right(&W_mand_plot, &W_mand_output_window);
    pack_widget_upper_right(&W_mand_plot, &W_mand_output_picture);
    pack_widget_upper_right(&W_mand_output_picture, &W_mand_output_picture_size_leftarrow);
    pack_widget_upper_right(&W_mand_output_picture_size_leftarrow, &W_mand_output_picture_size_label);
    pack_widget_upper_right(&W_mand_output_picture_size_label, &W_mand_output_picture_size_rightarrow);
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
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_movie_with_mandelbrot);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_uv_words_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_find_half_words_button);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_half_depth_title);
      pack_widget_upper_right(&W_mand_path_half_depth_title, &W_mand_path_half_depth_leftarrow);
      pack_widget_upper_right(&W_mand_path_half_depth_leftarrow, &W_mand_path_half_depth_label);
      pack_widget_upper_right(&W_mand_path_half_depth_label, &W_mand_path_half_depth_rightarrow);
      pack_widget_upper_right(&W_mand_plot, &W_mand_path_half_start_stop_title);
      pack_widget_upper_right(&W_mand_path_half_start_stop_title, &W_mand_path_half_start_leftarrow);
      pack_widget_upper_right(&W_mand_path_half_start_leftarrow, &W_mand_path_half_start_label);
      pack_widget_upper_right(&W_mand_path_half_start_label, &W_mand_path_half_start_rightarrow);
      pack_widget_upper_right(&W_mand_path_half_start_rightarrow, &W_mand_path_half_end_leftarrow);
      pack_widget_upper_right(&W_mand_path_half_end_leftarrow, &W_mand_path_half_end_label);
      pack_widget_upper_right(&W_mand_path_half_end_label, &W_mand_path_half_end_rightarrow);
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
  pack_widget_upper_right(NULL, &W_point_trap_check);
  pack_widget_upper_right(&W_point_trap_check, &W_point_trap_leftarrow);
  pack_widget_upper_right(&W_point_trap_leftarrow, &W_point_trap_depth_label);
  pack_widget_upper_right(&W_point_trap_depth_label, &W_point_trap_rightarrow);
  pack_widget_upper_right(&W_point_trap_rightarrow, &W_point_trap_status);
  //this is the one spare 20px row the bottom panel's height budget sets
  //aside (see the comment on height_rest above); if a later row is ever
  //added here too, pack_widget_upper_right now wraps it into a new column
  //instead of silently dropping it, but it is worth keeping this panel
  //inside its budget rather than relying on that.
  pack_widget_upper_right(NULL, &W_point_param_entry);
  pack_widget_upper_right(&W_point_param_entry, &W_point_uv_entry);

  //the status line is the last row, under everything else.  If there was no
  //room for it, leave set_status inert rather than blitting to a stale position.
  status_widget_ready = false;
  pack_widget_upper_right(NULL, &W_status);
  status_widget_ready = (!widgets.empty() && widgets.back() == &W_status);


  
  //invalidate the grids and stuff
  if (window_mode != LIMIT) {
    mand_reset_mesh();
  }
  
  //draw all the widgets
  for (int i=0; i<(int)widgets.size(); ++i) {
    widgets[i]->initial_draw();
  }
  
  //plot the limit set


  /* One-line hover help per control (see Widget::help).  Assigned here, after the
     widgets are constructed and before they are packed, so the list stays in one
     place rather than scattered through the construction code. */
  W_switch_to_limit.help = "show only the limit set Lambda_s of the selected parameter";
  W_switch_to_mandelbrot.help = "show only parameter space";
  W_switch_to_combined.help = "show the limit set and parameter space side by side";
  W_limit_plot.help = "the limit set Lambda_s; click to recenter";
  W_limit_depth_title.help = "word length n: the picture is the image of a point (or disk) under all words in f,g of length n";
  W_limit_depth_auto.help = "adjust n automatically so the picture stays sharp while zooming";
  W_limit_chunky.help = "draw disks instead of points; then the picture provably CONTAINS the limit set";
  W_limit_colors.help = "color f(Lambda) and g(Lambda) differently; where they overlap g is drawn on top";
  W_limit_uv_graph.help = "label the balls at the given depth and join closest f/g pairs; paths in this graph are short-hop paths";
  W_limit_trap.help = "draw the trap certifying this parameter: the two words u,v whose clusters interleave, in magenta and green";
  W_limit_trap_depth_leftarrow.help = "shorter trap words: faster, but a trap may exist only at greater length";
  W_limit_trap_depth_rightarrow.help = "longer trap words: finds traps closer to the boundary of M, and costs more";
  W_limit_trap_zoom.help = "take the view to the trap, which is only a few pixels across at this scale";
  W_limit_nifs.help = "plot half the difference set Lambda-Lambda (the 3-generator IFS with translations -1,0,1); s is in M exactly when the marked point 1/s lies in it";
  W_limit_2d.help = "the affine limit set; only meaningful for real part of s between 0.6 and 0.68";
  W_mand_plot.help = "parameter space: left-click picks a parameter (drag to sweep), right-click picks and zooms in";
  W_mand_recenter.help = "recenter the window on the selected parameter";
  /* The zoom buttons had no hover help at all, in either pane -- the one omission left in the
     help list, and worth filling because the mandelbrot zoom is where the landmark and root
     layers change behaviour: below a quarter across they switch to the targeted search. */
  W_mand_zoom_in.help = "halve the width of the parameter window, keeping the selected parameter centred; below about a quarter across the landmark and root layers switch to their targeted searches";
  W_mand_zoom_out.help = "double the width of the parameter window; it will not go past |s| = 1, where there is nothing";
  W_limit_zoom_in.help = "magnify the limit set about the centre of the pane";
  W_limit_zoom_out.help = "shrink the limit set about the centre of the pane";
  W_mand_mesh_title.help = "how many screen pixels make up one computed cell; smaller is sharper and slower";
  W_mand_connected_check.help = "plot M, the set where Lambda_s is connected; the shade shows how hard a touching pair was to find";
  W_mand_contains_half_check.help = "plot M_0, the parameters whose limit set contains the point 0 (the program works internally in a normalization where that point is 1/2, hence contains_half in the source)";
  W_mand_trap_check.help = "search for CKW traps, which certify a whole disk of parameters inside M; needs a well zoomed-in window and works best near the boundary of M";
  W_mand_scale_check.help = "draw a scale bar and the coordinate axes; a deep zoom is otherwise impossible to place";
  W_mand_circle_half_check.help = "draw the circle |s| = 1/2; below it the attractor is a Cantor set, so dM lies outside";
  W_mand_circle_sqrt2_check.help = "draw the circle |s| = 1/sqrt2; at or above it the attractor is robustly connected, so dM lies inside";
  W_mand_landmarks_check.help = "mark the renormalization points of complexity a+b <= N, where the limit-trap mechanism applies; click one to select it exactly";
  W_mand_landmarks_leftarrow.help = "lower the complexity bound N (fewer, simpler landmark points)";
  W_mand_landmarks_rightarrow.help = "raise the complexity bound N (up to 12); above 9 needs a zoomed-in window, where the landmarks are found by a targeted search instead of by enumerating all of them";
  W_mand_annulus_button.help = "at the selected LANDMARK (green -- not a root, whose coincidence is exact rather than asymptotic), certify limit traps over one fundamental annulus of E_sigma and write the figure as png/eps/pdf; Escape aborts";
  W_mand_annulus_result.help = "result of the last annulus run; it stays here after the status line moves on";
  W_mand_legend.help = "priority order of the color layers: the first one with a value wins; white means disconnected, or |s| > 1 where there is no attractor";
  W_mand_output_window.help = "print the current window corners to the terminal";
  W_mand_output_picture.help = "recompute this view at the chosen resolution and write it as png, eps and pdf plus a .txt sidecar";
  W_mand_path_create_by_drawing_button.help = "draw a path in parameter space by clicking a sequence of points";
  W_mand_path_create_by_boundary_button.help = "trace the boundary of the component of the complement of M containing the selected parameter";
  W_point_param_entry.help = "type a parameter and press Enter: 're im', 're+imi', a bare angle in degrees on |s|=1/sqrt2, or 'deg@r'";
  W_point_uv_entry.help = "with brackets, prefix(period) for each, Enter marks the LANDMARK points where u,v agree asymptotically, e.g. fg(fffffggggg) gf(gggggfffff); without brackets, two finite words, Enter marks the ROOTS where u(0)=v(0) exactly, e.g. fgff gfgg";
  W_mand_roots_check.help = "mark every s with u(0) = v(0) for some pair of FINITE words -- the roots of the {0,+-1} polynomials -- up to the degree at right; amber, against the landmarks' green, and click one to select it.  Only Im s > 0 is marked, one per conjugate pair, so the lower half is bare by convention and not for want of roots";
  W_mand_roots_leftarrow.help = "lower the degree bound (fewer roots)";
  W_mand_roots_rightarrow.help = "raise the degree bound; in a zoomed window the search is targeted and cheap, at a wide view it enumerates all 3^d polynomials (degree 10 is 95776 roots and some seconds)";
  W_point_point.help = "the selected parameter: real, imaginary, |s|, arg in degrees -- select and copy it";
  W_point_connected_check.help = "test whether Lambda_s is connected at the selected parameter";
  W_point_contains_half_check.help = "test whether Lambda_s contains the point 0, the fixed point of the symmetry of the pair f(z)=sz-1, g(z)=sz+1";
  W_point_uv_words_check.help = "find the closest pair of words u,v with u(Lambda) meeting v(Lambda)";
  W_point_trap_check.help = "search for a trap at the selected parameter, certifying a disk of parameters around it inside M";
  W_status.help = "status line: progress, results, and this help";
  /* SCHOTTKY_DUMP_WIDGETS: print the packed rectangle of every widget, with its help
     line as a name.  There is no other way to find out where a control ended up --
     the packer places everything relative to everything else -- so a test script
     otherwise has to guess pixel coordinates off a screenshot, which silently rots
     the moment a label changes width.  Three lines, gated on getenv, like
     SCHOTTKY_TIMING in draw_mand. */
  if (getenv("SCHOTTKY_DUMP_WIDGETS"))
    for (int i=0; i<(int)widgets.size(); ++i)
      fprintf(stderr, "[WIDGET] %3d  x %4d y %4d w %4d h %3d  %s\n", i,
              widgets[i]->ul.x, widgets[i]->ul.y, widgets[i]->width,
              widgets[i]->height, widgets[i]->help.c_str());
  if (window_mode != MANDELBROT) draw_limit();
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
      //an entry with focus eats every key, including 'q' -- otherwise there
      //would be no way to type the letter q into a parameter entry
      if (focus_widget != NULL) { focus_widget->handle_key(&e); continue; }
      if(XLookupKeysym(&e.xkey, 0) == XK_q){
        break;
      }

    //if it involves the mouse, we find the appropriate
    //widget to send it off to
    } else if (e.type == ButtonPress || e.type == MotionNotify) {
      /* Hover help.  Its own scan rather than part of the dispatch loop below,
         because that loop skips widgets with no click_signal -- exactly the plain
         labels and readouts a newcomer most wants explained.  Only report when the
         widget under the pointer CHANGES: motion events arrive in floods and
         set_status calls XFlush. */
      if (e.type == MotionNotify) {
        Widget* over = NULL;
        for (int i=0; i<(int)widgets.size(); ++i)
          if (widgets[i]->contains_pixel(e.xbutton.x, e.xbutton.y)) { over = widgets[i]; break; }
        if (over != hover_widget) {
          hover_widget = over;
          if (over != NULL && !over->help.empty()) set_status(over->help);
        }
      }
      for (int i=0; i<(int)widgets.size(); ++i) {
        if (widgets[i]->contains_pixel( e.xbutton.x, e.xbutton.y) &&
            widgets[i]->click_signal != NULL) {
          //focus only ever moves on a real click, never on mere motion --
          //otherwise the pointer passing over an entry while dragging
          //elsewhere would steal keystrokes out from under the user.
          if (e.type == ButtonPress) {
            if (focus_widget != NULL && focus_widget != widgets[i]) {
              ((WidgetEntry*)focus_widget)->focused = false; focus_widget->redraw();
              focus_widget = NULL;
            }
            if (widgets[i]->wants_keys()) {
              focus_widget = widgets[i];
              ((WidgetEntry*)focus_widget)->focused = true;
            }
          }
          (this->*(widgets[i]->click_signal))(&e);
          break;
        }
      }
    
    } else if (e.type == Expose) {
      Point2d<int> expose_ul( e.xexpose.x, e.xexpose.y );
      int ewidth = e.xexpose.width;
      int eheight = e.xexpose.height;
      for (int i=0; i<(int)widgets.size(); ++i) {
        if (widgets[i]->intersects_rectangle(expose_ul, ewidth, eheight)) {
          widgets[i]->redraw();
        }
      }

    //the window manager's close button; quit the same way 'q' does
    } else if (e.type == ClientMessage) {
      if ((Atom)e.xclient.data.l[0] == wm_delete_window) break;
    }
    
  }
}




void IFSGui::launch(IFSWindowMode m, const cpx& c) {
  
  //set the initial window mode
  window_mode = m;
  
  //set the ifs
  IFS.set_params(c,c);
  
  //set the initial settings for limit and mandelbrot
  limit_ll = cpx(-1, -1.5);
  limit_ur = cpx(2, 1.5);
  limit_depth = 12;
  limit_auto_depth = false;
  limit_chunky = true;
  limit_colors = true;
  limit_uv_graph = false;
  limit_uv_graph_depth = 3;
  limit_trap = false;
  /* 16 rather than 12: with the corrected TLB build depth, the first trap on the circle
     between 40 and 90 degrees needs up to 14 letters, so 12 reported "not found" at
     parameters where a trap exists. */
  limit_trap_depth = 16;
  trap_cache_valid = false;
  trap_cache_ok = false;
  trap_cache_depth = -1;
  limit_trap_located = false;
  limit_trap_radius = 0.0;
  limit_nifs = false;
  nifs_saved_ll = cpx(0.0, 0.0);   //written before read (the toggle sets them), but an
  nifs_saved_ur = cpx(0.0, 0.0);   //uninitialised window is a trap for any future caller
  limit_gifs = false;
  limit_2d = false;
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
  /* Annulus-figure defaults.  These MUST live in the constructor: they were briefly
     inside mand_rebuild_landmarks, which only runs if the user toggles Landmarks on,
     so pressing the annulus button first would have used uninitialised values. */
  mand_annulus_res = 360;
  /* Base tail length, before the descent term and the per-period term; see
     S_mand_annulus.  30 rather than 14 because 14 is not enough to close the annulus at
     s0, a PROVED interior point: there the descent costs 13 and the two periods 2, so a
     base of 14 gives cmax 29 and 96 uncovered pixels, while 30 gives cmax 45 and none.
     The cost of the extra levels is mild -- the search prunes, so cmax 45 is about 2.6x
     cmax 29, not (3/2)^16 -- and cmax is capped at 45 by the raster byte anyway. */
  mand_annulus_cmax = 30;
  mand_annulus_ball_depth = 20;
  mand_annulus_last = "annulus: not run yet";
  mand_scale = true;          /* on by default: knowing the scale is never unwanted */
  mand_circle_half = false;
  mand_circle_sqrt2 = false;
  mand_landmarks = false;
  mand_roots = false;
  mand_roots_deg = 8;            /* 7520 roots, about a second; see S_mand_roots_increase */
  mand_root_list_deg = -1;
  mand_root_list_on = false;
  mand_landmarks_N = 7;              /* 2421 points, 0.12s -- see S_mand_landmarks_increase */
  mand_landmark_list_N = -1;
  mand_landmark_selected = -1;
  mand_landmark_list_targeted = false;
  mand_landmark_list_on = false;
  mand_root_list_targeted = false;
  mand_root_list_ll = cpx(0.0, 0.0);
  mand_root_list_ur = cpx(0.0, 0.0);
  mand_lm_truncated = false;
  mand_rt_truncated = false;
  mand_drawn_lm = 0;
  mand_drawn_rt = 0;
  mand_landmark_list_from_uv = false;
  mand_landmark_list_ll = cpx(0.0, 0.0);
  mand_landmark_list_ur = cpx(0.0, 0.0);
  mand_output_picture_size = 1000;
  
  point_connected_check = true;
  point_connected_depth = 18;
  point_is_connected = false;
  point_contains_half_check = false;
  point_contains_half_depth = 18;
  point_is_contains_half = false;
  point_uv_words_check = false;
  point_uv_words_depth = 18;
  point_uv_words.resize(0);
  point_trap_check = false;
  point_trap_depth = 16;
  point_trap_words.resize(0);

  currently_drawing_path = false;

  //no widget exists yet at all, let alone one with focus
  focus_widget = NULL;
  hover_widget = NULL;

  //W_status has no pixmap until reset_and_pack_window builds it, so set_status
  //must be a no-op until then
  status_widget_ready = false;

  //set up the graphics
  display = XOpenDisplay(NULL);
  if (display == NULL) {
    //Without this check every Xlib call below dereferences a null Display* and the program
    //dies with a segfault and no explanation.  Happens over plain ssh, with DISPLAY unset,
    //or on macOS without XQuartz.
    std::cerr << "schottky: cannot open an X display";
    const char* dsp = getenv("DISPLAY");
    if (dsp == NULL || dsp[0] == 0) std::cerr << " (DISPLAY is not set)";
    else                            std::cerr << " (DISPLAY=" << dsp << ")";
    std::cerr << ".\nThe interactive program needs a display.  Under ssh use 'ssh -X', on\n"
              << "macOS install XQuartz, or use the headless tools certify_arc and funddom,\n"
              << "which need no display.\n";
    exit(1);
  }
  screen = DefaultScreen(display);
  if (display == NULL) {
    std::cout << "Failed to open display\n";
    return;
  }
  main_window_initialized = false;
  main_window_x = 20;              //until the WM tells us otherwise
  main_window_y = 20;
  the_gui_font = NULL;
  
  //reset (set) the window
  limit_sidebar_size = 145;   /* was 130, which clipped "Switch to mandelbrot" (137px) */
  mand_sidebar_size = 170;
  reset_and_pack_window();
  
  //go for it
  main_loop();
  
  //clean up at the end
  XCloseDisplay(display);
  
}