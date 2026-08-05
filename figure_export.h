#ifndef __FIGURE_EXPORT__
#define __FIGURE_EXPORT__

#include <string>
#include <vector>

#include "point.h"

/*****************************************************************************
 * figure_export -- write a mathematical figure as PNG, EPS or PDF.
 *
 * The problem this solves.  A picture of a limit set or of parameter space is
 * fundamentally a raster: one color per pixel, computed by a slow recursion.
 * But the things drawn ON TOP of it -- marked points, the uv-graph, a trap, the
 * circle |s| = 1/sqrt2, axes and labels -- are curves, and in a paper they
 * should stay curves, so that they are sharp at any magnification and can be
 * restyled by the person writing the paper.  So a Figure here is a raster
 * PLUS a list of vector overlays, and the EPS and PDF back ends keep the two
 * separate: the raster is embedded as an image, the overlays are emitted as
 * real paths.  PNG, being a pure raster format, rasterises the overlays.
 *
 * Coordinates.  Everything in a Figure -- the raster's extent and every
 * overlay -- is in MATHEMATICAL coordinates: the complex plane, y increasing
 * upwards.  The window is given by (x0,y0) lower-left and (x1,y1) upper-right.
 * Nothing here knows about pixels except Raster, whose row 0 is the BOTTOM row,
 * matching the convention XGraphics presents (X11 flips y; XGraphics flips it
 * back).  This is deliberate: a caller never has to think about y flips.
 *
 * Dependencies: none beyond the C++ standard library.  In particular the PNG
 * and PDF back ends need DEFLATE, which is implemented here (fixed Huffman
 * with LZ77 matching) rather than by linking zlib, so that this file does not
 * add a build dependency to a program whose only one is X11.
 *
 * Typical use:
 *
 *     figexp::Figure F;
 *     F.set_window(ll.real(), ll.imag(), ur.real(), ur.imag());
 *     F.raster.set_size(1000, 1000);
 *     ... F.raster.set_pixel(i, j, r, g, b) ...          // row 0 = bottom
 *     F.add_disk(cpx_center, radius, figexp::Style::stroke(1,0,0));
 *     figexp::Options opt;                               // sensible defaults
 *     opt.format = figexp::Options::PDF;
 *     std::string err;
 *     if (!figexp::write_figure(F, opt, "figure.pdf", &err)) ... err ...
 *
 * or, for the common case of "whatever the extension says":
 *
 *     figexp::write_auto(F, opt, "figure.eps", &err);
 *
 * ***************************************************************************/

namespace figexp {

/*---------------------------------------------------------------- raster ----*/

//an RGB image.  row 0 is the BOTTOM row, as in XGraphics
struct Raster {
  int w, h;
  std::vector<unsigned char> rgb;         //3*w*h, row-major from the bottom

  Raster() { w = h = 0; }
  bool empty() const { return w <= 0 || h <= 0; }
  void set_size(int W, int H, unsigned char fill = 255);

  //no bounds checking in the inline setter: this is called per pixel
  void set_pixel(int i, int j, unsigned char r, unsigned char g, unsigned char b) {
    unsigned char* p = &rgb[3*(j*w + i)];
    p[0] = r; p[1] = g; p[2] = b;
  }
  void set_pixel(int i, int j, long packed_rgb) {
    set_pixel(i, j, (unsigned char)((packed_rgb >> 16) & 0xFF),
                    (unsigned char)((packed_rgb >>  8) & 0xFF),
                    (unsigned char)( packed_rgb        & 0xFF));
  }
  void get_pixel(int i, int j, unsigned char* r, unsigned char* g, unsigned char* b) const {
    const unsigned char* p = &rgb[3*(j*w + i)];
    *r = p[0]; *g = p[1]; *b = p[2];
  }
};

/*----------------------------------------------------------------- style ----*/

struct Style {
  double r, g, b;          //0..1
  double line_width;       //in points; 0 means "use Options::line_width_pt"
  bool filled;             //fill instead of stroke
  bool dashed;

  Style();
  static Style stroke(double R, double G, double B, double lw = 0.0);
  static Style fill(double R, double G, double B);
  static Style dash(double R, double G, double B, double lw = 0.0);
};

/*------------------------------------------------------------ primitives ----*/

//a polyline or polygon in mathematical coordinates
struct Path {
  std::vector<Point2d<double> > pts;
  bool closed;
  Style st;
  Path() { closed = false; }
};

//a circle, given by center and radius in mathematical coordinates
struct Disk {
  Point2d<double> c;
  double r;
  Style st;
  Disk() { r = 0; }
};

//a dot of fixed SCREEN size, for marking a point (radius in points, not maths)
struct Dot {
  Point2d<double> c;
  double r_pt;
  Style st;
  Dot() { r_pt = 2.0; }
};

struct Text {
  Point2d<double> p;
  std::string s;
  double size_pt;
  int anchor;              //0 left, 1 centered, 2 right
  Style st;
  Text() { size_pt = 9.0; anchor = 0; }
};

/*---------------------------------------------------------------- figure ----*/

struct Figure {
  double x0, y0, x1, y1;               //the mathematical window
  Raster raster;                       //may be empty
  std::vector<Path> paths;
  std::vector<Disk> disks;
  std::vector<Dot>  dots;
  std::vector<Text> texts;
  std::string title;                   //goes in the file's metadata, not the image

  Figure();
  void set_window(double X0, double Y0, double X1, double Y1);

  //convenience builders
  void add_path(const std::vector<Point2d<double> >& p, const Style& s, bool closed = false);
  void add_segment(double ax, double ay, double bx, double by, const Style& s);
  void add_disk(double cx, double cy, double r, const Style& s);
  void add_dot(double cx, double cy, const Style& s, double r_pt = 2.0);
  void add_text(double px, double py, const std::string& s, const Style& st,
                double size_pt = 9.0, int anchor = 0);
  //a circle of radius r about the origin, as a smooth vector path
  void add_circle(double cx, double cy, double r, const Style& s, int nseg = 256);

  bool has_overlays() const {
    return !paths.empty() || !disks.empty() || !dots.empty() || !texts.empty();
  }
};

/*--------------------------------------------------------------- options ----*/

struct Options {
  enum Format { PNG, EPS, PDF };
  Format format;

  double width_pt;         //physical width of vector output, in points (72 = 1 inch)
  int    raster_px;        //width in pixels of PNG output, and of the embedded image;
                           //0 means "use the raster's own size"
  bool   vector_overlays;  //EPS/PDF: emit overlays as paths (true) or rasterise (false)
  //PNG only: emit a tRNS chunk making pixels of exactly (bg_r,bg_g,bg_b)
  //transparent.  It is a single color key, not an alpha channel, so it is
  //useful for a figure drawn on a plain background and useless otherwise.
  //EPS and PDF ignore it; PostScript has no transparency at all.
  bool   transparent;
  unsigned char bg_r, bg_g, bg_b;   //background, used where there is no raster
  bool   draw_frame;       //a hairline box around the figure
  bool   draw_axes;        //the real and imaginary axes, if they cross the window
  double line_width_pt;    //default stroke width
  bool   embed_raster;     //EPS/PDF: include the raster at all

  Options();
};

/*---------------------------------------------------------------- output ----*/

//returns false and sets *err on failure.  err may be NULL.
bool write_figure(const Figure& F, const Options& opt, const std::string& path,
                  std::string* err);

//as write_figure, but the format is taken from the file extension
//(.png/.eps/.epsf/.ps/.pdf); an unknown extension is an error.
bool write_auto(const Figure& F, Options opt, const std::string& path, std::string* err);

//".png" -> Options::PNG etc.  returns false if unrecognised
bool format_from_extension(const std::string& path, Options::Format* out);
const char* format_name(Options::Format f);

/*------------------------------------------------------ exposed utilities ----*/

//raw DEFLATE (RFC 1951) and zlib (RFC 1950) streams, so that other code in the
//tree can write PNGs without linking zlib
void deflate_raw(const unsigned char* data, size_t n, std::vector<unsigned char>& out);
void zlib_wrap(const unsigned char* data, size_t n, std::vector<unsigned char>& out);
unsigned long adler32(const unsigned char* data, size_t n);
unsigned long crc32_of(const unsigned char* data, size_t n);

//write a bare PNG from a bottom-up RGB buffer.  used by write_figure, exposed
//because it is independently useful.
bool write_png(const std::string& path, int w, int h,
               const std::vector<unsigned char>& rgb_bottom_up,
               std::string* err);

/* The viridis color ramp, t in [0,1] clamped, out as three bytes.
   Lives here rather than in the interactive program because two separate things color
   these rasters -- the GUI's annulus export and render_funddom.py -- and they used to
   disagree: the C++ was a straight blue-to-orange interpolation while the Python was a
   six-stop ramp, despite a comment claiming they matched.  The stop table below is
   duplicated verbatim in render_funddom.py (it cannot be shared across languages), and
   if you change either table, change both: they are compared sample by sample and are
   expected to agree to within one rounding step.
   Viridis over the older ramp: it is monotone in lightness, so a deeper tail length always
   reads as brighter rather than merely different, and it stays legible in grayscale. */
void viridis(double t, unsigned char rgb[3]);

} //namespace figexp

#endif
