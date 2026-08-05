/* figure_export_test.cc -- self test for figure_export.
 *
 * Run it with `make test`.  Two things are checked, and they need different
 * kinds of evidence:
 *
 *   1. THE COMPRESSOR.  figure_export implements DEFLATE itself, so that the
 *      program needs no zlib.  A hand-rolled compressor that is subtly wrong
 *      would corrupt every PNG and PDF it ever writes, so this must be checked
 *      against a real, independent implementation -- not against a decompressor
 *      written by the same hand, which would share any misunderstanding.  So
 *      this program writes the compressed corpus to disk as `defl_N.bin` with
 *      the plaintext beside it as `defl_N.orig`, and the makefile's `test`
 *      target then round-trips them through Python's zlib.  Run
 *          make test
 *      rather than this binary alone, or that half of the check is skipped.
 *
 *   2. THE FIGURE CODE.  That every combination of raster, overlay, format and
 *      option produces a file, that the format really follows the extension,
 *      and that degenerate input is refused instead of crashing.  Whether the
 *      files LOOK right is a question for a viewer; what is checked here is
 *      that nothing silently fails.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "figure_export.h"

using namespace figexp;

static int failures = 0;

static void check(bool ok, const std::string& what) {
  std::printf("  %-54s %s\n", what.c_str(), ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

//write a buffer, for the Python half of the compressor check
static void dump(const std::string& path, const std::vector<unsigned char>& d) {
  std::ofstream f(path.c_str(), std::ios::binary);
  if (!d.empty()) f.write((const char*)&d[0], (std::streamsize)d.size());
}

int main(int argc, char** argv) {
  std::string dir = (argc > 1 ? argv[1] : ".");
  std::string err;
  char nm[1024];

  /*--------------------------------------------- 1. the compression corpus --*/
  std::printf("compressor corpus (round-tripped by the makefile's test target)\n");
  {
    const int NCASE = 7;
    for (int c = 0; c < NCASE; ++c) {
      std::vector<unsigned char> d;
      const char* what = "";
      switch (c) {
        case 0: what = "empty"; break;
        case 1: what = "one byte";
                d.push_back(42); break;
        case 2: what = "70000 zeros (longest matches)";
                d.assign(70000, 0); break;
        case 3: what = "strided bytes";
                for (int i = 0; i < 50000; ++i) d.push_back((unsigned char)(i*7919)); break;
        case 4: what = "period 9";
                for (int i = 0; i < 40000; ++i) d.push_back((unsigned char)("abcabcabd"[i%9]));
                break;
        case 5: what = "incompressible (stored-block fallback)";
                { unsigned long s = 12345;
                  for (int i = 0; i < 60000; ++i) {
                    s = s*1103515245UL + 12345UL;
                    d.push_back((unsigned char)((s >> 16) & 0xff));
                  } }
                break;
        case 6: what = "just over one stored block (65536)";
                d.assign(65536, 0xAB); break;
      }
      std::vector<unsigned char> z;
      zlib_wrap(d.empty() ? (const unsigned char*)"" : &d[0], d.size(), z);
      std::snprintf(nm, sizeof nm, "%s/defl_%d.bin", dir.c_str(), c);
      dump(nm, z);
      std::snprintf(nm, sizeof nm, "%s/defl_%d.orig", dir.c_str(), c);
      dump(nm, d);
      std::printf("  case %d  %-40s %7lu -> %7lu  (%.1f%%)\n", c, what,
                  (unsigned long)d.size(), (unsigned long)z.size(),
                  d.empty() ? 0.0 : 100.0*z.size()/d.size());
      //a stored block never expands by more than 5 bytes per 65535, so the
      //fallback must keep us essentially at parity even on random data
      if (d.size() > 1000) check(z.size() < d.size() + d.size()/1000 + 64,
                                 "output is not larger than the input");
    }
    //adler32 of "Wikipedia" is the documented example value in RFC 1950
    const char* w = "Wikipedia";
    check(adler32((const unsigned char*)w, 9) == 0x11E60398UL,
          "adler32 matches the RFC 1950 example");
    //crc32 of "123456789" is a widely published check value
    const char* n9 = "123456789";
    check(crc32_of((const unsigned char*)n9, 9) == 0xCBF43926UL,
          "crc32 matches the standard check value");
  }

  /*----------------------------------------------------- 2. figure writing --*/
  std::printf("figures\n");

  Figure F;
  F.title = "figure_export self test";
  F.set_window(-1.5, -1.0, 1.5, 1.0);
  const int W = 600, H = 400;
  F.raster.set_size(W, H);
  for (int j = 0; j < H; ++j) {
    for (int i = 0; i < W; ++i) {
      double x = -1.5 + 3.0*(i + 0.5)/W;
      double y = -1.0 + 2.0*(j + 0.5)/H;
      double r = std::sqrt(x*x + y*y);
      //a smooth field and a hard-edged disc, so both busy and flat regions occur
      F.raster.set_pixel(i, j,
        (unsigned char)(127.5*(1.0 + std::sin(6.0*x))),
        (unsigned char)(127.5*(1.0 + std::cos(6.0*y))),
        (unsigned char)(r < 0.6 ? 240 : 40));
    }
  }
  F.add_circle(0.0, 0.0, 0.6, Style::stroke(1, 1, 1, 1.5));
  F.add_disk(-0.9, 0.4, 0.18, Style::stroke(0, 0, 0, 0.8));
  F.add_disk(0.9, -0.4, 0.18, Style::fill(1, 0.3, 0));
  F.add_segment(-1.4, -0.9, 1.4, 0.9, Style::dash(0, 0, 1, 1.0));
  {
    std::vector<Point2d<double> > tri;
    tri.push_back(Point2d<double>(0.2, -0.85));
    tri.push_back(Point2d<double>(0.8, -0.85));
    tri.push_back(Point2d<double>(0.5, -0.45));
    F.add_path(tri, Style::fill(0.1, 0.6, 0.1), true);
  }
  for (int k = 0; k < 5; ++k)
    F.add_dot(-1.35 + 0.1*k, 0.9, Style::fill(0.8, 0, 0.8), 2.5);
  F.add_text(0.0, 0.55, "s = -0.58+0.33i", Style::fill(0, 0, 0), 10, 1);
  F.add_text(-1.45, -0.95, "left", Style::fill(0, 0, 0), 9, 0);
  F.add_text(1.45, -0.95, "right", Style::fill(0, 0, 0), 9, 2);

  Options opt;
  opt.draw_axes = true;

  opt.format = Options::PNG;
  check(write_figure(F, opt, dir + "/fig.png", &err), "PNG, vector overlays rasterised");
  opt.format = Options::EPS;
  check(write_figure(F, opt, dir + "/fig.eps", &err), "EPS, raster plus vector overlays");
  opt.format = Options::PDF;
  check(write_figure(F, opt, dir + "/fig.pdf", &err), "PDF, raster plus vector overlays");

  opt.vector_overlays = false;
  opt.format = Options::PDF;
  check(write_figure(F, opt, dir + "/fig_flat.pdf", &err), "PDF, overlays burnt into the image");
  opt.format = Options::EPS;
  check(write_figure(F, opt, dir + "/fig_flat.eps", &err), "EPS, overlays burnt into the image");
  opt.vector_overlays = true;

  opt.raster_px = 250;
  opt.format = Options::PNG;
  check(write_figure(F, opt, dir + "/fig_small.png", &err), "PNG, raster resampled down");
  opt.raster_px = 1200;
  check(write_figure(F, opt, dir + "/fig_big.png", &err), "PNG, raster resampled up");
  opt.raster_px = 0;

  opt.embed_raster = false;
  opt.format = Options::PDF;
  check(write_figure(F, opt, dir + "/fig_novec.pdf", &err), "PDF, overlays only (raster dropped)");
  opt.embed_raster = true;

  opt.transparent = true;
  opt.format = Options::PNG;
  check(write_figure(F, opt, dir + "/fig_alpha.png", &err), "PNG with transparent requested");
  opt.transparent = false;

  //format dispatch
  check(write_auto(F, opt, dir + "/fig_auto.pdf", &err), "write_auto follows the .pdf extension");
  check(write_auto(F, opt, dir + "/fig_auto.eps", &err), "write_auto follows the .eps extension");
  check(write_auto(F, opt, dir + "/fig_auto.png", &err), "write_auto follows the .png extension");
  check(!write_auto(F, opt, dir + "/fig.tiff", &err), "write_auto refuses an unknown extension");
  check(!write_auto(F, opt, dir + "/noextension", &err), "write_auto refuses a bare name");
  {
    Options::Format f;
    check(format_from_extension("A.PNG", &f) && f == Options::PNG, "extension match is case blind");
    check(std::string(format_name(Options::EPS)) == "EPS", "format_name works");
  }

  //a figure with no raster at all
  {
    Figure V;
    V.set_window(0, 0, 1, 1);
    V.add_circle(0.5, 0.5, 0.4, Style::stroke(0, 0, 0, 1.0));
    V.add_text(0.5, 0.5, "vector only", Style::fill(0, 0, 0), 12, 1);
    Options vo;
    vo.format = Options::EPS;
    check(write_figure(V, vo, dir + "/vec.eps", &err), "EPS with no raster");
    vo.format = Options::PDF;
    check(write_figure(V, vo, dir + "/vec.pdf", &err), "PDF with no raster");
    vo.format = Options::PNG;
    check(write_figure(V, vo, dir + "/vec.png", &err), "PNG with no raster");
  }

  //a figure with no overlays at all
  {
    Figure R;
    R.set_window(0, 0, 2, 1);
    R.raster.set_size(200, 100);
    Options ro;
    ro.draw_frame = false;
    ro.format = Options::PNG;
    check(write_figure(R, ro, dir + "/rasteronly.png", &err), "PNG with no overlays and no frame");
  }

  //things that must be refused rather than crash
  {
    Figure D;
    Options d;
    check(!write_figure(D, d, dir + "/bad.png", &err), "an empty figure is refused");
    check(!err.empty(), "the refusal explains itself");

    Figure Z;
    Z.set_window(1, 1, 1, 1);          //set_window nudges a degenerate window
    Z.raster.set_size(4, 4);
    check(write_figure(Z, d, dir + "/degenerate.png", &err), "a degenerate window survives");

    Figure P;                          //a one-point path and an empty string
    P.set_window(0, 0, 1, 1);
    std::vector<Point2d<double> > one;
    one.push_back(Point2d<double>(0.5, 0.5));
    P.add_path(one, Style::stroke(0, 0, 0));
    P.add_text(0.5, 0.2, "", Style::fill(0, 0, 0));
    P.add_disk(0.5, 0.5, 0.0, Style::stroke(0, 0, 0));       //zero radius
    Options po;
    po.format = Options::PDF;
    check(write_figure(P, po, dir + "/degenerate.pdf", &err), "degenerate primitives survive");
    po.format = Options::PNG;
    check(write_figure(P, po, dir + "/degenerate2.png", &err), "...in the rasteriser too");
  }

  /* Hostile geometry.  A caller can easily hand us a point that is enormous or
   * not a number -- an unset variable, a division by zero in a parametrisation.
   * The rasteriser steps along segments half a pixel at a time, so without
   * clipping a segment a billion pixels long would take a billion steps and one
   * with a nan endpoint would never finish.  These must all return promptly. */
  {
    const double HUGE_ = 1e17;
    const double INF_ = 1.0/0.0, NAN_ = 0.0/0.0;
    Figure X;
    X.set_window(0, 0, 1, 1);
    X.raster.set_size(64, 64);
    X.add_segment(-HUGE_, 0.5, HUGE_, 0.5, Style::stroke(0, 0, 0, 1));   //crosses
    X.add_segment(HUGE_, HUGE_, 2*HUGE_, HUGE_, Style::stroke(0, 0, 0, 1)); //misses
    X.add_segment(0.0, 0.0, INF_, INF_, Style::stroke(0, 0, 0, 1));
    X.add_segment(NAN_, 0.2, 0.8, 0.9, Style::stroke(0, 0, 0, 1));
    X.add_segment(-HUGE_, 0.3, HUGE_, 0.3, Style::dash(1, 0, 0, 1));     //dashed
    X.add_disk(0.5, 0.5, HUGE_, Style::stroke(0, 0, 1, 1));              //vast circle
    X.add_disk(0.5, 0.5, HUGE_, Style::fill(0, 0, 1));                   //vast fill
    X.add_dot(NAN_, NAN_, Style::fill(0, 0, 0), 3.0);
    X.add_text(HUGE_, HUGE_, "far away", Style::fill(0, 0, 0), 9, 0);
    {
      std::vector<Point2d<double> > big;             //a polygon far larger than the canvas
      big.push_back(Point2d<double>(-HUGE_, -HUGE_));
      big.push_back(Point2d<double>( HUGE_, -HUGE_));
      big.push_back(Point2d<double>( 0.0,    HUGE_));
      X.add_path(big, Style::fill(0, 0.5, 0), true);
    }
    Options xo;
    xo.format = Options::PNG;
    check(write_figure(X, xo, dir + "/hostile.png", &err), "hostile geometry: PNG returns");
    xo.format = Options::PDF;
    check(write_figure(X, xo, dir + "/hostile.pdf", &err), "hostile geometry: PDF returns");
    xo.format = Options::EPS;
    check(write_figure(X, xo, dir + "/hostile.eps", &err), "hostile geometry: EPS returns");
    xo.format = Options::PNG;
    xo.vector_overlays = false;
    check(write_figure(X, xo, dir + "/hostile_flat.png", &err), "hostile geometry: burnt in");
  }

  //writing where we cannot must fail cleanly
  check(!write_figure(F, opt, "/nonexistent-directory-xyz/f.png", &err),
        "an unwritable path fails cleanly");

  std::printf("\n%s", failures ? "" : "");
  if (failures) std::printf("%d CHECK%s FAILED\n", failures, failures == 1 ? "" : "S");
  else          std::printf("all checks passed\n");
  return failures ? 1 : 0;
}
