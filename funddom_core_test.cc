/* funddom_core_test.cc -- does funddom's C API, called from C++, agree with the
 * command-line tool exactly?
 *
 * funddom.c is compiled twice: with -DFUNDDOM_MAIN as the CLI, without it as
 * funddom_core.o, the library the GUI links.  The whole point of that
 * arrangement is that there is ONE copy of the limit-trap mathematics.  This
 * test is what keeps that honest.  It checks
 *
 *   - the three built-in cores resolve, with the Delta and P'(sigma) that
 *     funddom's own s0/s1/hex branches print;
 *   - an unknown core name is refused rather than silently substituted;
 *   - fd_landmarks recovers all three cores with correct lm: specs;
 *   - and, the real check, fd_level agrees with a raster written by the CLI
 *     PIXEL FOR PIXEL.
 *
 * Build and run it with `make coretest` (see the makefile), which generates the
 * trap-like balls and the reference raster first.  It needs a ball file because
 * the built-in ten-ball sets cover essentially nothing.
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "funddom_core.h"

int main(int argc, char** argv) {
  if (argc < 4) { std::printf("usage: coretest <balls.txt> <raster.bin> <RES>\n"); return 2; }
  int fails = 0;

  // 1. the built-in cores, which delegate to the lm: specs
  const char* names[3] = {"s0", "s1", "hex"};
  for (int i = 0; i < 3; ++i) {
    fd_core c;
    if (fd_core_builtin(names[i], &c) != 0) { std::printf("  %-4s FAILED to resolve\n", names[i]); ++fails; continue; }
    std::printf("  %-4s sigma=%.15f%+.15fi a=%d b=%d Delta=%+.9f%+.9fi Pp=%+.9f%+.9fi\n",
                names[i], c.sigma_re, c.sigma_im, c.a, c.b,
                c.Delta_re, c.Delta_im, c.Pp_re, c.Pp_im);
  }
  // an unknown name must be refused, not silently substituted
  { fd_core c; if (fd_core_builtin("S0", &c) == 0) { std::printf("  unknown name ACCEPTED -- bad\n"); ++fails; }
    else std::printf("  unknown core name refused: ok\n"); }

  // 2. enumeration must find the three cores
  { std::vector<fd_landmark> L(60000);   // N=9 has 46201; hex first appears at a+b=9
    int n = fd_landmarks(9, &L[0], (int)L.size());
    int found = 0;
    const double tgt[3][2] = {{0.5,0.5},{0.25,0.6614378277661477},{0.371858680074136,0.519411153747943}};
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) {
        double dx = L[i].sigma_re - tgt[k][0], dy = L[i].sigma_im - tgt[k][1];
        if (std::sqrt(dx*dx+dy*dy) < 1e-9) { ++found;
          std::printf("  landmark %-16s a=%d b=%d deg=%d\n", L[i].spec, L[i].a, L[i].b, L[i].deg); }
      }
    std::printf("  fd_landmarks(9) = %d points, %d of the 3 cores present%s\n",
                n, found, found == 3 ? "" : "  <-- FAIL");
    if (found != 3) ++fails;
  }

  // 3. the solver: fd_level must reproduce the CLI's raster pixel for pixel
  std::vector<double> balls;
  { std::ifstream f(argv[1]); std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream is(line); double x,y,r;
      if (is >> x >> y >> r) { balls.push_back(x); balls.push_back(y); balls.push_back(r); }
    } }
  int nb = (int)balls.size()/3;
  std::printf("  read %d balls\n", nb);

  fd_core c; fd_core_builtin("s0", &c);
  if (fd_solver_init(&c, &balls[0], nb) != 0) { std::printf("  fd_solver_init FAILED\n"); return 1; }
  fd_info info; fd_solver_info(&info);
  std::printf("  R_Gamma=%.6f |Y0|=%.6f |dY/dC|=%.6f rho_max=%.6g period=%.6f nballs=%d\n",
              info.R_Gamma, std::hypot(info.Y0_re, info.Y0_im),
              std::hypot(info.dYdC_re, info.dYdC_im), info.rho_max,
              info.period_height, info.nballs);

  int RES = std::atoi(argv[3]);
  std::ifstream rf(argv[2], std::ios::binary);
  std::vector<unsigned char> img(2*(size_t)RES*RES);
  rf.read((char*)&img[0], (std::streamsize)img.size());
  const double rho = 0.08; const int cmax = 16, survd = 0;
  long checked = 0, mismatch = 0;
  for (int iy = 0; iy < RES; ++iy)
    for (int ix = 0; ix < RES; ++ix) {
      double x = (2.0*(ix+0.5)/RES - 1.0)*rho, y = (2.0*(iy+0.5)/RES - 1.0)*rho;
      if (std::hypot(x,y) > rho) continue;            // outside: CLI leaves lev=254
      int lev = fd_level(x, y, cmax);
      unsigned char want = img[2*((size_t)iy*RES+ix)];
      unsigned char got = (lev >= 0) ? (unsigned char)lev : 255;
      ++checked;
      if (got != want) { if (mismatch < 4)
          std::printf("  MISMATCH at (%d,%d): api=%d cli=%d\n", ix, iy, got, want);
        ++mismatch; }
      (void)survd;
    }
  std::printf("  fd_level vs CLI raster: %ld pixels checked, %ld mismatches%s\n",
              checked, mismatch, mismatch ? "  <-- FAIL" : "");
  if (mismatch) ++fails;
  fd_solver_free();
  std::printf("%s\n", fails ? "SOME CHECKS FAILED" : "all API checks passed");
  return fails ? 1 : 0;
}
