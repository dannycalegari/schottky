// certify_arc.cc
//
// Headless proof-of-concept driver for certifying that points on the circle
// |s| = 1/sqrt(2) lie in the interior of the Barnsley-Harrington Mandelbrot set M,
// using schottky's certified "trap-like ball" (TLB) machinery.
//
// Reuses ifs::TLB_for_region (build region-valid trap witnesses over a small
// parameter box) and ifs::check_TLB (per-parameter trap + certified disk radius).
//
// Convention: math parameter s == code z == w (diagonal). f:x->zx, g:x->z(x-1)+1.
//
// Build (from the schottky/ directory):
//   make certify_arc
// No X11 is needed, neither the headers nor libX11: the shared trap code is
// compiled with -DIFS_NO_GRAPHICS, which drops its debug windows.  Explicitly,
//   clang++ -O3 -DIFS_NO_GRAPHICS -c ifs.cc       -o ifs_nogfx.o
//   clang++ -O3 -DIFS_NO_GRAPHICS -c trap_grid.cc -o trap_grid_nogfx.o
//   clang++ -O3 -DIFS_NO_GRAPHICS -c movie.cc     -o movie_nogfx.o
//   clang++ -O3 -DIFS_NO_GRAPHICS certify_arc.cc ifs_nogfx.o trap_grid_nogfx.o \
//       movie_nogfx.o -lm -o certify_arc
//
// Usage: run ./certify_arc with no arguments for the authoritative list of all
// subcommands and their arguments.
//
// Except for dumptlbmanyz, which takes a parameter directly, all points are taken
// on the circle |s| = 1/sqrt(2).

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <fstream>
#include <ctime>
#include "ifs.h"
#include "rigor.h"

static const double R0 = 0.7071067811865476; // 1/sqrt(2)

struct CertResult {
  bool tlb_ok;       // TLB_for_region succeeded (region witnesses built)
  bool certified;    // check_TLB found a trap
  int  m;            // certifying word length (difficulty) ; -1 if none
  double eps;        // certified parameter-disk radius (ball of traps) around s
  double d;          // box half-width used
  int    tlb_size;   // number of trap-like balls produced
};

// Certify a single parameter s using a box of half-width d around s.
// beam==0 -> exhaustive check_TLB; beam>0 -> guided beam search check_TLB_bestfirst.
CertResult certify_point(ifs& F, cpx s, int n_depth, int uv_depth, double d, int beam) {
  CertResult res; res.tlb_ok=false; res.certified=false; res.m=-1; res.eps=0.0; res.d=d; res.tlb_size=0;

  cpx ll = s - cpx(d,d);
  cpx ur = s + cpx(d,d);

  std::vector<Ball> TLB;
  double C=0.0, Z=0.0;
  F.set_params(s, s);
  bool ok = F.TLB_for_region(TLB, ll, ur, n_depth, &C, &Z, 0);
  res.tlb_ok = ok;
  res.tlb_size = (int)TLB.size();
  if (!ok || TLB.size()==0) return res;

  F.set_params(s, s); // TLB_for_region restores z; make sure z=s for the check
  double trap_radius = 0.0;
  int m;
  if (beam > 0)
    m = F.check_TLB_bestfirst(TLB, &C, &Z, trap_radius,
                              (std::vector<std::pair<Bitword,Bitword> >*)NULL, uv_depth, beam);
  else
    m = F.check_TLB(TLB, &C, &Z, trap_radius,
                    (std::vector<std::pair<Bitword,Bitword> >*)NULL, uv_depth);
  if (m >= 0) {
    res.certified = true;
    res.m = m;
    res.eps = trap_radius;
  }
  return res;
}

// ---- convex-hull regression test -------------------------------------------
// Guards the robust monotone-chain convex_hull (ifs_trap_like.cc) against the
// degenerate inputs that made the old hull's merge step hang.
static double signed_area(const std::vector<cpx>& P) {
  double a=0.0; int n=(int)P.size();
  for (int i=0;i<n;++i){ const cpx& u=P[i]; const cpx& v=P[(i+1)%n];
    a += u.real()*v.imag() - v.real()*u.imag(); }
  return 0.5*a;
}
static double cross3(const cpx&O,const cpx&A,const cpx&B){
  return (A.real()-O.real())*(B.imag()-O.imag())-(A.imag()-O.imag())*(B.real()-O.real());
}
// returns true if every point of X is inside-or-on the CCW hull `ch`
static bool all_points_inside(const std::vector<int>& ch, const std::vector<cpx>& X){
  int h=(int)ch.size(); if(h<3) return true;
  double tol=1e-9;
  for(size_t p=0;p<X.size();++p){
    for(int i=0;i<h;++i){
      const cpx& A=X[ch[i]]; const cpx& B=X[ch[(i+1)%h]];
      if(cross3(A,B,X[p]) < -tol) return false;  // strictly right of a CCW edge => outside
    }
  }
  return true;
}
static int one_hull_case(const char* name, const std::vector<cpx>& X, bool expect_ge3){
  std::vector<int> ch;
  convex_hull(ch, X);               // must TERMINATE (old code hung here)
  std::vector<cpx> hp(ch.size());
  for(size_t i=0;i<ch.size();++i) hp[i]=X[ch[i]];
  double A = ch.size()>=3 ? signed_area(hp) : 0.0;
  bool inside = all_points_inside(ch, X);
  bool ccw = (ch.size()<3) || (A>0);
  bool ok = inside && ccw && (!expect_ge3 || ch.size()>=3);
  std::printf("  [%-18s] hull=%2d  signed_area=%+.4e  CCW=%d  all_inside=%d  => %s\n",
              name,(int)ch.size(),A,(int)ccw,(int)inside, ok?"PASS":"FAIL");
  return ok?0:1;
}
static int run_hulltest(ifs& F){
  int fails=0;
  // (1) the exact input that provably hung: depth-6 (64-pt) cover at e^{i pi/4}/sqrt2
  {
    double th=M_PI/4.0; cpx s(R0*std::cos(th),R0*std::sin(th));
    F.set_params(s,s);
    double min_r; F.minimal_enclosing_radius(min_r);
    Ball init(0.5,(F.z-1.0)/2.0,(1.0-F.w)/2.0,min_r);
    std::vector<Ball> balls; F.depth=6; F.compute_balls_right(balls,init,6);
    std::vector<cpx> C(balls.size());
    for(size_t i=0;i<balls.size();++i) C[i]=balls[i].center;
    fails += one_hull_case("circle depth6 cover", C, true);
  }
  // (2) synthetic degenerate cases
  { std::vector<cpx> X={cpx(0,0),cpx(1,0),cpx(2,0),cpx(3,0),cpx(1.5,0)};
    fails += one_hull_case("all collinear", X, false); }        // -> 2 endpoints
  { std::vector<cpx> X={cpx(0,0),cpx(0,0),cpx(1,0),cpx(1,0),cpx(0,1),cpx(0,1)};
    fails += one_hull_case("duplicates", X, true); }
  { std::vector<cpx> X={cpx(0,0),cpx(1,0),cpx(0,1)};
    fails += one_hull_case("triangle", X, true); }
  { std::vector<cpx> X;
    for(int i=0;i<400;++i){ double a=2*M_PI*i/400.0; X.push_back(cpx(std::cos(a),std::sin(a))); }
    fails += one_hull_case("circle 400pts", X, true); }         // CCW check on a known convex set
  { std::vector<cpx> X={cpx(0.5,0.5)};
    fails += one_hull_case("single point", X, false); }
  std::printf("hulltest: %s (%d failure%s)\n", fails==0?"ALL PASS":"FAILURES", fails, fails==1?"":"s");
  return fails;
}

// Covering can use the SEGMENT-PERSISTENCE certificate (relaxed 4-arc, sound over the
// box) instead of the conservative all-arcs one, via env RIGOR_SEG.  Default: all-arcs.
static int certify_dispatch(const std::string&u,const std::string&v,int nd,double lo,double hi,double&mm){
  static int use_seg = (std::getenv("RIGOR_SEG")!=nullptr)?1:0;
  return use_seg ? rig::certify_box_seg(u,v,nd,lo,hi,mm) : rig::certify_box(u,v,nd,lo,hi,mm);
}


/* Minimum argc for each command, checked centrally before any branch reads argv[2..], so an
   incomplete command line reports what it needs instead of running off the end of argv.
   Each entry is 2 + the number of REQUIRED positional arguments; optional trailing arguments
   with defaults must NOT be counted, or a correct command line gets refused.  Keep this table
   and the usage block below in step with the branches themselves. */
/* Trap words on the command line.  Internally a word is a bit string, '0' = f and
   '1' = g (ifs.h's convention, and what Bitword::str() prints), and every decoder
   here and in rigor.h tests `c == '0'`.  So a word spelled with the LETTERS f and g
   -- which the usage text used to invite -- was read as all-g, and the commands went
   on to report a confident wrong answer instead of an error.  Accept either spelling,
   normalize to bits, and reject anything else. */
static bool read_word(const char* a, std::string& out, const char* what)
{
  out.clear();
  for (const char* p = a; *p; ++p) {
    switch (*p) {
      case '0': case 'f': case 'F': out += '0'; break;
      case '1': case 'g': case 'G': out += '1'; break;
      default:
        std::fprintf(stderr, "%s: '%s' is not a word -- use bits 0/1 (0=f, 1=g) or the "
                     "letters f/g, outermost letter first.\n", what, a);
        return false;
    }
  }
  if (out.empty()) {
    std::fprintf(stderr, "%s: empty word.\n", what);
    return false;
  }
  return true;
}

static int min_argc_for(const std::string& c)
{
  struct { const char* n; int k; } T[] = {
    {"point",6},{"pointz",7},{"sweep",8},{"cover",7},{"rigorcert",7},{"rigorcertseg",7},{"rigorcover",6},
    {"rigorthin",8},{"dumptlb",5},{"dumptlbmany",8},{"dumptlbmanyz",9},{"mix",6},{"mixfind",6},
    {"mixpair",6},{"mixscan",5},{"findtrap",5},{"ckwfind",6},{"hulltest",2},{0,0}
  };
  for (int i = 0; T[i].n; ++i) if (c == T[i].n) return T[i].k;
  return 2;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr,
      "usage: %s <command> [args]\n\n"
      "Arguments in <> are required, those in [] are optional and show their default.\n\n"
      "  point        <deg> <n_depth> <uv_depth> <box_halfwidth> [beam=0]\n"
      "  pointz       <re> <im> <n_depth> <uv_depth> <box_halfwidth> [beam=0]\n"
      "  sweep        <deg_lo> <deg_hi> <steps> <n_depth> <uv_depth> <box_halfwidth> [beam=0]\n"
      "  cover        <deg_lo> <deg_hi> <n_depth> <uv_depth> <d0>\n"
      "                 [beam=3000] [min_gap_deg=5e-4] [prefix=cover]\n"
      "  rigorcert    <deg_lo> <deg_hi> <u> <v> <n_depth>\n"
      "  rigorcertseg <deg_lo> <deg_hi> <u> <v> <n_depth>\n"
      "  rigorcover   <chunk_id> <deg_lo> <deg_hi> <n_depth>\n"
      "                 [dmax=0.005] [dmin=1e-6] [ndmax=n_depth+6] [findnd=12]\n"
      "  rigorthin    <deg_lo> <deg_hi> <wA> <wB> <wC> <wD> [theta0_deg=10]\n"
      "  dumptlb      <deg> <n_depth> <box_halfwidth>\n"
      "  dumptlbmany  <deg> <n_depth> <box_halfwidth> <ngaps> <ntrials> <nradial>\n"
      "  dumptlbmanyz <re> <im> <n_depth> <box_halfwidth> <ngaps> <ntrials> <nradial>\n"
      "  mix          <deg> <n_depth> <uv_depth> <box_halfwidth> [kmax=0]\n"
      "  mixfind      <deg> <max_uv> <n_depth> <kmax>\n"
      "                 [ratio_goal=1.5] [min_uv=4] [max_pixels=1024] [dmult=1.2] [verbose=0]\n"
      "  mixpair      <deg> <u_bits> <v_bits> <n_depth> [max_pixels=512] [dmult=6] [verbose=1]\n"
      "  mixscan      <max_uv> <n_depth> <kmax>\n"
      "                 [ratio_goal=1.5] [min_uv=4] [max_pixels=1024] [dmult=1.2]\n"
      "                 (reads one angle in degrees per line from stdin)\n"
      "  findtrap     <deg> <max_uv> <max_n> [max_pixels=512] [verbose=1]\n"
      "  ckwfind      <deg> <max_uv> <n_depth> <kmax>\n"
      "                 [ratio_goal=1.5] [min_uv=4] [max_pixels=1024] [dmult=1.2] [verbose=0]\n"
      "  hulltest\n\n"
      "Angles are degrees on the circle |s| = 1/sqrt2, which is the arc the paper works on.\n"
      "M itself is two-dimensional, and pointz and dumptlbmanyz take the parameter directly\n"
      "instead, so they reach the rest of it -- including points like the CKW hexahole.\n"
      "A command given too few arguments says what it needs.\n\n"
      "u, v, wA..wD are words, outermost letter first, written either as bits (0 = f,\n"
      "1 = g -- the form every dump and certificate log prints) or as the letters f/g.\n"
      "rigorcert and rigorcertseg VERIFY an explicitly given (u,v) over a box of parameters\n"
      "in interval arithmetic; rigorcover SEARCHES, covering an arc with certified boxes and\n"
      "checkpointing to cov_<chunk_id>.jsonl so it can be resumed.\n", argv[0]);
    return 2;
  }
  std::string cmd = argv[1];
  if (argc < min_argc_for(cmd)) {
    std::fprintf(stderr, "%s: needs %d argument%s, got %d.  Run %s with no arguments for the "
                 "full list.\n", cmd.c_str(), min_argc_for(cmd) - 2,
                 min_argc_for(cmd) == 3 ? "" : "s", argc - 2, argv[0]);
    return 2;
  }
  ifs F;
  F.initialize(cpx(R0,0.0), cpx(R0,0.0), 512, 0);

  if (cmd == "hulltest") {
    return run_hulltest(F);
  }

  if (cmd == "point") {
    double deg = std::atof(argv[2]);
    int n_depth = std::atoi(argv[3]);
    int uv_depth = std::atoi(argv[4]);
    double d = std::atof(argv[5]);
    int beam = (argc > 6) ? std::atoi(argv[6]) : 0;
    double th = deg*M_PI/180.0;
    cpx s(R0*std::cos(th), R0*std::sin(th));
    CertResult r = certify_point(F, s, n_depth, uv_depth, d, beam);
    std::printf("s = %.15f + %.15f i   (|s|=%.6f, arg=%.4f deg)  beam=%d\n",
                s.real(), s.imag(), std::abs(s), deg, beam);
    std::printf("  TLB_ok=%d  tlb_size=%d  certified=%d  m=%d  eps=%.6e  box_halfwidth=%.6e  eps>=d? %d\n",
                (int)r.tlb_ok, r.tlb_size, (int)r.certified, r.m, r.eps, r.d, (int)(r.eps>=r.d));
    return 0;
  }

  if (cmd == "dumptlb") {
    double deg = std::atof(argv[2]); int n_depth = std::atoi(argv[3]); double d = std::atof(argv[4]);
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    std::vector<Ball> TLB; double C=0.0, Z=0.0; F.set_params(s, s);
    bool ok = F.TLB_for_region(TLB, s-cpx(d,d), s+cpx(d,d), n_depth, &C, &Z, 0);
    std::printf("# analytic trap-like balls (T_0) at arg=%.4f ok=%d C=%.6f Z=%.6f count=%d\n",
                deg,(int)ok,C,Z,(int)TLB.size());
    // Record WHICH parameter these balls belong to.  A ball set is only valid at the s it
    // was computed at, and funddom cross-checks this line against its own sigma; without it
    // a mismatched file silently yields a plausible, meaningless coverage number.
    std::printf("# s = %.17g %.17g\n", s.real(), s.imag());
    std::printf("# center_re center_im radius\n");
    for (int i=0;i<(int)TLB.size();++i)
      std::printf("%.9f %.9f %.9f\n", TLB[i].center.real(), TLB[i].center.imag(), TLB[i].radius);
    return 0;
  }

  if (cmd == "dumptlbmany" || cmd == "dumptlbmanyz") {
    // dumptlbmany <deg> <n_depth> <d> <ngaps> <ntrials> <nradial>
    // Same construction as CKW Def 8.2.3 / ifs::trap_like_balls_from_balls, but emitting
    // EVERY admissible ball instead of one best ball from each of the 5 largest hull gaps:
    // for each of the <ngaps> largest gaps of the convex hull of Sigma_n, for <ntrials>
    // points p along the gap segment and <nradial> depths along the inward normal, the
    // ball of radius alpha = dist(q,Sigma_n) - subtraction about q - x_i (i=1,2, the two
    // hull contact points).  The union of these is a much larger certified trap-like set.
    double deg = 0.0; int n_depth; double d; int ngaps, ntr, nrad; cpx s;
    if (cmd == "dumptlbmanyz") {          // s given explicitly (off the circle |s|=R0)
      double re = std::atof(argv[2]), im = std::atof(argv[3]);
      n_depth = std::atoi(argv[4]); d = std::atof(argv[5]);
      ngaps = std::atoi(argv[6]); ntr = std::atoi(argv[7]); nrad = std::atoi(argv[8]);
      s = cpx(re, im); deg = std::arg(s)*180.0/M_PI;
    } else {
      deg = std::atof(argv[2]); n_depth = std::atoi(argv[3]); d = std::atof(argv[4]);
      ngaps = std::atoi(argv[5]); ntr = std::atoi(argv[6]); nrad = std::atoi(argv[7]);
      double th = deg*M_PI/180.0; s = cpx(R0*std::cos(th), R0*std::sin(th));
    }
    F.set_params(s, s);
    cpx ll = s - cpx(d,d), ur = s + cpx(d,d);
    // The construction lives in ifs::trap_like_balls_many so that the GUI can use it too;
    // this subcommand is now just its command-line face.  Balls come back as (center,
    // radius) pairs, already the two translates per admissible ball.
    std::vector<Ball> TLB;
    if (!F.trap_like_balls_many(TLB, ll, ur, n_depth, ngaps, ntr, nrad, false, 0)) {
      std::printf("# no trap-like balls (no minimal radius, or degenerate hull)\n");
      return 1;
    }
    std::printf("# rich trap-like balls arg=%.6f n_depth=%d #balls=%d\n",
                deg, n_depth, (int)TLB.size());
    // see the note in dumptlb: the parameter is part of the data.  dumptlbmanyz reaches
    // points off |s|=1/sqrt2, where the angle alone does not even identify s.
    std::printf("# s = %.17g %.17g\n", F.z.real(), F.z.imag());
    std::printf("# center_re center_im radius\n");
    for (int i=0;i<(int)TLB.size();++i)
      std::printf("%.9f %.9f %.9f\n", TLB[i].center.real(), TLB[i].center.imag(), TLB[i].radius);
    std::fprintf(stderr, "emitted %d balls\n", (int)TLB.size());
    return 0;
  }

  if (cmd == "mix") {
    // validation/exploration: compare check_TLB vs check_TLB_mixed(kmax) at one arg.
    double deg = std::atof(argv[2]);
    int n_depth = std::atoi(argv[3]);
    int uv = std::atoi(argv[4]);
    double d = std::atof(argv[5]);
    int kmax = (argc > 6) ? std::atoi(argv[6]) : 0;
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    std::vector<Ball> TLB; double C=0.0, Z=0.0;
    F.set_params(s, s);
    bool ok = F.TLB_for_region(TLB, s-cpx(d,d), s+cpx(d,d), n_depth, &C, &Z, 0);
    std::printf("arg=%.4f  TLB_ok=%d  tlb_size=%d\n", deg, (int)ok, (int)TLB.size());
    if (!ok || TLB.size()==0) return 0;
    double e0=0.0, e1=0.0;
    F.set_params(s, s);
    int m0 = F.check_TLB(TLB, &C, &Z, e0, (std::vector<std::pair<Bitword,Bitword> >*)NULL, uv);
    F.set_params(s, s);
    int m1 = F.check_TLB_mixed(TLB, &C, &Z, e1, (std::vector<std::pair<Bitword,Bitword> >*)NULL, uv, kmax);
    std::printf("  check_TLB           : m=%d  eps=%.6e\n", m0, e0);
    std::printf("  check_TLB_mixed(k=%d): m=%d  eps=%.6e\n", kmax, m1, e1);
    std::printf("  [k=0 expect match]  m %s ; eps %s\n",
                (m0==m1?"MATCH":"DIFFER"),
                (std::fabs(e0-e1)<1e-9?"MATCH":(e1>=e0?"mixed>=":"mixed<")));
    return 0;
  }

  if (cmd == "findtrap") {
    // Call the ORIGINAL same-length find_trap (sanity check of the geometric
    // interleaving machinery in this headless build).
    //   findtrap <deg> <max_uv> <max_n> [max_pix=512] [verbose=1]
    double deg = std::atof(argv[2]);
    int max_uv = std::atoi(argv[3]);
    int max_n  = std::atoi(argv[4]);
    int max_pix= (argc>5)? std::atoi(argv[5]) : 512;
    int verbose= (argc>6)? std::atoi(argv[6]) : 1;
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    F.set_params(s, s);
    double eps=0, diff=0;
    bool got = F.find_trap(max_uv, max_n, max_pix, 1.0, &eps, &diff, verbose);
    std::printf("%.6f %s eps=%.6g difficulty=%.6g\n", deg, got?"TRAP":"none", eps, diff);
    return got?0:1;
  }

  if (cmd == "mixpair") {
    // Debug: build copies uL, vL from EXPLICIT outermost-first bit strings and
    // run the geometric interleaving check (isolates it from the seed search).
    //   mixpair <deg> <u_bits> <v_bits> <n_depth> [max_pixels=512] [verbose=1]
    double deg = std::atof(argv[2]);
    std::string us, vs;
    if (!read_word(argv[3], us, "mixpair u") || !read_word(argv[4], vs, "mixpair v")) return 2;
    int n_depth = std::atoi(argv[5]);
    int max_pix = (argc > 6) ? std::atoi(argv[6]) : 512;
    double dmult = (argc > 7) ? std::atof(argv[7]) : 6.0;
    int verbose = (argc > 8) ? std::atoi(argv[8]) : 1;
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    F.set_params(s, s);
    double mir; F.minimal_enclosing_radius(mir);
    Ball init(0.5, (s-1.0)/2.0, (1.0-s)/2.0, mir*dmult);
    // build ball for a word: act_on_left(first) then act_on_right(rest)
    Ball zb = F.act_on_left(us[0]-'0', init);
    for (int i=1;i<(int)us.size();++i) zb = F.act_on_right(us[i]-'0', zb);
    Ball wb = F.act_on_left(vs[0]-'0', init);
    for (int i=1;i<(int)vs.size();++i) wb = F.act_on_right(vs[i]-'0', wb);
    int lz=zb.word_len, lw=wb.word_len;
    double rmin=std::min(zb.radius,wb.radius);
    double ratio=std::abs(zb.center-wb.center)/rmin;
    std::printf("u=%s(|%d|,gen%d) v=%s(|%d|,gen%d) ratio=%.4f centers %.5f%+.5fi  %.5f%+.5fi\n",
      us.c_str(),lz,zb.last_gen_index(), vs.c_str(),lw,wb.last_gen_index(),
      ratio, zb.center.real(),zb.center.imag(), wb.center.real(),wb.center.imag());
    int Lfinal = std::max(lz,lw)+n_depth;
    std::vector<Ball> ZB,WB,balls;
    F.compute_balls_right(ZB, zb, Lfinal-lz);
    F.compute_balls_right(WB, wb, Lfinal-lw);
    int nZ=(int)ZB.size();
    balls.swap(ZB); balls.insert(balls.end(),WB.begin(),WB.end());
    std::printf("gridding %d z-balls + %d w-balls (Lfinal=%d)\n", nZ, (int)WB.size(), Lfinal);
    int nzc=0,nwc=0;
    bool topo = F.trap_interleaves_topological(balls, max_pix, &nzc, &nwc);
    double mtd=-1;
    bool robust = F.find_trap_given_balls(balls, max_pix, &mtd, verbose>1?verbose:0);
    std::printf("%.6f  TOPOLOGICAL=%s (zcut=%d wcut=%d)  ROBUST=%s (mtd=%.4g)  u=%s v=%s\n",
       deg, topo?"TRAP":"none", nzc, nwc, robust?"TRAP":"none", mtd, us.c_str(), vs.c_str());
    return topo?0:1;
  }

  if (cmd == "rigorcert") {
    // Interval-arithmetic certificate over an s-box for an explicit (u,v).
    //   rigorcert <deg_lo> <deg_hi> <u> <v> <nd>
    double lo=std::atof(argv[2]), hi=std::atof(argv[3]);
    std::string u, v;
    if (!read_word(argv[4], u, "rigorcert u") || !read_word(argv[5], v, "rigorcert v")) return 2;
    int nd=std::atoi(argv[6]);
    double mm=-1; int rc=rig::certify_box(u,v,nd,lo,hi,mm);
    std::printf("%s min_margin=%.6e arc=[%.5f,%.5f] (%s)\n",
      (rc==1?"CERTIFIED":(rc==0?"not-a-trap":"interval-FAILED")), mm, lo, hi,
      (rc==1?"certified":(rc==0?"float-not-trap":"box-too-wide")));
    return rc==1?0:1;
  }

  if (cmd == "rigorcertseg") {
    // SEGMENT-PERSISTENCE interval certificate over an s-box (relaxed 4-arc criterion).
    //   rigorcertseg <deg_lo> <deg_hi> <u> <v> <nd>
    double lo=std::atof(argv[2]), hi=std::atof(argv[3]);
    std::string u, v;
    if (!read_word(argv[4], u, "rigorcertseg u") || !read_word(argv[5], v, "rigorcertseg v")) return 2;
    int nd=std::atoi(argv[6]);
    double mm=-1; int rc=rig::certify_box_seg(u,v,nd,lo,hi,mm);
    std::printf("%s min_margin=%.6e arc=[%.5f,%.5f] (%s)\n",
      (rc==1?"CERTIFIED":(rc==0?"not-a-trap":"interval-FAILED")), mm, lo, hi,
      (rc==1?"certified":(rc==0?"float-not-trap":"box-too-wide")));
    return rc==1?0:1;
  }

  if (cmd == "rigorthin") {
    // THIN-TRAP interval certificate: verify the two chords (f-chord wA-wB, g-chord wC-wD)
    // cross transversely over the arg-box, with acute crossing angle >= theta0.
    //   rigorthin <deg_lo> <deg_hi> <wA> <wB> <wC> <wD> [theta0_deg=10]
    double lo=std::atof(argv[2]), hi=std::atof(argv[3]);
    std::string wA,wB,wC,wD;
    if (!read_word(argv[4], wA, "rigorthin wA") || !read_word(argv[5], wB, "rigorthin wB") ||
        !read_word(argv[6], wC, "rigorthin wC") || !read_word(argv[7], wD, "rigorthin wD")) return 2;
    double th0=(argc>8)?std::atof(argv[8]):10.0;
    double mm=-1; int rc=rig::certify_thin_trap(wA,wB,wC,wD,lo,hi,th0,mm);
    std::printf("%s min_margin=%.6e arc=[%.5f,%.5f] theta0=%.1f (%s)\n",
      (rc==1?"CERTIFIED":"interval-FAILED"), mm, lo, hi, th0,
      (rc==1?"transverse-crossing-over-box":"not-certified-over-box"));
    return rc==1?0:1;
  }

  if (cmd == "rigorcover") {
    // In-process rigorous interval covering worker.
    //   rigorcover <cid> <lo> <hi> <nd> [dmax=0.005] [dmin=1e-6]
    // Writes cov_<cid>.json (live) and cov_<cid>.jsonl (durable, resumable) in CWD.
    int cid=std::atoi(argv[2]); double LO=std::atof(argv[3]), HI=std::atof(argv[4]); int nd=std::atoi(argv[5]);
    double dmax=(argc>6)?std::atof(argv[6]):0.005;
    double dmin=(argc>7)?std::atof(argv[7]):1e-6;   // box floor (near-tangencies need ~2e-6)
    int ndmax=(argc>8)?std::atoi(argv[8]):nd+6;      // certify-escalation cap (deep pass: 18)
    int findnd=(argc>9)?std::atoi(argv[9]):12;       // find_trap_mixed search depth (deep pass: 16)
    const double SKIP=1e-4;                          // skip past a genuine hard point
    std::string statef="cov_"+std::to_string(cid)+".json";
    std::string ckptf ="cov_"+std::to_string(cid)+".jsonl";
    std::vector<std::pair<double,double> > iv;   // merged certified intervals
    // resume from checkpoint
    { std::ifstream in(ckptf.c_str()); std::string ln; std::vector<std::pair<double,double> > bx;
      while(std::getline(in,ln)){ double a,b; if(std::sscanf(ln.c_str(),"{\"a\":%lf,\"b\":%lf",&a,&b)==2) bx.push_back(std::make_pair(a,b)); }
      std::sort(bx.begin(),bx.end());
      for(size_t k=0;k<bx.size();++k){ if(!iv.empty() && std::fabs(iv.back().second-bx[k].first)<1e-9) iv.back().second=bx[k].second; else iv.push_back(bx[k]); } }
    // resume at the end of the contiguous certified prefix from LO; reset iv to
    // that prefix so merge_push stays ordered.  Re-covering anything past the
    // first gap is acceptable (the durable .jsonl keeps every certified box).
    std::sort(iv.begin(),iv.end());
    double x=LO; { double end=LO;
      for(size_t k=0;k<iv.size();++k){ if(iv[k].first<=end+1e-9) end=std::max(end,iv[k].second); else break; }
      x=std::max(LO,std::min(end,HI));
      iv.clear(); if(x>LO) iv.push_back(std::make_pair(LO,x)); }
    FILE* ck=std::fopen(ckptf.c_str(),"a");
    std::string stuckf="cov_"+std::to_string(cid)+".stuck";
    FILE* sk=std::fopen(stuckf.c_str(),"a");
    double last_good=dmax;   // carry-forward box width (hard bands crawl then grow back)
    long tinyrun=0;          // consecutive fine boxes -- crawl guard vs stalling near a hard point
    long pairs=0, stuck=0; time_t t0=time(NULL); double last_flush=0;
    auto flush=[&](double pos,bool done){ double now=(double)time(NULL);
      if(!done && now-last_flush<0.3) return; last_flush=now;
      double cov=0; for(size_t k=0;k<iv.size();++k) cov+=iv[k].second-iv[k].first;
      FILE* f=std::fopen(statef.c_str(),"w");
      std::fprintf(f,"{\"id\":%d,\"lo\":%.6f,\"hi\":%.6f,\"pos\":%.6f,\"intervals\":[",cid,LO,HI,pos);
      for(size_t k=0;k<iv.size();++k) std::fprintf(f,"%s[%.6f,%.6f]",(k?",":""),iv[k].first,iv[k].second);
      std::fprintf(f,"],\"covered\":%.6f,\"boxes\":%d,\"pairs\":%ld,\"stuck\":%ld,\"elapsed\":%.0f,\"done\":%s}",
        cov,(int)iv.size(),pairs,stuck,now-(double)t0,done?"true":"false");
      std::fclose(f); };
    auto merge_push=[&](double a,double b){ if(!iv.empty() && std::fabs(iv.back().second-a)<1e-9) iv.back().second=b; else iv.push_back(std::make_pair(a,b)); };
    flush(x,false);
    while(x<HI-1e-12){
      double probe=std::min(x+dmin, HI-1e-9);
      double th=probe*M_PI/180.0; cpx sp(R0*std::cos(th),R0*std::sin(th));
      F.set_params(sp,sp);
      std::string u,v; double mtd=-1;
      // search DEEP, and collect several candidate pairs so we can pick the
      // MOST ROBUST one (largest point-certificate margin) -> widest boxes.
      std::vector<std::pair<std::string,std::string> > cands;
      bool got=F.find_trap_mixed(22,findnd,1024,0,1.5,4,1.2,&u,&v,&mtd,0,&cands);
      if(!got){ // no trap pair here -> genuine hard spot; flag & skip
        std::fprintf(sk,"%.9g nopair\n",x); std::fflush(sk); stuck++; x+=SKIP; last_good=dmax; flush(x,false); continue; }
      { // Pick the CHEAPEST-nd pair (keeps per-box cost low) and, among pairs
        // that certify at that nd, the most ROBUST (largest point margin =>
        // widest boxes).  Only escalate nd if no candidate certifies lower.
        double bestm=-1e18; bool found=false;
        for(int ndt=nd; ndt<=std::min(ndmax,12) && !found; ndt+=2){
          for(size_t ci=0; ci<cands.size(); ++ci){ double mm;
            if(certify_dispatch(cands[ci].first,cands[ci].second,ndt,probe,probe,mm)==1){
              found=true; if(mm>bestm){ bestm=mm; u=cands[ci].first; v=cands[ci].second; } } } }
      }
      pairs++; bool used=false;
      while(x<HI-1e-12){
        double start=std::min(dmax, last_good*2.0);  // carry forward (mild)
        bool ok=false; double good_d=0; double mm;
        // primary: base nd, shrink width until it certifies (cheap; the bulk
        // certifies at nd=6, so we never pay for higher nd here)
        for(double d=start; d>=dmin; d*=0.5)
          if(certify_dispatch(u,v,nd,x,std::min(x+d,HI),mm)==1){ ok=true; good_d=d; break; }
        // escalate nd ONLY if base nd fails at every width (deep spot)
        for(int ndt=nd+2; ndt<=ndmax && !ok; ndt+=2)
          for(double d=start; d>=dmin; d*=0.5)
            if(certify_dispatch(u,v,ndt,x,std::min(x+d,HI),mm)==1){ ok=true; good_d=d; break; }
        if(!ok) break;   // this pair fails even at floor -> re-find fresh pair
        double b=std::min(x+good_d,HI);
        std::fprintf(ck,"{\"a\":%.9g,\"b\":%.9g,\"u\":\"%s\",\"v\":\"%s\"}\n",x,b,u.c_str(),v.c_str()); std::fflush(ck);
        merge_push(x,b); x=b; last_good=good_d; used=true; flush(x,false);
        // crawl guard: if we've taken many consecutive ultra-fine boxes we are
        // approaching a hard point (margin -> 0); flag the neighborhood and jump
        // past it rather than grinding indefinitely.
        if(good_d<1e-5){ if(++tinyrun>500){ std::fprintf(sk,"%.9g crawl\n",x); std::fflush(sk); stuck++; x+=2e-3; last_good=dmax; tinyrun=0; break; } }
        else tinyrun=0;
      }
      if(!used){ // pair found but certifies no box even at floor -> near a hard point
        std::fprintf(sk,"%.9g nocert\n",x); std::fflush(sk); stuck++; x+=SKIP; last_good=dmax; flush(x,false); }
    }
    flush(HI,true); std::fclose(ck); std::fclose(sk);
    std::printf("chunk %d done: covered %.5f of [%.3f,%.3f], pairs=%ld stuck=%ld\n",
      cid, [&]{double c=0;for(size_t k=0;k<iv.size();++k)c+=iv[k].second-iv[k].first;return c;}(), LO,HI,pairs,stuck);
    return 0;
  }

  if (cmd == "ckwfind") {
    // Find a trap (topological), then test the CKW round-disk POINT certificate.
    //   ckwfind <deg> <max_uv> <n_depth> <kmax> [ratio_goal=1.5] [min_uv=4]
    //           [max_pixels=1024] [dmult=1.2] [verbose=0]
    // Prints: <deg> <topo TRAP|none> <u> <v> <ckw OK|FAIL> <min_margin>
    double deg  = std::atof(argv[2]);
    int max_uv  = std::atoi(argv[3]);
    int n_depth = std::atoi(argv[4]);
    int kmax    = std::atoi(argv[5]);
    double rg   = (argc > 6) ? std::atof(argv[6]) : 1.5;
    int min_uv  = (argc > 7) ? std::atoi(argv[7]) : 4;
    int max_pix = (argc > 8) ? std::atoi(argv[8]) : 1024;
    double dm   = (argc > 9) ? std::atof(argv[9]) : 1.2;
    int verbose = (argc > 10) ? std::atoi(argv[10]) : 0;
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    F.set_params(s, s);
    std::string u="-", v="-"; double mtd=-1.0;
    bool topo = F.find_trap_mixed(max_uv, n_depth, max_pix, kmax, rg, min_uv, dm,
                                  &u, &v, &mtd, verbose);
    bool ckw=false; double mm=-1.0;
    if (topo) { F.set_params(s,s); ckw = F.ckw_point_certificate(u, v, n_depth, max_pix, dm, &mm, verbose); }
    std::printf("%.6f %s %s %s %s %.6g\n", deg, (topo?"TRAP":"none"),
                (topo?u.c_str():"-"), (topo?v.c_str():"-"),
                (topo?(ckw?"ckwOK":"ckwFAIL"):"-"), mm);
    std::fflush(stdout);
    return (topo && ckw) ? 0 : 1;
  }

  if (cmd == "mixfind") {
    // Geometric mixed-length trap finder at one or many args on |s|=1/sqrt2.
    //   mixfind <deg> <max_uv> <n_depth> <kmax> [ratio_goal=1.5] [min_uv=4]
    //           [max_pixels=1024] [dmult=1.2] [verbose=0]
    // Prints: <deg> <TRAP|none> <u> <v> <|u|> <|v|> <#cut_comps>
    // Exit code 0 if a trap was found, 1 otherwise (handy for scripting).
    double deg   = std::atof(argv[2]);
    int max_uv   = std::atoi(argv[3]);
    int n_depth  = std::atoi(argv[4]);
    int kmax     = std::atoi(argv[5]);
    double rg    = (argc > 6) ? std::atof(argv[6]) : 1.5;
    int min_uv   = (argc > 7) ? std::atoi(argv[7]) : 4;
    int max_pix  = (argc > 8) ? std::atoi(argv[8]) : 1024;
    double dmult = (argc > 9) ? std::atof(argv[9]) : 1.2;
    int verbose  = (argc > 10) ? std::atoi(argv[10]) : 0;
    double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
    F.set_params(s, s);
    std::string u="-", v="-"; double mtd=-1.0;
    bool got = F.find_trap_mixed(max_uv, n_depth, max_pix, kmax, rg, min_uv, dmult,
                                 &u, &v, &mtd, verbose);
    std::printf("%.6f %s %s %s %d %d %.6g\n",
                deg, (got?"TRAP":"none"),
                (got?u.c_str():"-"), (got?v.c_str():"-"),
                (got?(int)u.size():0), (got?(int)v.size():0), mtd);
    std::fflush(stdout);
    return got ? 0 : 1;
  }

  if (cmd == "mixscan") {
    // Batch mixfind over a list of args read from stdin (one deg per line).
    //   mixscan <max_uv> <n_depth> <kmax> [ratio_goal=1.5] [min_uv=4] [max_pixels=1024] [dmult=1.2]
    int max_uv   = std::atoi(argv[2]);
    int n_depth  = std::atoi(argv[3]);
    int kmax     = std::atoi(argv[4]);
    double rg    = (argc > 5) ? std::atof(argv[5]) : 1.5;
    int min_uv   = (argc > 6) ? std::atoi(argv[6]) : 4;
    int max_pix  = (argc > 7) ? std::atoi(argv[7]) : 1024;
    double dmult = (argc > 8) ? std::atof(argv[8]) : 1.2;
    double deg;
    while (std::scanf("%lf", &deg) == 1) {
      double th = deg*M_PI/180.0; cpx s(R0*std::cos(th), R0*std::sin(th));
      F.set_params(s, s);
      std::string u="-", v="-"; double mtd=-1.0;
      bool got = F.find_trap_mixed(max_uv, n_depth, max_pix, kmax, rg, min_uv, dmult,
                                   &u, &v, &mtd, 0);
      std::printf("%.6f %s %s %s %d %d %.6g\n",
                  deg, (got?"TRAP":"none"),
                  (got?u.c_str():"-"), (got?v.c_str():"-"),
                  (got?(int)u.size():0), (got?(int)v.size():0), mtd);
      std::fflush(stdout);
    }
    return 0;
  }

  if (cmd == "pointz") {
    /* `point` off the circle.  M is a two-dimensional set and most of it is nowhere near
       |s| = 1/sqrt2; the angle-only interface was an artefact of the arc the paper needed,
       not of the mathematics.  Same computation as `point`, with s given directly -- as
       dumptlbmanyz already does for the ball dumps. */
    double re = std::atof(argv[2]);
    double im = std::atof(argv[3]);
    int n_depth = std::atoi(argv[4]);
    int uv_depth = std::atoi(argv[5]);
    double d = std::atof(argv[6]);
    int beam = (argc > 7) ? std::atoi(argv[7]) : 0;
    cpx s(re, im);
    if (std::abs(s) >= 1.0) {
      std::fprintf(stderr, "pointz: |s| = %.6f >= 1, so the maps are not contractions and "
                   "there is no attractor.\n", std::abs(s));
      return 2;
    }
    CertResult r = certify_point(F, s, n_depth, uv_depth, d, beam);
    std::printf("s = %.15f + %.15f i   (|s|=%.6f, arg=%.4f deg)  beam=%d\n",
                s.real(), s.imag(), std::abs(s),
                std::atan2(s.imag(), s.real())*180.0/M_PI, beam);
    std::printf("  TLB_ok=%d  tlb_size=%d  certified=%d  m=%d  eps=%.6e  "
                "box_halfwidth=%.6e  eps>=d? %d\n",
                r.tlb_ok?1:0, r.tlb_size, r.certified?1:0, r.m, r.eps, r.d,
                (r.certified && r.eps >= r.d) ? 1 : 0);
    return r.certified ? 0 : 1;
  }

  if (cmd == "sweep") {
    double lo = std::atof(argv[2]);
    double hi = std::atof(argv[3]);
    int steps = std::atoi(argv[4]);
    int n_depth = std::atoi(argv[5]);
    int uv_depth = std::atoi(argv[6]);
    double d = std::atof(argv[7]);
    int beam = (argc > 8) ? std::atoi(argv[8]) : 0;
    std::printf("# sweep |s|=1/sqrt2  n_depth=%d uv_depth=%d box_halfwidth=%.3e beam=%d\n", n_depth, uv_depth, d, beam);
    std::printf("# %8s %6s %5s %4s %12s %12s %8s\n",
                "arg_deg","tlbOK","cert","m","eps","min(eps,d)","tlbN");
    for (int i=0;i<=steps;++i) {
      double deg = lo + (hi-lo)*double(i)/double(steps>0?steps:1);
      double th = deg*M_PI/180.0;
      cpx s(R0*std::cos(th), R0*std::sin(th));
      CertResult r = certify_point(F, s, n_depth, uv_depth, d, beam);
      double eff = (r.eps < d ? r.eps : d);
      std::printf("  %8.3f %6d %5d %4d %12.4e %12.4e %8d\n",
                  deg, (int)r.tlb_ok, (int)r.certified, r.m, r.eps, (r.certified?eff:0.0), r.tlb_size);
      std::fflush(stdout);
    }
    return 0;
  }

  if (cmd == "cover") {
    // Adaptive-d covering of an arc of |s|=1/sqrt2 by a union of certified
    // parameter-disks (balls of traps).  One TLB is built per box and reused
    // for many balls (amortizing the expensive TLB build).  A certified disk
    // B_rho(s0) covers exactly the circle points with |s - s0| < rho, i.e. the
    // arc-interval [theta0 - h, theta0 + h] with h = 2*asin(rho/(2 R0)).
    // rho = min(eps, chebyshev distance from s0 to the box boundary), so the
    // covered disk is validly inside the box the TLB was certified on.
    double lo = std::atof(argv[2]);        // deg
    double hi = std::atof(argv[3]);        // deg
    int n_depth = std::atoi(argv[4]);
    int uv_depth = std::atoi(argv[5]);
    double d0 = std::atof(argv[6]);        // starting box half-width
    int beam = (argc > 7) ? std::atoi(argv[7]) : 3000;
    /* Unlike point/sweep, cover ALWAYS calls check_TLB_bestfirst, whose beam
       truncation is gated on beam_width > 0.  So beam=0 here does not mean
       "exhaustive": it means no pruning at all, and the frontier grows without
       bound.  Refuse rather than appear to hang. */
    if (beam <= 0) {
      std::fprintf(stderr, "cover: beam must be positive (cover has no exhaustive mode; "
                   "beam=0 disables pruning entirely and will not finish).  Try 3000.\n");
      return 2;
    }
    double min_gap_deg = (argc > 8) ? std::atof(argv[8]) : 5e-4; // step when stuck
    std::string prefix = (argc > 9) ? argv[9] : "cover";        // output file prefix
    const double DEG = M_PI/180.0;
    double frontier = lo*DEG, hi_r = hi*DEG, min_gap = min_gap_deg*DEG;

    std::vector<std::pair<double,double> > gaps;   // uncertified sub-intervals (rad)
    long balls=0, boxes=0, checks=0;
    std::string certname = prefix + "_cert.log", gapname = prefix + "_gaps.log";
    FILE* cert = std::fopen(certname.c_str(),"w");
    std::fprintf(cert,"# arg_deg Re_s Im_s m eps d_box rho halfwidth_deg u_word v_word   "
                      "(trap (u,v) [0=f,1=g] certifies the disk B_eps(s), i.e. arc [arg-hw,arg+hw])\n");

    // current reusable box
    bool have_box=false; std::vector<Ball> TLB; double C=0,Z=0,dbox=0; cpx s_c=0;
    int fail_run=0;                                // consecutive certification failures
    double eps_last=0.0;                           // last accepted eps (for box refresh)
    const double QUANTUM = 0.01*DEG;               // stuck-detector progress unit
    const long STUCK_BUDGET = 150;                 // checks allowed without a quantum of progress
    double mark = frontier; long checks_at_mark = 0; int stuck_run = 0;
    auto ch_dist = [](cpx s, cpx c, double d){        // box room around s (>0 inside)
      return d - std::max(std::fabs(s.real()-c.real()), std::fabs(s.imag()-c.imag())); };

    std::printf("# cover [%.4f,%.4f] deg  n=%d uv=%d d0=%.2e beam=%d min_gap=%.1e deg\n",
                lo,hi,n_depth,uv_depth,d0,beam,min_gap_deg);
    long prog=0;
    while (frontier < hi_r) {
      if (boxes/50 > prog) { prog=boxes/50;
        std::fprintf(stderr,"  ..progress frontier=%.4f deg  boxes=%ld balls=%ld gaps=%zu\n",
                     frontier/DEG, boxes, balls, gaps.size()); }
      cpx s = cpx(R0,0.0)*std::exp(cpx(0.0,frontier));
      // (re)build the box if we have none, or the frontier nears its edge, or the
      // box room is about to limit rho below the trap radius eps (so placed balls
      // keep rho=eps -- avoids a tiny-ball tail and seam gaps at box boundaries).
      bool need = !have_box || ch_dist(s, s_c, dbox) <= 0.12*dbox
                  || (eps_last>0.0 && ch_dist(s, s_c, dbox) < 3.0*eps_last);
      if (need) {
        // fixed-d box: build the TLB at d0 (shrink only if the build itself fails).
        // Big d0 amortizes the TLB build over many balls; hard-notch args simply
        // fail check_TLB below and are skipped -- to be recovered by a later pass
        // with a smaller d0 (the complementary-arc recursion).
        have_box=false;
        for (double d=d0; d>=d0/30.0; d/=3.0) {
          double th_c = frontier + 0.55*d/R0;                 // center box ahead
          cpx c = cpx(R0,0.0)*std::exp(cpx(0.0,th_c));
          std::vector<Ball> T; double cc=0,zz=0;
          F.set_params(c,c);
          if (F.TLB_for_region(T, c-cpx(d,d), c+cpx(d,d), n_depth, &cc, &zz, 0) && T.size()>0) {
            TLB.swap(T); C=cc; Z=zz; dbox=d; s_c=c; have_box=true; ++boxes; break;
          }
        }
        if (!have_box) {  // cannot even build a TLB here: backoff-skip
          double skip = min_gap; for (int k=0;k<fail_run && skip<0.05*DEG;++k) skip*=2.0;
          gaps.push_back(std::make_pair(frontier,frontier+skip)); frontier += skip; ++fail_run;
          continue;
        }
        s = cpx(R0,0.0)*std::exp(cpx(0.0,frontier));
      }
      double bd = ch_dist(s, s_c, dbox);
      if (bd <= 0) { have_box=false; continue; }
      F.set_params(s,s);
      double eps=0.0;
      std::vector<std::pair<Bitword,Bitword> > tw;    // captures the trap words (u,v)
      int m = F.check_TLB_bestfirst(TLB, &C, &Z, eps, &tw, uv_depth, beam); ++checks;
      double rho = (m>=0 && eps>0.0) ? std::min(eps, bd) : 0.0;
      double h   = (rho>0.0) ? 2.0*std::asin(std::min(1.0, rho/(2.0*R0))) : 0.0;
      if (m>=0 && eps>0.0 && h > 0.0) {
        // trap words for this ball (0=f, 1=g); the trap (u,v) certifies the whole
        // parameter-disk, hence the arc-range [arg-h, arg+h].
        std::string uw = tw.empty()? std::string("?") : tw[0].first.str();
        std::string vw = tw.empty()? std::string("?") : tw[0].second.str();
        std::fprintf(cert,"%.9f %.12f %.12f %d %.6e %.3e %.6e %.6e %s %s\n",
                     frontier/DEG, s.real(), s.imag(), m, eps, dbox, rho, h/DEG,
                     uw.c_str(), vw.c_str());
        ++balls; frontier += h; fail_run=0; eps_last=eps;     // contiguous overlap
      } else {
        // outright certification failure: exponential-backoff skip so a hard
        // region costs O(log width).  Gaps re-attacked by a second cover pass.
        double skip = min_gap; for (int k=0;k<fail_run && skip<0.05*DEG;++k) skip*=2.0;
        gaps.push_back(std::make_pair(frontier,frontier+skip)); frontier+=skip; ++fail_run;
      }
      // Stuck detector: a coverable-but-hard band (e.g. around 45 deg) is covered
      // only by microscopic balls, grinding forever.  If the frontier fails to
      // advance one QUANTUM within a check budget, declare it a hard band, record
      // a gap, and skip past it (escalating the skip on consecutive stalls).
      if (frontier >= mark + QUANTUM) { mark=frontier; checks_at_mark=checks; if(stuck_run>0) stuck_run--; }
      else if (checks - checks_at_mark > STUCK_BUDGET) {
        double sk = QUANTUM; for (int k=0;k<stuck_run && sk<0.06*DEG;++k) sk*=2.0;
        gaps.push_back(std::make_pair(frontier,frontier+sk)); frontier+=sk; ++stuck_run;
        mark=frontier; checks_at_mark=checks; have_box=false; fail_run=0;
      }
    }
    std::fclose(cert);

    // merge gaps (already in increasing order) and report
    std::vector<std::pair<double,double> > mg;
    for (size_t i=0;i<gaps.size();++i) {
      if (!mg.empty() && gaps[i].first <= mg.back().second + 1e-12) {
        if (gaps[i].second > mg.back().second) mg.back().second = gaps[i].second;
      } else mg.push_back(gaps[i]);
    }
    double total = (hi_r - lo*DEG), gapsum=0.0;
    for (size_t i=0;i<mg.size();++i) {              // clamp gaps to [lo,hi]
      if (mg[i].first < lo*DEG) mg[i].first = lo*DEG;
      if (mg[i].second > hi_r)  mg[i].second = hi_r;
      if (mg[i].second > mg[i].first) gapsum += mg[i].second-mg[i].first;
    }
    std::printf("boxes=%ld  balls=%ld  checks=%ld\n", boxes, balls, checks);
    std::printf("covered %.6f of %.6f deg  (%.4f%%);  %zu gap-interval(s), total %.6f deg\n",
                (total-gapsum)/DEG, total/DEG, 100.0*(total-gapsum)/total, mg.size(), gapsum/DEG);
    for (size_t i=0;i<mg.size() && i<40;++i)
      std::printf("   GAP [%.5f, %.5f] deg  (%.5f deg)\n",
                  mg[i].first/DEG, mg[i].second/DEG, (mg[i].second-mg[i].first)/DEG);
    FILE* gf = std::fopen(gapname.c_str(),"w");
    std::fprintf(gf,"# gap_lo_deg  gap_hi_deg\n");
    for (size_t i=0;i<mg.size();++i)
      if (mg[i].second>mg[i].first) std::fprintf(gf,"%.6f %.6f\n", mg[i].first/DEG, mg[i].second/DEG);
    std::fclose(gf);
    std::printf("certificate log: %s ; gaps: %s\n", certname.c_str(), gapname.c_str());
    return 0;
  }

  std::fprintf(stderr, "unknown command '%s'\n", cmd.c_str());
  return 2;
}
