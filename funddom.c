/* funddom.c -- limit-trap parameters C at a coincidence core sigma, and how much of a
 * fundamental domain of the elliptic curve E_sigma = C^* / <sigma^b> they cover.
 *
 * Setup (OUR normalization: f(x)=sx-1, g(x)=sx+1, base point 0, p_w(s)=w(s,0)).
 * sigma is a renormalization point for u,v (length a) and s,t (length b):
 *   pi(u s^inf, sigma) = pi(v t^inf, sigma),  P(z) = p_{u s^inf}(z) - p_{v t^inf}(z),  P(sigma)=0.
 * CKW Lemma 9.2.5: if, for words x,y of length c,
 *   V(C;x,y) = sigma^{-a-c}(p_u-p_v) + sigma^{-c}(p_x-p_y) + sigma^{-a-c} C P'(sigma)
 * is trap-like for sigma, then u s^n x, v t^n y is an ordinary trap for sigma + C sigma^{bn}
 * for all large n; C "admits a limit trap".
 *
 * REDUCTION USED HERE.  Put Delta = p_u-p_v and
 *      Y(C) = -(Delta + C P'(sigma))/sigma^a          (so Y(0) = the coincidence difference).
 * Since p_x-p_y = sum_{j<c} sigma^j e_j with e_j = x_j-y_j in {0,+2,-2}, we get
 *      V(C;x,y) in R   <=>   Y_c := (Y - sum_{j<c} sigma^j e_j)/sigma^c  in  R,
 * where R is the trap-like set at sigma (R = -R, CKW Lemma 8.2.5).  So C admits a
 * *certified* limit trap iff the backward orbit tree of Y(C) under the three inverse
 * branches of the DIFFERENCE-SET IFS {z -> sigma z, sigma z +- 2} reaches R.  The
 * reachable set is  S = union_{c>=0} Phi^c(R),  Phi = Hutchinson operator of that IFS.
 * Y_j must stay in |Y_j| <= 2/(1-|sigma|) (= the radius of the difference set Gamma),
 * which prunes the 3-ary tree to ~(3/2)^j nodes.
 *
 * Exact self-similarity: C -> sigma^b C is Y-Y0 -> sigma^b (Y-Y0), realized on words by
 * prepending (s,t); so C certified at level c => sigma^b C certified at level c+b.  Hence
 * the true set of limit-trap parameters is invariant under <sigma^b> and descends to
 * E_sigma = C^* / <sigma^b>; covering one fundamental domain gives ALL C, i.e. M contains a
 * punctured neighborhood of sigma, i.e. sigma in int(M).
 *
 * usage:  run with no arguments for the authoritative list; in outline,
 *         funddom <core> <ann|log> <cmax> <RESX> <RESY> <rho> <survdepth> <out.bin>
 *                 [nper] [balls.txt]
 *         where <core> is s0, s1, hex, lm:<A>:<B>[:index], or coin:<u>:<v>[:k].
 *   ann : pixels are C in the square [-rho, rho]^2
 *   log : x = theta = arg C in [0,2pi); y = xi = log|C| from log rho (top) down through
 *         nper fundamental periods, each of height b L, L = log(1/|sigma|).  One period
 *         is a fundamental domain of E_sigma (the "annulus" one), glued by
 *         (theta, xi) ~ (theta + b arg sigma, xi - b L).
 * out.bin holds 2 bytes/pixel: (level, gamma) with level = least c certifying C (255 = none),
 * gamma = 1 if Y(C) survives to depth survdepth (numerical over-approximation of Y in Gamma,
 * i.e. of C in T_sigma = the asymptotic limit-trap region of CKW Thm 9.2.2).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "funddom_core.h"

typedef double complex cx;

#define MAXM 64   /* longest coincidence word accepted by the coin: selector */

static cx  SIG, ISIG;
static double RG, RG2, IR, II, RBALLS2;
static int  NB;
static cx  *BC;
static double *BR, *BCR, *BCI, *BR2;
static int  BEST;
static cx  Y0_S, SCALEC_S;   /* Y(C) = Y0_S + SCALEC_S*C; set by fd_solver_init */

/* uniform-grid index over the trap-like balls (they are many when read from file) */
static double GX0, GY0, GH;
static int GNX, GNY, *GSTART, *GLIST;

static void build_grid(void)
{
  double xlo = 1e30, xhi = -1e30, ylo = 1e30, yhi = -1e30, mr = 0;
  for (int i = 0; i < NB; ++i) {
    if (BCR[i] - BR[i] < xlo) xlo = BCR[i] - BR[i];
    if (BCR[i] + BR[i] > xhi) xhi = BCR[i] + BR[i];
    if (BCI[i] - BR[i] < ylo) ylo = BCI[i] - BR[i];
    if (BCI[i] + BR[i] > yhi) yhi = BCI[i] + BR[i];
    if (BR[i] > mr) mr = BR[i];
  }
  GH = mr / 8.0; if (GH <= 0) GH = 1e-3;
  GX0 = xlo - GH; GY0 = ylo - GH;
  GNX = (int)((xhi - xlo) / GH) + 3; GNY = (int)((yhi - ylo) / GH) + 3;
  if ((double)GNX * GNY > 4e7) { GH *= 4; GNX = GNX/4 + 3; GNY = GNY/4 + 3; }
  long NC = (long)GNX * GNY, tot = 0;
  GSTART = calloc(NC + 1, sizeof(int));
  int *cur = NULL;
  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 1) {
      for (long c = 1; c <= NC; ++c) GSTART[c] += GSTART[c - 1];
      tot = GSTART[NC];
      GLIST = malloc(sizeof(int) * (size_t)tot);
      cur = malloc(sizeof(int) * (size_t)NC);
      for (long c = 0; c < NC; ++c) cur[c] = GSTART[c];
    }
    for (int i = 0; i < NB; ++i) {
      int i0 = (int)((BCR[i] - BR[i] - GX0) / GH), i1 = (int)((BCR[i] + BR[i] - GX0) / GH);
      int j0 = (int)((BCI[i] - BR[i] - GY0) / GH), j1 = (int)((BCI[i] + BR[i] - GY0) / GH);
      if (i0 < 0) i0 = 0; if (j0 < 0) j0 = 0;
      if (i1 >= GNX) i1 = GNX - 1; if (j1 >= GNY) j1 = GNY - 1;
      for (int jj = j0; jj <= j1; ++jj)
        for (int ii = i0; ii <= i1; ++ii) {
          long c = (long)jj * GNX + ii;
          if (pass == 0) ++GSTART[c + 1];
          else GLIST[cur[c]++] = i;
        }
    }
  }
  free(cur);
  fprintf(stderr, "  grid %dx%d h=%.5f entries=%ld (%.1f per ball)\n",
          GNX, GNY, GH, tot, (double)tot / NB);
}

static inline int in_R(double yr, double yi)
{
  int ii = (int)((yr - GX0) / GH), jj = (int)((yi - GY0) / GH);
  if (ii < 0 || jj < 0 || ii >= GNX || jj >= GNY) return 0;
  long c = (long)jj * GNX + ii;
  for (int k = GSTART[c]; k < GSTART[c + 1]; ++k) {
    int i = GLIST[k];
    double dr = yr - BCR[i], di = yi - BCI[i];
    if (dr * dr + di * di < BR2[i]) return 1;
  }
  return 0;
}

/* least c <= cmax with (Y - sum_{j<c} sigma^j e_j)/sigma^c in R, e_j in {0,+-2} */
static void dfs(double yr, double yi, int j, int cmax)
{
  if (j >= BEST) return;
  double q = yr * yr + yi * yi;
  if (q > RG2) return;
  if (q < RBALLS2 && in_R(yr, yi)) { BEST = j; return; }
  if (j == cmax) return;
  dfs(yr * IR - yi * II, yr * II + yi * IR, j + 1, cmax);
  double ar = yr - 2.0, br_ = yr + 2.0;
  dfs(ar * IR - yi * II, ar * II + yi * IR, j + 1, cmax);
  dfs(br_ * IR - yi * II, br_ * II + yi * IR, j + 1, cmax);
}

/* Y within |sigma|^D R_Gamma of the level-D difference set: over-approximates Y in Gamma */
static int surv(double yr, double yi, int j, int D)
{
  if (yr * yr + yi * yi > RG2) return 0;
  if (j == D) return 1;
  if (surv(yr * IR - yi * II, yr * II + yi * IR, j + 1, D)) return 1;
  double ar = yr - 2.0, br_ = yr + 2.0;
  if (surv(ar * IR - yi * II, ar * II + yi * IR, j + 1, D)) return 1;
  if (surv(br_ * IR - yi * II, br_ * II + yi * IR, j + 1, D)) return 1;
  return 0;
}


/* ------------------------------------------------- polynomial roots (real coeffs) --
 * Durand-Kerner followed by a Newton polish.  Used by the `coin:` core selector to
 * solve for sigma from a coincidence, so the caller need not know sigma in advance.  */
static cx polyval_(const double *c, int D, cx z)
{ cx r = c[D]; for (int j = D - 1; j >= 0; --j) r = r * z + c[j]; return r; }
static cx polyder_(const double *c, int D, cx z)
{ cx r = D * c[D]; for (int j = D - 1; j >= 1; --j) r = r * z + j * c[j]; return r; }

static void poly_roots_(const double *c, int D, cx *root)
{
  for (int i = 0; i < D; ++i) root[i] = cpow(0.4 + 0.9 * I, i + 1);
  for (int it = 0; it < 500; ++it) {
    double mv = 0;
    for (int i = 0; i < D; ++i) {
      cx num = polyval_(c, D, root[i]) / c[D], den = 1.0;
      for (int j = 0; j < D; ++j) if (j != i) den *= (root[i] - root[j]);
      cx dz = num / den; root[i] -= dz;
      if (cabs(dz) > mv) mv = cabs(dz);
    }
    if (mv < 1e-15) break;
  }
  for (int i = 0; i < D; ++i)
    for (int it = 0; it < 60; ++it) {
      cx d = polyder_(c, D, root[i]);
      if (cabs(d) < 1e-300) break;
      cx dz = polyval_(c, D, root[i]) / d;
      root[i] -= dz;
      if (cabs(dz) < 1e-17) break;
    }
}

/* Set up the renormalization data from a FINITE coincidence u(0) = v(0).
 *
 * With eps_j = -1 for f and +1 for g and d_j = (eps^u_j - eps^v_j)/2 in {0,+-1}, the
 * condition u(0) = v(0) is exactly sum_j d_j sigma^j = 0, so sigma is a root of a
 * polynomial with coefficients in {0,+-1}.  Then u and v are the SAME affine map, so
 * taking s = t (any word) gives pi(u s^inf) = pi(v t^inf): sigma is a renormalization
 * point for (u, v, f, f), with a = |u| and b = 1.  In that case
 *      P(z) = p_u(z) - p_v(z) = 2 sum_j d_j z^j,   P(sigma) = 0,
 *      Delta = p_u - p_v = 0,   P'(sigma) = 2 sum_j j d_j sigma^{j-1}.
 * Returns the number of admissible roots (1/2 < |z| < 1, Im z > 0), fills them in rts,
 * and reports d and its degree.  Returns -1 if (u,v) is not a legal coincidence pair.  */
static int coincidence_data(const char *us, const char *vs, int *dout, int *degout,
                            cx *rts, int maxr)
{
  int m = (int)strlen(us);
  if (m < 2 || (int)strlen(vs) != m) return -1;
  if (m > MAXM) return -2;                  /* d[], c[] and all[] are MAXM on the stack */
  int d[MAXM];
  for (int j = 0; j < m; ++j) {
    if ((us[j] != 'f' && us[j] != 'g') || (vs[j] != 'f' && vs[j] != 'g')) return -1;
    d[j] = ((us[j] == 'g' ? 1 : -1) - (vs[j] == 'g' ? 1 : -1)) / 2;
  }
  if (d[0] == 0) return -1;                 /* u, v must start with different letters */
  if (d[0] == 1) for (int j = 0; j < m; ++j) d[j] = -d[j];   /* normalize d_0 = -1     */
  int D = m - 1;
  while (D > 0 && d[D] == 0) --D;           /* strip a common suffix                  */
  if (D < 1) return -1;
  double c[MAXM];
  for (int j = 0; j <= D; ++j) c[j] = d[j];
  cx all[MAXM];
  poly_roots_(c, D, all);
  int n = 0;
  for (int i = 0; i < D && n < maxr; ++i) {
    double ab = cabs(all[i]);
    if (ab <= 0.5 || ab >= 1.0) continue;
    if (cimag(all[i]) <= 1e-9) continue;    /* conjugates give mirror-image pictures  */
    if (cabs(polyval_(c, D, all[i])) > 1e-10) continue;
    rts[n++] = all[i];
  }
  for (int j = 0; j < m; ++j) dout[j] = d[j];
  *degout = D;
  return n;
}


#ifdef FUNDDOM_MAIN

static void usage(const char *prog)
{
  fprintf(stderr,
"%s -- limit-trap parameters C at a coincidence core, and how much of a fundamental\n"
"domain of the elliptic curve E_sigma = C^* / <sigma^b> they cover.\n"
"\n"
"usage:\n"
"  %s <core> <mode> <cmax> <RESX> <RESY> <rho> <survdepth> <out.bin> [nper] [balls.txt]\n"
"\n"
"  core       s0   sigma = 1/2 + i/2                  (twindragon, a=2, b=1)\n"
"             coin:<u>:<v>[:k]  derive everything from a finite coincidence u(0)=v(0):\n"
"                  give u and v as equal-length f/g words starting with different\n"
"                  letters; sigma is solved for as a root of sum_j d_j z^j with\n"
"                  d_j=(eps^u_j-eps^v_j)/2, and a=|u|, b=1, s=t=f.  If the polynomial\n"
"                  has several admissible roots, k picks one.  Run '%s coin:<u>:<v>'\n"
"                  with no further arguments to list them.  Requires a ball file.\n"
"                  The three named cores above are INFINITE coincidences and cannot be\n"
"                  reached this way -- use s0/s1/hex for those.\n"
"             s1   sigma = 1/4 + (sqrt7/4) i          (tame twindragon, a=1, b=4)\n"
"             hex  the CKW hexahole ~0.37186+0.51941i (a=8, b=1; on the boundary of M,\n"
"                  so a useful control -- it must NOT reach full coverage).\n"
"                  hex requires a ball file.\n"
"             lm:<A>:<B>[:index]  any LANDMARK POINT, given by the sign strings of A\n"
"                  and B over {-,0,+} (meaning coefficients -1, 0, +1) in\n"
"                  Q(z) = A(z)(1-z^b) + z^a B(z), with A starting '-'; a=|A|, b=|B|.\n"
"                  This is exactly the 'spec' column that '%s landmarks' prints, so a\n"
"                  point found there can be handed straight back here.  index picks\n"
"                  among the admissible roots.  Run '%s lm:<A>:<B>' with no further\n"
"                  arguments to see Q and its roots.  Requires a ball file.\n"
"                  Unlike coin:, this reaches the infinite coincidences: s0 is\n"
"                  lm:-+:-:0, s1 is lm:-:+--+:0, hex is lm:-+---+++:-:1.\n"
"  mode       ann  pixels are C in the square [-rho,rho]^2\n"
"             log  x = arg C in [0,2pi), y = log|C| downwards from log rho through nper\n"
"                  fundamental periods, each of height b*log(1/|sigma|).  One period is a\n"
"                  fundamental domain of E_sigma.\n"
"  cmax       largest tail length c to try (cost grows like (3/2)^cmax)\n"
"  RESX RESY  raster size\n"
"  rho        ann: half-width of the square.  log: outer radius |C| = rho\n"
"  survdepth  depth of the survival test that decides whether C lies in T_sigma at all.\n"
"             0 does not skip the test: it leaves the single check |Y(C)| <= R_Gamma,\n"
"             which is sound but so weak that the T_sigma column tends to read 100%%.\n"
"             Each extra level prunes more (at s0, rho=0.16: 100%% at 0 against 95.5%%\n"
"             at 8), so use 8 or more when that column is meant to mean something.\n"
"  out.bin    two bytes per pixel: (level, gamma).  level = least c certifying that pixel,\n"
"             255 = never certified, 254 = outside the plotted disc; gamma = C in T_sigma.\n"
"             Render it with:  render_funddom.py out.bin RESX RESY out.png\n"
"  nper       number of fundamental periods in log mode (default 1)\n"
"  balls.txt  trap-like balls, one per line as 'center_re center_im radius', in schottky's\n"
"             own normalization f(x)=zx, g(x)=z(x-1)+1 -- i.e. exactly what\n"
"             'certify_arc dumptlb' and 'certify_arc dumptlbmany' print.  Without this the\n"
"             built-in 10-ball sets for s0 and s1 are used, which are far too small to\n"
"             cover anything; generate a proper set with\n"
"                 certify_arc dumptlbmany 45 22 1e-9 40 14 8 > tlb.txt\n"
"                 python3 prune_tlb.py tlb.txt tlb_pruned.txt\n"
"\n"
"There is also a listing mode, which takes no other arguments:\n"
"  %s landmarks <Nmax>\n"
"             list the LANDMARK POINTS -- the renormalization points, i.e. the sigma at\n"
"             which the limit-trap mechanism applies -- of complexity a+b <= Nmax, one\n"
"             per line as\n"
"                 sigma_re sigma_im |sigma| arg_deg a b deg Delta_re Delta_im Pp_re Pp_im\n"
"             so each line carries everything a coverage run needs.  These are the roots\n"
"             in 1/2 < |sigma| < 1, Im sigma > 0, of Q(z) = A(z)(1-z^b) + z^a B(z) over\n"
"             all A, B with coefficients in {0,+-1} and A_0 = -1.  Unlike coin:, this\n"
"             reaches the INFINITE coincidences, so s0, s1 and hex all appear (marked).\n"
"             Counts: 1, 18, 99, 533, 2421, 10958, 46201 for Nmax = 3..9.  Cost grows\n"
"             like (N-1)*3^(N-1) polynomials: Nmax 8 is under a second, 9 a few seconds,\n"
"             10 about 16s and enough to reach the 100000-point cap.\n"
"\n"
"example: how much of one fundamental domain at the twindragon core is covered?\n"
"  certify_arc dumptlbmany 45 22 1e-9 40 14 8 > tlb45.txt\n"
"  python3 prune_tlb.py tlb45.txt tlb45_pruned.txt\n"
"  %s s0 log 24 700 155 0.08 0 s0.bin 1 tlb45_pruned.txt\n",
    prog, prog, prog, prog, prog, prog, prog);
}

/*==========================================================================*
 * Landmark points: enumerating the renormalization points of low complexity
 *
 * sigma is a renormalization point when pi(u s^inf, sigma) = pi(v t^inf, sigma)
 * for words u, v of length a and s, t of length b -- that is exactly the
 * hypothesis of CKW Lemma 9.2.5, so these are the parameters at which the
 * limit-trap mechanism applies.  Writing d_j = (eps^u_j - eps^v_j)/2 in
 * {0,+-1} for the difference of the two coding sequences, the sequence d is
 * some A of length a followed by some B of length b repeated forever, and
 *
 *     sum_j d_j z^j  =  A(z) + z^a B(z)/(1 - z^b),
 *
 * so the coincidence condition, cleared of its denominator, is
 *
 *     Q(z) := A(z)(1 - z^b) + z^a B(z) = 0,      A_0 = -1,  A,B in {0,+-1}.
 *
 * Everything funddom needs then comes straight out of Q:
 *
 *     Delta = p_u - p_v = 2 A(sigma),      P(z) = 2 Q(z)/(1 - z^b),
 *     and since Q(sigma) = 0,   P'(sigma) = 2 Q'(sigma)/(1 - sigma^b).
 *
 * VERIFIED against all three hard-coded cores: this recovers s0 at (a,b)=(2,1)
 * with A=(-1,1), B=(-1); s1 at (1,4) with A=(-1), B=(1,-1,-1,1); and the CKW
 * hexahole at (8,1) with A=(-1,1,-1,-1,-1,1,1,1), B=(-1) -- in each case by a
 * unique (A,B) -- and the Delta and P'(sigma) above agree with the literals in
 * the s0/s1/hex branches below.  If you change any of this, re-run
 * `funddom landmarks 8` and check the three ==> marks still appear.
 *
 * Note the enumeration is by a+b, not by the degree of Q (which is a+b-1 once
 * the leading coefficient survives).  Complexity a+b is the honest knob: it is
 * the total length of the defining word data.
 *==========================================================================*/

#endif /* FUNDDOM_MAIN: usage() is only needed by the CLI */

#define LM_MAXN 12       /* a+b ceiling: cost grows like (N-1)*3^(N-1) polynomials */

/* ia and ib are the base-3 packings of A and B, kept so that Delta and P'(sigma)
   can be rebuilt on output without storing a coefficient vector per point. */
typedef struct { cx sigma; int a, b, deg, ia, ib; } landmark;

/* unpack A (with A_0 = -1) and B, and build Q's coefficients */
static int lm_build(int a, int b, int ia, int ib, int *A, int *B, double *c)
{
  A[0] = -1;
  { int t = ia; for (int j = 1; j < a; ++j) { A[j] = t % 3 - 1; t /= 3; } }
  { int t = ib; for (int k = 0; k < b; ++k) { B[k] = t % 3 - 1; t /= 3; } }
  for (int j = 0; j <= a + b; ++j) c[j] = 0.0;
  for (int j = 0; j < a; ++j) { c[j] += A[j]; c[j + b] -= A[j]; }
  for (int k = 0; k < b; ++k) c[a + k] += B[k];
  int D = a + b;
  while (D > 0 && fabs(c[D]) < 1e-12) --D;
  return D;
}

static int lm_cmp_key(const void *p, const void *q)
{ double x = *(const double*)p, y = *(const double*)q; return x < y ? -1 : (x > y ? 1 : 0); }

/* Enumerate landmark points with a+b <= N.  Returns the count, or -1 on a bad N.
   Points are appended to out[] up to maxout.  Roots on |z| = 1 (which the factor
   (1 - z^b) contributes whenever B is identically zero) and outside the annulus
   1/2 < |z| < 1 are dropped, as are conjugates: Im sigma > 0 only, since a
   conjugate pair gives mirror-image pictures. */
static int landmarks(int N, landmark *out, int maxout, int verbose)
{
  if (N < 2 || N > LM_MAXN) return -1;
  int nout = 0;
  for (int tot = 2; tot <= N; ++tot) {
    for (int a = 1; a < tot; ++a) {
      int b = tot - a;
      int nA = 1, nB = 1;                       /* 3^(a-1) choices of A, 3^b of B */
      for (int i = 0; i < a - 1; ++i) nA *= 3;
      for (int i = 0; i < b; ++i) nB *= 3;
      for (int ia = 0; ia < nA; ++ia) {
        for (int ib = 0; ib < nB; ++ib) {
          int A[LM_MAXN], B[LM_MAXN];
          double c[2 * LM_MAXN + 1];
          int D = lm_build(a, b, ia, ib, A, B, c);
          if (D < 1) continue;
          cx root[2 * LM_MAXN + 1];
          poly_roots_(c, D, root);
          for (int i = 0; i < D; ++i) {
            double ab = cabs(root[i]);
            /* the margins matter: without them the roots of unity coming from
               (1 - z^b) sneak in at |z| = 1 - 1e-16 and are reported as landmarks */
            if (ab <= 0.5 + 1e-9 || ab >= 1.0 - 1e-9) continue;
            if (cimag(root[i]) <= 1e-9) continue;
            if (cabs(polyval_(c, D, root[i])) > 1e-9) continue;
            int dup = 0;
            for (int m = 0; m < nout; ++m)
              if (cabs(out[m].sigma - root[i]) < 1e-8) { dup = 1; break; }
            if (dup) continue;                  /* keep the first, i.e. least a+b */
            if (nout >= maxout) {
              if (verbose) fprintf(stderr, "landmarks: hit the %d-point cap\n", maxout);
              return nout;
            }
            out[nout].sigma = root[i]; out[nout].a = a; out[nout].b = b;
            out[nout].deg = D; out[nout].ia = ia; out[nout].ib = ib; ++nout;
          }
        }
      }
    }
    if (verbose) fprintf(stderr, "  a+b <= %d : %d points\n", tot, nout);
  }
  (void)lm_cmp_key;
  return nout;
}


/* Parse a sign string over {-,0,+} into coefficients {-1,0,+1}.  Returns the length,
   or -1 on a bad character or an over-long string. */
static int lm_signs(const char *t, int *out)
{
  int n = 0;
  for (; t[n]; ++n) {
    if (n >= LM_MAXN) return -1;
    if      (t[n] == '-') out[n] = -1;
    else if (t[n] == '0') out[n] =  0;
    else if (t[n] == '+') out[n] = +1;
    else return -1;
  }
  return n;
}

/* Derive the renormalization data of the landmark point given by sign strings A and B.
   Picks root `pick` among the admissible ones.  Returns the number of admissible roots,
   0 if none, or -1 on malformed input. */
static int lm_data(const char *As, const char *Bs, int pick, cx *sigma_out,
                   int *a_out, int *b_out, cx *Delta_out, cx *Pp_out, int verbose)
{
  int A[LM_MAXN], B[LM_MAXN];
  int a = lm_signs(As, A), b = lm_signs(Bs, B);
  if (a < 1 || b < 1) return -1;
  if (A[0] != -1) return -1;                  /* normalization: A_0 = -1 */
  double c[2 * LM_MAXN + 1];
  for (int j = 0; j <= a + b; ++j) c[j] = 0.0;
  for (int j = 0; j < a; ++j) { c[j] += A[j]; c[j + b] -= A[j]; }
  for (int k = 0; k < b; ++k) c[a + k] += B[k];
  int D = a + b;
  while (D > 0 && fabs(c[D]) < 1e-12) --D;
  if (D < 1) return 0;
  cx root[2 * LM_MAXN + 1];
  poly_roots_(c, D, root);
  cx adm[2 * LM_MAXN + 1]; int n = 0;
  for (int i = 0; i < D; ++i) {
    double ab = cabs(root[i]);
    if (ab <= 0.5 + 1e-9 || ab >= 1.0 - 1e-9) continue;
    if (cimag(root[i]) <= 1e-9) continue;
    if (cabs(polyval_(c, D, root[i])) > 1e-9) continue;
    adm[n++] = root[i];
  }
  if (verbose) {
    fprintf(stderr, "lm: A=%s B=%s  a=%d b=%d  Q(z) =", As, Bs, a, b);
    for (int j = 0; j <= D; ++j) if (c[j] != 0.0) fprintf(stderr, " %+g z^%d", c[j], j);
    fprintf(stderr, "\n    admissible roots (1/2<|z|<1, Im z>0): %d\n", n);
    for (int i = 0; i < n; ++i)
      fprintf(stderr, "     [%d] %.15f%+.15fi  |z|=%.9f arg=%.5f deg\n", i,
              creal(adm[i]), cimag(adm[i]), cabs(adm[i]), carg(adm[i]) * 180.0 / M_PI);
  }
  if (n == 0) return 0;
  if (pick < 0 || pick >= n) return -2;
  cx z = adm[pick];
  cx Delta = 0.0;
  for (int j = a - 1; j >= 0; --j) Delta = Delta * z + A[j];
  Delta *= 2.0;                                /* Delta = p_u - p_v = 2 A(sigma) */
  *sigma_out = z; *a_out = a; *b_out = b;
  *Delta_out = Delta;
  *Pp_out = 2.0 * polyder_(c, D, z) / (1.0 - cpow(z, b));
  return n;
}

/*===========================================================================*
 * The public C API declared in funddom_core.h.  funddom.c is compiled twice:
 * with -DFUNDDOM_MAIN it is the command-line tool, without it is the library
 * object the GUI links.  So there is one copy of the mathematics.
 *===========================================================================*/

/* The three built-in cores are themselves landmark points, and `funddom
   landmarks` recovers each with exactly the Delta and P'(sigma) hard-coded in the
   s0/s1/hex branches -- verified, including a byte-identical coverage raster for
   hex.  So rather than duplicate those literals here, delegate. */
int fd_core_builtin(const char *name, fd_core *out)
{
  if (!name || !out) return -1;
  if (!strcmp(name, "s0"))  return fd_core_from_lm("-+", "-", 0, out) > 0 ? 0 : -1;
  if (!strcmp(name, "s1"))  return fd_core_from_lm("-", "+--+", 0, out) > 0 ? 0 : -1;
  if (!strcmp(name, "hex")) return fd_core_from_lm("-+---+++", "-", 1, out) > 0 ? 0 : -1;
  return -1;
}

int fd_core_from_lm(const char *A, const char *B, int pick, fd_core *out)
{
  cx sigma, Delta, Pp; int a, b;
  if (!A || !B || !out) return -1;
  int n = lm_data(A, B, pick, &sigma, &a, &b, &Delta, &Pp, 0);
  if (n <= 0) return n;
  out->sigma_re = creal(sigma); out->sigma_im = cimag(sigma);
  out->a = a; out->b = b;
  out->Delta_re = creal(Delta); out->Delta_im = cimag(Delta);
  out->Pp_re = creal(Pp);       out->Pp_im = cimag(Pp);
  return n;
}

int fd_landmarks(int N, fd_landmark *out, int maxout)
{
  if (!out || maxout < 1) return -1;
  landmark *L = malloc(sizeof(landmark) * (size_t)maxout);
  if (!L) return -1;
  int n = landmarks(N, L, maxout, 0);
  for (int i = 0; i < n && i < maxout; ++i) {
    int A[LM_MAXN], B[LM_MAXN]; double c[2 * LM_MAXN + 1];
    int D = lm_build(L[i].a, L[i].b, L[i].ia, L[i].ib, A, B, c);
    char As[LM_MAXN + 1], Bs[LM_MAXN + 1];
    for (int j = 0; j < L[i].a; ++j) As[j] = (A[j] < 0 ? '-' : (A[j] > 0 ? '+' : '0'));
    As[L[i].a] = 0;
    for (int k = 0; k < L[i].b; ++k) Bs[k] = (B[k] < 0 ? '-' : (B[k] > 0 ? '+' : '0'));
    Bs[L[i].b] = 0;
    int pick = 0;                       /* which admissible root this one is */
    { cx rr[2 * LM_MAXN + 1]; poly_roots_(c, D, rr);
      for (int m = 0; m < D; ++m) {
        double ab = cabs(rr[m]);
        if (ab <= 0.5 + 1e-9 || ab >= 1.0 - 1e-9) continue;
        if (cimag(rr[m]) <= 1e-9) continue;
        if (cabs(polyval_(c, D, rr[m])) > 1e-9) continue;
        if (cabs(rr[m] - L[i].sigma) < 1e-8) break;
        ++pick;
      } }
    out[i].sigma_re = creal(L[i].sigma); out[i].sigma_im = cimag(L[i].sigma);
    out[i].a = L[i].a; out[i].b = L[i].b; out[i].deg = L[i].deg;
    snprintf(out[i].spec, sizeof out[i].spec, "lm:%s:%s:%d", As, Bs, pick);
  }
  free(L);
  return n;
}

int fd_solver_init(const fd_core *c, const double *balls3, int nballs)
{
  if (!c || !balls3 || nballs < 1) return -1;
  cx sigma = c->sigma_re + c->sigma_im * I;
  cx Delta = c->Delta_re + c->Delta_im * I;
  cx Pp    = c->Pp_re    + c->Pp_im    * I;
  fd_solver_free();
  SIG = sigma; ISIG = 1.0 / sigma;
  RG  = 2.0 / (1.0 - cabs(sigma));
  /* schottky's normalization f(x)=zx, g(x)=z(x-1)+1 is conjugate to ours by
     x -> (2x-1)/(1-z), so trap-like VECTORS scale by 2/(1-sigma).  Applied here
     so no caller can forget it. */
  cx conv = 2.0 / (1.0 - sigma);
  double tlbscale = cabs(conv);
  NB = nballs;
  BC = malloc(NB * sizeof(cx)); BR = malloc(NB * sizeof(double));
  BCR = malloc(NB * sizeof(double)); BCI = malloc(NB * sizeof(double));
  BR2 = malloc(NB * sizeof(double));
  if (!BC || !BR || !BCR || !BCI || !BR2) { fd_solver_free(); return -1; }
  RG2 = RG * RG; IR = creal(ISIG); II = cimag(ISIG);
  double rb = 0.0;
  for (int i = 0; i < NB; ++i) {
    BC[i] = conv * (balls3[3*i] + balls3[3*i+1] * I);
    BR[i] = tlbscale * balls3[3*i+2];
    BCR[i] = creal(BC[i]); BCI[i] = cimag(BC[i]); BR2[i] = BR[i] * BR[i];
    if (cabs(BC[i]) + BR[i] > rb) rb = cabs(BC[i]) + BR[i];
  }
  RBALLS2 = rb * rb;
  build_grid();
  Y0_S     = -Delta / cpow(sigma, c->a);
  SCALEC_S = -Pp    / cpow(sigma, c->a);       /* Y = Y0 + scaleC * C */
  return 0;
}

void fd_solver_free(void)
{
  free(BC);  BC  = NULL;  free(BR);  BR  = NULL;
  free(BCR); BCR = NULL;  free(BCI); BCI = NULL;
  free(BR2); BR2 = NULL;
  free(GSTART); GSTART = NULL;  free(GLIST); GLIST = NULL;
  NB = 0;
}

void fd_solver_info(fd_info *out)
{
  if (!out) return;
  out->R_Gamma = RG;
  out->Y0_re = creal(Y0_S);      out->Y0_im = cimag(Y0_S);
  out->dYdC_re = creal(SCALEC_S); out->dYdC_im = cimag(SCALEC_S);
  /* NOT R_Gamma/|dY/dC|: Y is centered at Y0, not 0, and the two differ by a
     factor of about 1.7 at s0 and 1.8 at the hexahole. */
  out->rho_max = (cabs(SCALEC_S) > 0.0)
               ? (RG - cabs(Y0_S)) / cabs(SCALEC_S) : 0.0;
  out->period_height = 0.0;      /* filled below if sigma is sane */
  if (cabs(SIG) > 0.0 && cabs(SIG) < 1.0) out->period_height = -log(cabs(SIG));
  out->nballs = NB;
}

int fd_level(double C_re, double C_im, int cmax)
{
  if (NB < 1) return -1;
  cx Y = Y0_S + SCALEC_S * (C_re + C_im * I);
  BEST = cmax + 1;
  dfs(creal(Y), cimag(Y), 0, cmax);
  return (BEST <= cmax) ? BEST : -1;
}

int fd_survives(double C_re, double C_im, int depth)
{
  if (NB < 1) return 0;
  cx Y = Y0_S + SCALEC_S * (C_re + C_im * I);
  return surv(creal(Y), cimag(Y), 0, depth) ? 1 : 0;
}


#ifdef FUNDDOM_MAIN

int main(int argc, char **argv)
{
  /* `funddom landmarks <Nmax>` lists the renormalization points of complexity
     a+b <= Nmax, with the renormalization data each one needs. */
  if (argc >= 2 && !strcmp(argv[1], "landmarks")) {
    if (argc < 3) {
      fprintf(stderr, "usage: %s landmarks <Nmax>   (complexity a+b, 2..%d)\n"
                      "  prints: sigma_re sigma_im |sigma| arg_deg a b deg "
                      "Delta_re Delta_im Pp_re Pp_im\n", argv[0], LM_MAXN);
      return 1;
    }
    int N = atoi(argv[2]);
    static landmark L[100000];
    int n = landmarks(N, L, 100000, 1);
    if (n < 0) { fprintf(stderr, "landmarks: Nmax must be in 2..%d\n", LM_MAXN); return 1; }
    printf("# landmark points (renormalization points) with a+b <= %d\n", N);
    printf("# sigma_re sigma_im |sigma| arg_deg a b deg Delta_re Delta_im Pp_re Pp_im spec\n");
    printf("# 'spec' feeds straight back in: funddom <spec> <mode> ... to run a coverage\n");
    const cx s0 = 0.5 + 0.5 * I;
    const cx s1 = 0.25 + (sqrt(7.0) / 4.0) * I;
    const cx hx = 0.371858680074136 + 0.519411153747943 * I;
    for (int i = 0; i < n; ++i) {
      cx z = L[i].sigma;
      int A[LM_MAXN], B[LM_MAXN]; double c[2 * LM_MAXN + 1];
      int D = lm_build(L[i].a, L[i].b, L[i].ia, L[i].ib, A, B, c);
      cx Delta = 0.0;
      for (int j = L[i].a - 1; j >= 0; --j) Delta = Delta * z + A[j];
      Delta *= 2.0;                                  /* Delta = p_u - p_v = 2 A(sigma) */
      cx Pp = 2.0 * polyder_(c, D, z) / (1.0 - cpow(z, L[i].b));
      const char *tag = "";
      if (cabs(z - s0) < 1e-9) tag = "   <== s0 (twindragon)";
      else if (cabs(z - s1) < 1e-9) tag = "   <== s1 (tame twindragon)";
      else if (cabs(z - hx) < 1e-9) tag = "   <== hex (CKW hexahole, on dM)";
      /* the sign strings, and which admissible root this is, so the line can be fed
         straight back to the lm: core selector */
      char As[LM_MAXN + 1], Bs[LM_MAXN + 1];
      for (int j = 0; j < L[i].a; ++j) As[j] = (A[j] < 0 ? '-' : (A[j] > 0 ? '+' : '0'));
      As[L[i].a] = 0;
      for (int k = 0; k < L[i].b; ++k) Bs[k] = (B[k] < 0 ? '-' : (B[k] > 0 ? '+' : '0'));
      Bs[L[i].b] = 0;
      int pick = 0;
      { cx rr[2 * LM_MAXN + 1]; poly_roots_(c, D, rr);
        for (int m = 0; m < D; ++m) {
          double ab = cabs(rr[m]);
          if (ab <= 0.5 + 1e-9 || ab >= 1.0 - 1e-9) continue;
          if (cimag(rr[m]) <= 1e-9) continue;
          if (cabs(polyval_(c, D, rr[m])) > 1e-9) continue;
          if (cabs(rr[m] - z) < 1e-8) break;
          ++pick;
        } }
      printf("%.15f %.15f  %.9f %10.5f  %d %d %d  %+.9f %+.9f  %+.9f %+.9f  lm:%s:%s:%d%s\n",
             creal(z), cimag(z), cabs(z), carg(z) * 180.0 / M_PI,
             L[i].a, L[i].b, L[i].deg,
             creal(Delta), cimag(Delta), creal(Pp), cimag(Pp), As, Bs, pick, tag);
    }
    fprintf(stderr, "%d landmark points with a+b <= %d\n", n, N);
    return 0;
  }

  /* `funddom coin:<u>:<v>` on its own is a query: derive and print the renormalization
     data for that coincidence and stop, so one can inspect a coincidence before
     committing to a raster. */
  if (argc < 9 && argc >= 2 && !strncmp(argv[1], "lm:", 3)) {
    char buf[512]; snprintf(buf, sizeof buf, "%s", argv[1] + 3);
    char *As = buf, *Bs = strchr(buf, ':'); int pick = 0;
    if (!Bs) { fprintf(stderr, "lm: expected lm:<A>:<B>[:index], signs over {-,0,+}\n"); return 1; }
    *Bs++ = 0;
    { char *k = strchr(Bs, ':'); if (k) { *k = 0; pick = atoi(k + 1); } }
    cx z, Delta, Pp; int a, b;
    int n = lm_data(As, Bs, pick, &z, &a, &b, &Delta, &Pp, 1);
    if (n == -1) { fprintf(stderr, "lm: A and B must be nonempty strings over {-,0,+} with A "
                                   "starting '-'\n"); return 1; }
    if (n == 0)  { fprintf(stderr, "lm: no admissible root\n"); return 1; }
    if (n == -2) { fprintf(stderr, "lm: root index out of range\n"); return 1; }
    printf("sigma = %.15f%+.15fi   |sigma| = %.12f   arg = %.6f deg\n",
           creal(z), cimag(z), cabs(z), carg(z) * 180.0 / M_PI);
    printf("  a = %d, b = %d;  Delta = %.9f%+.9fi;  P'(sigma) = %.9f%+.9fi\n",
           a, b, creal(Delta), cimag(Delta), creal(Pp), cimag(Pp));
    return 0;
  }

  if (argc < 9 && argc >= 2 && !strncmp(argv[1], "coin:", 5)) {
    char buf[512]; snprintf(buf, sizeof buf, "%s", argv[1] + 5);
    char *us = buf, *vs = strchr(buf, ':');
    if (vs) { *vs++ = 0; char *k = strchr(vs, ':'); if (k) *k = 0; }
    if (us && vs) {
      int d[MAXM], D; cx rts[MAXM];
      int n = coincidence_data(us, vs, d, &D, rts, MAXM);
      if (n == -2) {
        fprintf(stderr, "coin: words longer than %d are not supported (got %d)\n",
                MAXM, (int)strlen(us));
        return 1;
      }
      if (n < 0) {
        fprintf(stderr, "coin: u and v must be equal-length f/g words starting with "
                        "different letters (got u=%s v=%s)\n", us, vs);
        return 1;
      }
      printf("coincidence u = %s , v = %s   (m = %d)\n", us, vs, (int)strlen(us));
      printf("  d = (");
      for (int j = 0; j < (int)strlen(us); ++j) printf("%s%+d", j ? "," : "", d[j]);
      printf(")\n  polynomial (degree %d):  ", D);
      { int first = 1;
        for (int j = 0; j <= D; ++j) {
          if (!d[j]) continue;
          if (d[j] > 0) printf(first ? "" : " + ");
          else          printf(first ? "-" : " - ");
          if (j == 0)      printf("1");
          else if (j == 1) printf("z");
          else             printf("z^%d", j);
          first = 0;
        }
        printf("  =  0\n"); }
      printf("  admissible roots (1/2 < |sigma| < 1, Im sigma > 0): %d\n", n);
      for (int i = 0; i < n; ++i) {
        cx sg = rts[i];
        cx Pp = 0; for (int j = 1; j <= D; ++j) Pp += 2.0 * j * d[j] * cpow(sg, j - 1);
        printf("   [%d] sigma = %.15f%+.15fi   |sigma| = %.12f   arg = %.6f deg\n",
               i, creal(sg), cimag(sg), cabs(sg), carg(sg) * 180.0 / M_PI);
        printf("       a = %d, b = 1, s = t = f;  Delta = 0;  P'(sigma) = %.9f%+.9fi\n",
               (int)strlen(us), creal(Pp), cimag(Pp));
      }
      if (n == 0)
        printf("  none: this coincidence has no root in the annulus where Lambda is "
               "connected with interior.\n");
      else
        printf("\nto compute the coverage of one fundamental domain, e.g.:\n"
               "  %s %s log 24 700 155 0.08 0 out.bin 1 balls.txt\n", argv[0], argv[1]);
      return 0;
    }
  }
  if (argc < 9) { usage(argv[0]); return 1; }
  const char *core = argv[1], *mode = argv[2];
  int cmax = atoi(argv[3]), RESX = atoi(argv[4]), RESY = atoi(argv[5]);
  double rho = atof(argv[6]);
  long ndropped = 0;
  int survd = atoi(argv[7]);
  const char *out = argv[8];

  /* [nper] and [balls.txt] are both optional trailing arguments.  Taking them
     strictly by position means that supplying a ball file but not nper parses the
     FILENAME as nper -- atoi gives 0 -- and then silently never loads the file, so
     the run reports the coverage of a single circle computed from the built-in
     balls.  So identify them by shape instead: an argument that does not look like
     a number is the ball file, wherever it appears. */
  const char *ballfile = NULL;
  int nper = 1, nper_given = 0;
  for (int i = 9; i < argc; ++i) {
    const char *A = argv[i];
    if (!strcmp(A, "-")) continue;                  /* explicit "no ball file"      */
    char *end = NULL;
    long val = strtol(A, &end, 10);
    if (end != A && *end == '\0') {                 /* a bare integer: that is nper */
      if (nper_given) { fprintf(stderr, "funddom: nper given twice ('%s')\n", A); return 1; }
      nper = (int)val; nper_given = 1;
    } else {
      if (ballfile) { fprintf(stderr, "funddom: two ball files ('%s' and '%s')\n",
                              ballfile, A); return 1; }
      ballfile = A;
    }
  }
  if (nper < 1) {
    fprintf(stderr, "funddom: nper must be at least 1 (got %d)\n", nper);
    return 1;
  }
  /* the raster stores the level in one byte with 254/255 reserved, and the
     per-period tallies below are fixed [64] arrays */
  if (cmax < 0 || cmax > 253) {
    fprintf(stderr, "funddom: cmax must be in 0..253 (got %d)\n", cmax);
    return 1;
  }
  if (nper > 64) {
    fprintf(stderr, "funddom: nper must be at most 64 (got %d)\n", nper);
    return 1;
  }
  if (RESX < 1 || RESY < 1) {
    fprintf(stderr, "funddom: RESX and RESY must be positive\n");
    return 1;
  }
  if (!(rho > 0.0)) {
    fprintf(stderr, "funddom: rho must be positive (got %g)\n", rho);
    return 1;
  }
  int a, b;
  cx Delta, Pp, sigma;
  double tlbscale;                 /* dumptlb balls live in the code normalization       */
                                   /* f(x)=zx, g(x)=z(x-1)+1; ours is that conjugated by */
                                   /* x -> (2x-1)/(1-z), so vectors scale by 2/(1-sigma). */
  static double raw[200000][3]; int nraw = 0;

  if (!strncmp(core, "coin:", 5)) {
    /* renormalization data derived from a finite coincidence u(0) = v(0) */
    char buf[512]; snprintf(buf, sizeof buf, "%s", core + 5);
    char *us = buf, *vs = strchr(buf, ':'); int pick = 0;
    if (!vs) { fprintf(stderr, "coin: expected coin:<u>:<v>[:index]\n"); return 1; }
    *vs++ = 0;
    { char *k = strchr(vs, ':'); if (k) { *k = 0; pick = atoi(k + 1); } }
    int d[MAXM], D; cx rts[MAXM];
    int n = coincidence_data(us, vs, d, &D, rts, MAXM);
    if (n == -2) { fprintf(stderr, "coin: words longer than %d are not supported (got %d)\n",
                           MAXM, (int)strlen(us)); return 1; }
    if (n < 0) { fprintf(stderr, "coin: illegal coincidence pair u=%s v=%s\n", us, vs); return 1; }
    if (n == 0) { fprintf(stderr, "coin: no root of the coincidence polynomial in "
                                  "1/2 < |sigma| < 1 with Im sigma > 0\n"); return 1; }
    if (pick < 0 || pick >= n) {
      fprintf(stderr, "coin: root index %d out of range (there are %d; run "
                      "'%s coin:%s:%s' to list them)\n", pick, n, argv[0], us, vs);
      return 1;
    }
    sigma = rts[pick];
    a = (int)strlen(us); b = 1;
    Delta = 0.0;                     /* u and v are the same affine map at sigma */
    Pp = 0; for (int j = 1; j <= D; ++j) Pp += 2.0 * j * d[j] * cpow(sigma, j - 1);
    nraw = 0;                        /* a ball file is required for coin: cores  */
    fprintf(stderr, "coin: u=%s v=%s root[%d]  sigma=%.12f%+.12fi  a=%d b=1  "
                    "P'=%.6f%+.6fi\n", us, vs, pick, creal(sigma), cimag(sigma), a,
            creal(Pp), cimag(Pp));
  } else if (!strncmp(core, "lm:", 3)) {
    /* a landmark point given by the sign strings of A and B; see `funddom landmarks` */
    char buf[512]; snprintf(buf, sizeof buf, "%s", core + 3);
    char *As = buf, *Bs = strchr(buf, ':'); int pick = 0;
    if (!Bs) { fprintf(stderr, "lm: expected lm:<A>:<B>[:index]\n"); return 1; }
    *Bs++ = 0;
    { char *k = strchr(Bs, ':'); if (k) { *k = 0; pick = atoi(k + 1); } }
    int n = lm_data(As, Bs, pick, &sigma, &a, &b, &Delta, &Pp, 1);
    if (n <= 0) { fprintf(stderr, "lm: could not derive a landmark from A=%s B=%s\n", As, Bs);
                  return 1; }
    nraw = 0;                        /* a ball file is required for lm: cores */
  } else if (!strcmp(core, "s0")) {
    sigma = 0.5 + 0.5 * I; a = 2; b = 1;
    Delta = 2.0 * (sigma - 1.0);                         /* p_fg - p_gf                  */
    Pp    = 2.0 * (4.0 * sigma - 2.0) / (sigma - 1.0);   /* P = 2(2z^2-2z+1)/(z-1)       */
    /* certify_arc dumptlb 45 22 1e-9 */
    double t[10][3] = {
      { 0.318538161,-0.191715745,0.060743600},{-0.191715745, 0.318538161,0.060743600},
      {-0.318538161, 0.191715745,0.060743600},{ 0.191715745,-0.318538161,0.060743600},
      { 0.333251953, 0.082766928,0.055621207},{-0.333251953, 0.082766928,0.055621207},
      {-0.333251953,-0.082766928,0.055621207},{ 0.333251953,-0.082766928,0.055621207},
      { 0.171623684, 0.234138035,0.031395016},{-0.098884128,-0.036369778,0.031395016}};
    memcpy(raw, t, sizeof t); nraw = 10;
  } else if (!strcmp(core, "hex")) {
    /* CONTROL: the CKW hexahole, proved to lie in the BOUNDARY of M (Thm 9.1.1).
     * Renormalization point for u=fgfffggg, v=gfgggfff (a=8), s=f, t=g (b=1);
     * P(z) = 2 p_omega(z)/(1-z) with p_omega = -1+2z-2z^2+2z^5-2z^8, so P(omega)=0.
     * Here the coincidence is UNIQUE, so CKW Thm 9.2.2(2) forces holes: the asymptotic
     * picture must have an invariant uncovered part.  Requires a ball file. */
    sigma = 0.371858680074136 + 0.519411153747943 * I; a = 8; b = 1;
    {
      cx z = sigma;
      Delta = 2.0 * (-1.0 + z - z*z - z*z*z - z*z*z*z
                     + z*z*z*z*z + z*z*z*z*z*z + z*z*z*z*z*z*z);
      cx dpo = 2.0 - 4.0*z + 10.0*z*z*z*z - 16.0*z*z*z*z*z*z*z;
      Pp = 2.0 * dpo / (1.0 - z);
    }
    nraw = 0;
  } else if (!strcmp(core, "s1")) {
    sigma = 0.25 + (sqrt(7.0) / 4.0) * I; a = 1; b = 4;
    Delta = -2.0;                                        /* p_f - p_g                    */
    {   /* P = 2(2z^4-z^3-z^2+z-1)/(1-z^4) */
      cx z = sigma, Np = 8.0*z*z*z - 3.0*z*z - 2.0*z + 1.0;
      Pp = 2.0 * Np / (1.0 - z*z*z*z);
    }
    /* certify_arc dumptlb 69.29518870 22 1e-9 */
    double t[10][3] = {
      { 0.508249253,-0.179628461,0.084159589},{-0.241750744, 0.481809365,0.084159589},
      {-0.508249253, 0.179628461,0.084159589},{ 0.241750744,-0.481809365,0.084159589},
      { 0.245889525, 0.291241421,0.057657069},{-0.379110475,-0.039477488,0.057657069},
      {-0.245889525,-0.291241421,0.057657069},{ 0.379110475, 0.039477488,0.057657069},
      { 0.129580122,-0.235650755,0.038736027},{ 0.067080127, 0.260427615,0.038736027}};
    memcpy(raw, t, sizeof t); nraw = 10;
  } else {
    /* Without this the old code fell through to s1, so a typo -- "S0", "hexahole",
       anything -- silently produced s1's answer with the bogus name echoed back as
       though it had been recognised, and exit 0. */
    fprintf(stderr, "funddom: unknown core '%s'\n\n", core);
    usage(argv[0]);
    return 1;
  }

  if (ballfile) {                     /* richer certified trap-like set, same normalization */
    FILE *f = fopen(ballfile, "r");
    if (!f) { perror(ballfile); return 1; }
    char line[512]; double x, y, r; nraw = 0;
    double fs_re = 0, fs_im = 0; int have_fs = 0;
    while (fgets(line, sizeof line, f)) {
      /* certify_arc's dumps carry "# s = <re> <im>": the parameter the balls were
         computed at.  Trap-like balls are only valid at their own s, and using a
         mismatched file is the one remaining way to get a confident, plausible,
         entirely meaningless coverage number -- so check it. */
      if (!have_fs && sscanf(line, "# s = %lf %lf", &fs_re, &fs_im) == 2) { have_fs = 1; continue; }
      if (line[0] != '#' && sscanf(line, "%lf %lf %lf", &x, &y, &r) == 3 && r > 0) {
        if (nraw < 200000) { raw[nraw][0] = x; raw[nraw][1] = y; raw[nraw][2] = r; ++nraw; }
        else ++ndropped;   /* silently ignoring these would understate the coverage */
      }
    }
    fclose(f);
    if (ndropped > 0)
      fprintf(stderr, "  warning: %s holds more than %d balls; %ld past the limit were\n"
                      "           ignored, so the coverage below is a LOWER bound.  Prune the\n"
                      "           file first (prune_tlb.py) rather than truncating it here.\n",
              ballfile, 200000, ndropped);
    if (have_fs) {
      double dd = cabs((fs_re + fs_im * I) - sigma);
      if (dd > 1e-9) {
        fprintf(stderr,
          "funddom: %s holds trap-like balls for s = %.12f%+.12fi, but this core has\n"
          "         sigma = %.12f%+.12fi (differing by %.3g).  Those balls are not valid\n"
          "         here and the coverage would be meaningless.  Regenerate them at sigma,\n"
          "         e.g. 'certify_arc dumptlbmanyz %.15f %.15f 20 1e-9 40 14 8'.\n",
          ballfile, fs_re, fs_im, creal(sigma), cimag(sigma), dd,
          creal(sigma), cimag(sigma));
        return 1;
      }
      fprintf(stderr, "  ball file parameter matches sigma to %.1e\n", dd);
    } else {
      fprintf(stderr, "  warning: %s has no '# s = re im' line, so it cannot be checked\n"
                      "  against sigma.  Regenerate it with a current certify_arc.\n", ballfile);
    }
    fprintf(stderr, "  read %d trap-like balls from %s\n", nraw, ballfile);
  }

  SIG = sigma; ISIG = 1.0 / sigma;
  RG  = 2.0 / (1.0 - cabs(sigma));
  cx conv = 2.0 / (1.0 - sigma);
  tlbscale = cabs(conv);
  NB = nraw;
  BC = malloc(NB * sizeof(cx)); BR = malloc(NB * sizeof(double));
  BCR = malloc(NB * sizeof(double)); BCI = malloc(NB * sizeof(double));
  BR2 = malloc(NB * sizeof(double));
  RG2 = RG * RG; IR = creal(ISIG); II = cimag(ISIG);
  double rb = 0.0;
  for (int i = 0; i < nraw; ++i) {
    BC[i] = conv * (raw[i][0] + raw[i][1] * I);
    BR[i] = tlbscale * raw[i][2];
    BCR[i] = creal(BC[i]); BCI[i] = cimag(BC[i]); BR2[i] = BR[i] * BR[i];
    if (cabs(BC[i]) + BR[i] > rb) rb = cabs(BC[i]) + BR[i];
  }
  if (NB == 0) {
    fprintf(stderr, "core %s has no built-in trap-like balls: pass a ball file as the last\n"
                    "argument.  Generate one with\n"
                    "  certify_arc dumptlbmany  <deg> <n_depth> <d> <ngaps> <ntrials> <nradial>   (on |s|=1/sqrt2)\n"
                    "  certify_arc dumptlbmanyz <re> <im> <n_depth> <d> <ngaps> <ntrials> <nradial>   (anywhere)\n"
                    "then prune it with  python3 prune_tlb.py in.txt out.txt\n", core);
    return 1;
  }
  RBALLS2 = rb * rb;
  build_grid();
  cx Y0 = -Delta / cpow(sigma, a);
  cx scaleC = -Pp / cpow(sigma, a);      /* Y = Y0 + scaleC * C */

  fprintf(stderr, "core %s: sigma=%.9f%+.9fi a=%d b=%d\n", core, creal(sigma), cimag(sigma), a, b);
  fprintf(stderr, "  Delta=%.6f%+.6fi  P'(sigma)=%.6f%+.6fi\n",
          creal(Delta), cimag(Delta), creal(Pp), cimag(Pp));
  fprintf(stderr, "  Y0 = %.6f%+.6fi  |Y0|=%.6f   dY/dC = %.6f%+.6fi (|.|=%.6f)\n",
          creal(Y0), cimag(Y0), cabs(Y0), creal(scaleC), cimag(scaleC), cabs(scaleC));
  fprintf(stderr, "  R_Gamma=%.6f  trap-like balls: %d, radii", RG, NB);
  for (int i = 0; i < NB && i < 6; ++i) fprintf(stderr, " %.4f", BR[i]);
  fprintf(stderr, "\n");

  unsigned char *img = malloc(2 * (size_t)RESX * RESY);
  if (!img) { fprintf(stderr, "funddom: out of memory for a %dx%d raster\n", RESX, RESY); return 1; }
  /* open the output before the raster loop: this computation can run for hours and
     it is cruel to discover an unwritable path afterwards */
  FILE *fp = fopen(out, "wb");
  if (!fp) { perror(out); free(img); return 1; }
  double L = -log(cabs(sigma));
  long ncov = 0, ngam = 0, ntot = 0;
  long hist[260]; memset(hist, 0, sizeof hist);
  long pcov[64], pgam[64], ptot[64];
  memset(pcov, 0, sizeof pcov); memset(pgam, 0, sizeof pgam); memset(ptot, 0, sizeof ptot);

  for (int iy = 0; iy < RESY; ++iy) {
    for (int ix = 0; ix < RESX; ++ix) {
      cx C;
      int inside = 1, per = 0;
      if (!strcmp(mode, "ann")) {
        double x = (2.0 * (ix + 0.5) / RESX - 1.0) * rho;
        double y = (2.0 * (iy + 0.5) / RESY - 1.0) * rho;
        C = x + y * I;
        if (cabs(C) > rho) inside = 0;
      } else {
        double frac = (iy + 0.5) / RESY;              /* 0 = outer edge |C|=rho */
        double xi = log(rho) - nper * b * L * frac;
        double th = 2.0 * M_PI * (ix + 0.5) / RESX;
        C = cexp(xi + th * I);
        per = (int)(frac * nper); if (per > 63) per = 63;
      }
      unsigned char lev = 254, gam = 0;
      if (inside) {
        cx Y = Y0 + scaleC * C;
        BEST = cmax + 1;
        dfs(creal(Y), cimag(Y), 0, cmax);
        lev = (BEST <= cmax) ? (unsigned char)BEST : 255;
        gam = (unsigned char)surv(creal(Y), cimag(Y), 0, survd);
        ++ntot; ++ptot[per];
        if (lev != 255) { ++ncov; ++hist[lev]; ++pcov[per]; }
        if (gam) { ++ngam; ++pgam[per]; }
      }
      img[2 * ((size_t)iy * RESX + ix)]     = lev;
      img[2 * ((size_t)iy * RESX + ix) + 1] = gam;
    }
    if (iy % 16 == 0) { fprintf(stderr, "\r  row %d/%d", iy, RESY); fflush(stderr); }
  }
  fprintf(stderr, "\r  done                    \n");
  if (fwrite(img, 1, 2 * (size_t)RESX * RESY, fp) != 2 * (size_t)RESX * RESY)
    { perror(out); fclose(fp); free(img); return 1; }
  if (fclose(fp)) { perror(out); free(img); return 1; }
  printf("%s %s cmax=%d rho=%g nper=%d : covered %.4f%%   in T_sigma %.4f%%   (n=%ld)\n",
         core, mode, cmax, rho, nper, 100.0 * ncov / ntot, 100.0 * ngam / ntot, ntot);
  /* Not gated on survd: with survdepth 0 the survival test is still the meaningful
     |Y| <= R_Gamma check, and survdepth 0 is what the README's own s0 example uses,
     so gating here made the warning silent in the commonest configuration.  Also
     warn when nothing was covered, which reads like a result but usually is not. */
  if (ngam == 0 || ncov == 0)
    fprintf(stderr,
      "warning: %s, so the figure above is vacuous rather than a\n"
      "         mathematical result.  Y(C) = Y0 + C dY/dC with |Y0| = %.4g and |dY/dC| = %.4g,\n"
      "         and the difference set containing Y has radius R_Gamma = %.4g, so C can only be\n"
      "         in range while |C| <= (R_Gamma - |Y0|)/|dY/dC| = %.3g.  You used rho = %g.\n"
      "         %s\n",
      (ngam == 0 ? "NO sampled C lies in T_sigma" : "NOT ONE sampled C was certified"),
      cabs(Y0), cabs(scaleC), RG, (RG - cabs(Y0)) / cabs(scaleC), rho,
      (rho > (RG - cabs(Y0)) / cabs(scaleC)
         ? "Retry with rho below that bound."
         : "rho is within that bound, so suspect too small a cmax or too few trap-like balls."));
  if (!strcmp(mode, "log"))
    for (int p = 0; p < nper; ++p)
      printf("  period %d (|C| ~ %.5f..%.5f): covered %.3f%%  in T_sigma %.3f%%\n", p,
             rho * exp(-b * L * (p + 1)), rho * exp(-b * L * p),
             100.0 * pcov[p] / ptot[p], 100.0 * pgam[p] / ptot[p]);
  printf("  first-certified level histogram (cumulative %%):");
  long run = 0;
  for (int c = 0; c <= cmax; ++c) { run += hist[c];
    if (hist[c]) printf(" c=%d:%.2f", c, 100.0 * run / ntot); }
  printf("\n");
  return 0;
}

#endif /* FUNDDOM_MAIN */
