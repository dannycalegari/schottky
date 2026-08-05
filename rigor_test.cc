/* rigor_test.cc -- known-answer tests for rigor.h.
 *
 * rigor.h is the one file here whose entire value is soundness: when certify_box says
 * CERTIFIED, a claim in the paper rests on it.  Everything else in this repository can be
 * wrong in a way that shows up as a bad picture; this can be wrong in a way that shows up
 * as a theorem that is not true.  So it needs tests whose expected answers are known
 * independently of the code being tested.
 *
 * Three kinds of check, in increasing order of how much they would catch:
 *
 *  1. HAND-COMPUTED ANSWERS.  Cases whose result can be worked out on paper -- the exact
 *     range of a product of two intervals, |[-3,2]| = [0,3], that PI_I really does contain
 *     the true pi (it does, but only just: the nearest double to pi is 1.22e-16 away and
 *     one ulp there is 4.44e-16, so the one-ulp widening is load-bearing, not decorative).
 *
 *  2. THE CONTAINMENT PROPERTY, against an independent oracle.  This is the soundness
 *     invariant: for every operation, the computed interval must contain the true result
 *     for every point of the inputs.  Checking that needs arithmetic more exact than the
 *     arithmetic under test, so this program DUMPS a corpus of (inputs, computed output)
 *     in hex float -- exactly, no decimal rounding -- and rigor_test.py re-derives each
 *     result in exact rational arithmetic and checks containment.  A self-written
 *     verifier would share any misunderstanding; Python's Fraction does not.
 *
 *  3. THE FORMULA AGAINST A DIFFERENT IMPLEMENTATION.  rigor.h computes a word's disk
 *     centre from a closed form.  ifs.cc computes it by composing the maps.  The two were
 *     written separately, so agreement is evidence the closed form is right -- which no
 *     amount of interval-arithmetic testing would tell us.
 *
 * Build:  make rigortest
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <complex>

#include "rigor.h"
#include "ifs.h"

static int checks = 0, failures = 0;
static FILE* corpus = NULL;

static void ok(bool cond, const std::string& what) {
  ++checks;
  if (!cond) { ++failures; std::printf("  FAIL  %s\n", what.c_str()); }
  else       { std::printf("  ok    %s\n", what.c_str()); }
}

/* Record one operation for the exact-arithmetic verifier.  %a prints a double with no
   loss, which matters: a decimal round-trip could hide exactly the one-ulp errors this
   is looking for. */
static void dump1(const char* op, const rig::R& a, const rig::R& r) {
  if (corpus) std::fprintf(corpus, "%s %a %a -> %a %a\n", op, a.lo, a.hi, r.lo, r.hi);
}
static void dump2(const char* op, const rig::R& a, const rig::R& b, const rig::R& r) {
  if (corpus) std::fprintf(corpus, "%s %a %a %a %a -> %a %a\n",
                           op, a.lo, a.hi, b.lo, b.hi, r.lo, r.hi);
}
static void dumpc(const char* op, const rig::C& a, const rig::C& b, const rig::R& r) {
  if (corpus) std::fprintf(corpus, "%s %a %a %a %a %a %a %a %a -> %a %a\n", op,
                           a.re.lo, a.re.hi, a.im.lo, a.im.hi,
                           b.re.lo, b.re.hi, b.im.lo, b.im.hi, r.lo, r.hi);
}
static void dumpc2(const char* op, const rig::C& a, const rig::C& b, const rig::C& r) {
  if (corpus) std::fprintf(corpus, "%s %a %a %a %a %a %a %a %a -> %a %a %a %a\n", op,
                           a.re.lo, a.re.hi, a.im.lo, a.im.hi,
                           b.re.lo, b.re.hi, b.im.lo, b.im.hi,
                           r.re.lo, r.re.hi, r.im.lo, r.im.hi);
}

/*=========================== 1. hand-computed answers ======================*/

static void test_hand_computed() {
  using namespace rig;
  std::printf("\n-- hand-computed answers --\n");

  ok(radd(R(1.0), R(2.0)).lo <= 3.0 && radd(R(1.0), R(2.0)).hi >= 3.0,
     "[1,1] + [2,2] contains 3");
  ok(rsub(R(1.0), R(2.0)).lo <= -1.0 && rsub(R(1.0), R(2.0)).hi >= -1.0,
     "[1,1] - [2,2] contains -1");

  /* the exact range of [-2,3]*[-5,7] is [-15,21]: the four corner products are
     10, -14, -15, 21 */
  R p = rmul(R(-2.0, 3.0), R(-5.0, 7.0));
  ok(p.lo <= -15.0 && p.hi >= 21.0, "[-2,3] * [-5,7] contains [-15,21]");
  ok(p.lo > -16.0 && p.hi < 22.0,   "[-2,3] * [-5,7] is not grossly overwide");

  R q = rabsval(R(-3.0, 2.0));
  ok(q.lo <= 0.0 && q.hi >= 3.0, "|[-3,2]| contains [0,3]");
  R q2 = rabsval(R(2.0, 5.0));
  ok(q2.lo <= 2.0 && q2.hi >= 5.0, "|[2,5]| contains [2,5]");

  /* An interval must never be inverted, whichever order it was built in. */
  ok(R(3.0, 1.0).lo == 1.0 && R(3.0, 1.0).hi == 3.0, "R(3,1) normalises to [1,3]");

  /* PI_I must contain the true pi.  Compared against a 20-digit decimal literal, which
     is far more precise than double, so this is a real check rather than a tautology. */
  const long double true_pi = 3.14159265358979323846L;
  ok((long double)PI_I.lo <= true_pi && (long double)PI_I.hi >= true_pi,
     "PI_I contains the true pi");
  ok(PI_I.lo < PI_I.hi, "PI_I is not a degenerate point interval");

  /* 1/sqrt2, likewise. */
  const long double true_inv_sqrt2 = 0.70710678118654752440L;
  R is2 = inv_sqrt2();
  ok((long double)is2.lo <= true_inv_sqrt2 && (long double)is2.hi >= true_inv_sqrt2,
     "inv_sqrt2() contains the true 1/sqrt2");

  /* cos over [0,pi] attains -1, and the enumeration of multiples of pi is what finds it:
     an implementation that only looked at the endpoints would return [-1,1] here by luck
     but would miss the minimum on, say, [pi/2, 3pi/2]. */
  R c1 = cos_iv(R(0.0, 3.2));
  ok(c1.lo <= -1.0 && c1.hi >= 1.0, "cos_iv([0,3.2]) contains [-1,1]");
  R c2 = cos_iv(R(1.6, 4.7));            /* contains pi but neither 0 nor 2pi */
  ok(c2.lo <= -1.0, "cos_iv([1.6,4.7]) finds the minimum -1 at pi (interior)");
  R s1 = sin_iv(R(1.0, 2.0));            /* contains pi/2 */
  ok(s1.hi >= 1.0, "sin_iv([1,2]) finds the maximum 1 at pi/2 (interior)");
  R s2 = sin_iv(R(4.0, 5.0));            /* contains 3pi/2 */
  ok(s2.lo <= -1.0, "sin_iv([4,5]) finds the minimum -1 at 3pi/2 (interior)");

  /* A degenerate box in s must still produce an interval containing the float answer. */
  C s = s_of_deg_box(45.0, 45.0);
  cpxf sf = s_of_deg_f(45.0);
  ok(s.re.lo <= sf.real() && sf.real() <= s.re.hi &&
     s.im.lo <= sf.imag() && sf.imag() <= s.im.hi,
     "s_of_deg_box(45,45) contains s_of_deg_f(45)");

  /* A wider box must contain a narrower one. */
  C wide = s_of_deg_box(40.0, 50.0);
  ok(wide.re.lo <= s.re.lo && wide.re.hi >= s.re.hi &&
     wide.im.lo <= s.im.lo && wide.im.hi >= s.im.hi,
     "s_of_deg_box(40,50) contains s_of_deg_box(45,45)");
}

/*=========== 2. corpus for the exact-arithmetic containment check ==========*/

static void build_corpus() {
  using namespace rig;
  std::printf("\n-- corpus for the exact-arithmetic verifier --\n");
  /* Deliberately awkward: values straddling zero, values of very different magnitude,
     denormal-ish smallness, exact powers of two (where nextafter steps change size), and
     numbers with no finite binary expansion. */
  static const double ends[] = {
    0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 0.1, -0.1, 1.0/3.0, -1.0/3.0,
    1e-8, -1e-8, 1e8, -1e8, 1e-300, 3.14159265358979, -2.718281828459045,
    0.7071067811865476, 1e16, 4503599627370497.0 /* 2^52+1 */
  };
  const int NE = (int)(sizeof(ends)/sizeof(ends[0]));
  std::vector<R> ivs;
  for (int i=0;i<NE;++i) ivs.push_back(R(ends[i]));                  /* points */
  for (int i=0;i<NE;++i) for (int j=i+1;j<NE;++j) if (((i*7+j)%5)==0) /* thinned pairs */
    ivs.push_back(R(ends[i], ends[j]));

  long n = 0;
  for (size_t i=0;i<ivs.size();++i) {
    dump1("rsqrt",  ivs[i], rsqrt(ivs[i]));
    dump1("rabs",   ivs[i], rabsval(ivs[i]));
    dump1("rneg",   ivs[i], rneg(ivs[i]));
    n += 3;
    for (size_t j=0;j<ivs.size();j+=7) {
      dump2("radd", ivs[i], ivs[j], radd(ivs[i], ivs[j]));
      dump2("rsub", ivs[i], ivs[j], rsub(ivs[i], ivs[j]));
      dump2("rmul", ivs[i], ivs[j], rmul(ivs[i], ivs[j]));
      n += 3;
      /* rdiv REQUIRES 0 not in the divisor; feeding it one would be testing a
         precondition violation, not the code.  See the note in the report below. */
      if (ivs[j].lo > 0.0 || ivs[j].hi < 0.0) {
        dump2("rdiv", ivs[i], ivs[j], rdiv(ivs[i], ivs[j]));
        ++n;
      }
    }
  }
  /* complex ops */
  for (size_t i=0;i<ivs.size();i+=3) for (size_t j=0;j<ivs.size();j+=11) {
    C a(ivs[i], ivs[j]), b(ivs[j], ivs[i]);
    dumpc2("cadd", a, b, cadd(a,b));
    dumpc2("csub", a, b, csub(a,b));
    dumpc2("cmul", a, b, cmul(a,b));
    dumpc ("cabs", a, b, cabsval(a));      /* b ignored by the verifier for cabs */
    dumpc ("cdist", a, b, cdist(a,b));
    n += 5;
  }
  /* the trig intervals, over boxes chosen to straddle the extremes */
  for (double lo = -7.0; lo <= 7.0; lo += 0.37) {
    for (double w = 0.0; w <= 3.0; w += 0.73) {
      R th(lo, lo+w);
      dump1("cos", th, cos_iv(th));
      dump1("sin", th, sin_iv(th));
      n += 2;
    }
  }
  std::printf("  wrote %ld operations to the corpus\n", n);
}

/*========== 3. the closed form against ifs.cc's ball construction =========*/

/* rigor.h's word_center_* use a closed form.  ifs.cc builds the same disk by composing
   the maps: act_on_left for the first letter, act_on_right for the rest, from the seed
   ball of radius 1.01*minimal_enclosing_radius centred at 1/2.  If the closed form is
   wrong, no interval test would notice -- both the float and interval versions would be
   consistently wrong together.  This is the only check here that could catch that. */
static void test_formula_vs_ifs() {
  using namespace rig;
  std::printf("\n-- the closed form vs ifs.cc's own ball construction --\n");
  static const char* words[] = {
    "0", "1", "01", "10", "010", "101", "0100", "1011",
    "010010001001", "111100110111", "0100111000111", NULL };
  static const double degs[] = { 45.0, 60.0, 83.6, 30.0, 75.0 };

  double worst = 0.0;
  int compared = 0;
  for (int d=0; d<5; ++d) {
    cpxf s = s_of_deg_f(degs[d]);
    ifs F;
    F.initialize(cpx(s.real(), s.imag()), cpx(s.real(), s.imag()), 512, 0);
    double min_r;
    if (!F.minimal_enclosing_radius(min_r)) continue;
    Ball seed(0.5, (cpx(s.real(),s.imag())-1.0)/2.0,
                   (1.0-cpx(s.real(),s.imag()))/2.0, 1.01*min_r);
    for (int w=0; words[w]; ++w) {
      std::string u = words[w];
      Ball b = F.act_on_left(u[0]-'0', seed);
      for (size_t i=1;i<u.size();++i) b = F.act_on_right(u[i]-'0', b);
      cpxf closed = word_center_f(u, s);
      double err = std::abs(cpxf(b.center.real(), b.center.imag()) - closed);
      if (err > worst) worst = err;
      ++compared;
    }
  }
  std::printf("  compared %d (word, s) pairs\n", compared);
  ok(worst < 1e-12, "closed form agrees with the composed maps (worst "
                    + std::string(1, ' ') + std::to_string(worst) + ")");

  /* And the interval version must contain the float version, at a degenerate box. */
  int contained = 0, total = 0;
  for (int d=0; d<5; ++d) {
    C S = s_of_deg_box(degs[d], degs[d]);
    cpxf s = s_of_deg_f(degs[d]);
    for (int w=0; words[w]; ++w) {
      C iv = word_center_iv(words[w], S);
      cpxf f = word_center_f(words[w], s);
      ++total;
      if (iv.re.lo <= f.real() && f.real() <= iv.re.hi &&
          iv.im.lo <= f.imag() && f.imag() <= iv.im.hi) ++contained;
    }
  }
  ok(contained == total, "word_center_iv contains word_center_f at a degenerate box ("
                         + std::to_string(contained) + "/" + std::to_string(total) + ")");

  /* Nesting: a wider s-box must give a wider centre box. */
  int nested = 0; total = 0;
  for (int w=0; words[w]; ++w) {
    C narrow = word_center_iv(words[w], s_of_deg_box(59.99, 60.01));
    C wide   = word_center_iv(words[w], s_of_deg_box(59.0,  61.0));
    ++total;
    if (wide.re.lo <= narrow.re.lo && wide.re.hi >= narrow.re.hi &&
        wide.im.lo <= narrow.im.lo && wide.im.hi >= narrow.im.hi) ++nested;
  }
  ok(nested == total, "a wider s-box gives a containing centre box ("
                      + std::to_string(nested) + "/" + std::to_string(total) + ")");
}

/*==================== 4. the certificates, known answers ==================*/

static void test_certificates() {
  using namespace rig;
  std::printf("\n-- certificates, on cases whose answer is known independently --\n");

  /* A pair taken from the project's own certificate data, on the box it was recorded
     for.  Known answer: CERTIFIED. */
  const std::string U = "010010001001", V = "111100110111";
  double mm = -1;
  int rc = certify_box(U, V, 6, 83.625, 83.626, mm);
  ok(rc == 1 && mm > 0.0, "certify_box CERTIFIES a box from the certificate data");

  double mm2 = -1;
  int rc2 = certify_box_seg(U, V, 6, 83.625, 83.626, mm2);
  ok(rc2 == 1 && mm2 > 0.0, "certify_box_seg certifies the same box");

  /* A word cannot link with itself: the two clouds coincide, so no alternation exists.
     Whatever else happens, this must NOT come back CERTIFIED. */
  double mm3 = -1;
  ok(certify_box(U, U, 6, 83.625, 83.626, mm3) != 1,
     "certify_box refuses u linked with itself");

  /* u and v must start with different letters for the construction to mean anything;
     two words on the same side cannot interleave. */
  double mm4 = -1;
  ok(certify_box("0100", "0101", 6, 83.625, 83.626, mm4) != 1,
     "certify_box refuses two words on the same side");

  /* SOUNDNESS UNDER NARROWING.  A narrower box is strictly easier: if the wide box
     certifies, every sub-box must certify too.  A failure here means the interval
     bookkeeping is not monotone, which would undermine every covering run. */
  int subs_ok = 0, subs = 0;
  double lo = 83.625, hi = 83.626;
  for (int k=0;k<5;++k) {
    double a = lo + (hi-lo)*k/10.0, b = hi - (hi-lo)*k/10.0;
    double m = -1;
    ++subs;
    if (certify_box(U, V, 6, a, b, m) == 1) ++subs_ok;
  }
  ok(subs_ok == subs, "every sub-box of a certified box also certifies ("
                      + std::to_string(subs_ok) + "/" + std::to_string(subs) + ")");

  /* WIDENING must not manufacture a certificate.  Over a huge arc the alternation
     certainly does not persist, so the answer must be a refusal, not a claim. */
  double mm5 = -1;
  ok(certify_box(U, V, 6, 0.0, 90.0, mm5) != 1,
     "certify_box refuses a 90-degree box (cannot certify what is not true)");

  /* The thin-trap certificate.  Known answer at the twindragon core, from the paper's
     own construction. */
  double tm = -1;
  ok(certify_thin_trap("01","10","00","11", 44.9, 45.1, 10.0, tm) == 1 && tm > 0.0,
     "certify_thin_trap certifies the crossing at the 45-degree core");

  /* Two chords that are the SAME segment cannot cross transversely: r x s = 0, so the
     divisor is not certainly nonzero and it must refuse rather than divide by zero. */
  double tm2 = -1;
  ok(certify_thin_trap("01","10","01","10", 44.9, 45.1, 10.0, tm2) != 1,
     "certify_thin_trap refuses a degenerate (self) crossing");

  /* Demanding an impossible crossing angle must refuse. */
  double tm3 = -1;
  ok(certify_thin_trap("01","10","00","11", 44.9, 45.1, 89.9, tm3) != 1,
     "certify_thin_trap refuses when the required angle is too large");

  /* Margins are meaningful: a narrower box should not have a SMALLER margin than a
     wider one (more room, not less). */
  double mw = -1, mn = -1;
  certify_thin_trap("01","10","00","11", 44.0, 46.0, 10.0, mw);
  certify_thin_trap("01","10","00","11", 44.9, 45.1, 10.0, mn);
  ok(mn >= mw, "the thin-trap margin does not shrink as the box narrows");

  /* THE CASE THAT MATTERS MOST, and which an earlier version of this test missed.
     Everything above reaches a refusal through the FLOAT pre-check -- the pair is not a
     trap even numerically, so the interval stage never runs.  The dangerous path is the
     other one: a pair that IS a float trap on a box too wide for the intervals to
     separate.  There the only correct answer is -1 (interval-FAILED); returning 1 would
     convert "I could not verify this" into "CERTIFIED", which is the one failure mode
     that could put a false statement into a paper.

     Measured on the real pair: width 0.002 about 83.6255 certifies, width 0.005 and up
     does not.  So 0.01 is comfortably inside the interval-failed regime while still being
     a genuine float trap.  A mutant that accepts a non-positive margin passes every other
     check in this file and fails these. */
  double mf = -1;
  int rcf = certify_box(U, V, 6, 83.6205, 83.6305, mf);
  ok(rcf == -1, "certify_box returns interval-FAILED (not CERTIFIED) on a too-wide box "
                "that is still a float trap");
  double mfs = -1;
  ok(certify_box_seg(U, V, 6, 83.6205, 83.6305, mfs) == -1,
     "certify_box_seg likewise refuses the too-wide box");
  double mft = -1;
  ok(certify_thin_trap("01","10","00","11", 40.0, 50.0, 10.0, mft) == -1,
     "certify_thin_trap returns interval-FAILED on a 10-degree box");

  /* And the invariant that ties the margin to the verdict, over a sweep of widths: a
     CERTIFIED verdict must always come with a strictly positive margin.  This is what
     makes the margin worth printing, and it is violated by exactly the mutation above. */
  int swept = 0, consistent = 0;
  for (int k = 0; k < 24; ++k) {
    double w = 0.0002 * (k + 1) * (k + 1);
    double m = -1;
    int rc = certify_box(U, V, 6, 83.6255 - w/2, 83.6255 + w/2, m);
    ++swept;
    if (rc != 1 || m > 0.0) ++consistent;
  }
  ok(consistent == swept, "a CERTIFIED verdict always carries a positive margin, over "
                          + std::to_string(swept) + " box widths");

  int swept2 = 0, consistent2 = 0;
  for (int k = 0; k < 20; ++k) {
    double w = 0.05 * (k + 1);
    double m = -1;
    int rc = certify_thin_trap("01","10","00","11", 45.0 - w/2, 45.0 + w/2, 10.0, m);
    ++swept2;
    if (rc != 1 || m > 0.0) ++consistent2;
  }
  ok(consistent2 == swept2, "the same for the thin-trap certificate, over "
                            + std::to_string(swept2) + " widths");
}

int main(int argc, char** argv) {
  const char* corpus_path = (argc > 1) ? argv[1] : "test_out/rigor_corpus.txt";
  corpus = std::fopen(corpus_path, "w");
  if (!corpus) { std::fprintf(stderr, "cannot write %s\n", corpus_path); return 2; }

  std::printf("rigor_test -- known-answer tests for rigor.h\n");
  test_hand_computed();
  build_corpus();
  test_formula_vs_ifs();
  test_certificates();
  std::fclose(corpus);

  std::printf("\n%d checks, %d failure%s\n", checks, failures,
              failures==1?"":"s");
  std::printf("corpus written to %s -- now run rigor_test.py on it for the\n"
              "exact-arithmetic containment check.\n", corpus_path);
  return failures == 0 ? 0 : 1;
}
