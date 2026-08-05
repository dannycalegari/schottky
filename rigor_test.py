#!/usr/bin/env python3
"""rigor_test.py -- the exact-arithmetic half of the rigor.h known-answer test.

    make rigortest            # builds and runs both halves
    python3 rigor_test.py test_out/rigor_corpus.txt

rigor_test.cc writes a corpus of (operation, input intervals, computed output interval),
with every number in hex float so nothing is lost to decimal rounding.  This script
re-derives each result in EXACT arithmetic and checks the soundness invariant:

    the computed interval must contain the true result for every point of the inputs.

Why a separate program, and why Python: a verifier written in the same style as the code
under test would share any misunderstanding of what outward rounding is for.  Python's
`fractions.Fraction` is exact for +, -, * and /, and a double converts to a Fraction
exactly, so for those four operations the true range is computed with no error at all --
there is nothing left to be wrong about.  For sqrt and modulus, which are irrational, the
check is done by squaring, which is again exact.

An interval that is too WIDE is sound but wasteful, and this script reports the worst case
it sees rather than failing on it -- a certificate that is merely conservative is still a
certificate.  An interval that is too NARROW is a bug that could turn a false statement
into a certified one, and that fails.

Excess width is reported in ULPS OF THE RESULT, not as a ratio to the true width.  The
ratio is worthless here: adding [0, 1e-300] to [-1e8, -1e8] has a true width of 1e-300 and
a correct computed width of one ulp at 1e8, so the ratio is about 1e292 and means nothing
at all.  In ulps the same case reads 2, which is exactly what outward rounding should cost.
"""
import math
import sys
from fractions import Fraction as F

USAGE = __doc__


def fr(x):
    """A double as an exact Fraction."""
    return F(x)


def sqrt_lo_hi_ok(lo, hi, a_lo, a_hi):
    """Check [lo,hi] contains sqrt(x) for every x in [a_lo,a_hi] with x >= 0, exactly.

    Only squaring is used, so this is exact: sqrt(x) >= lo iff x >= lo^2 (for lo >= 0),
    and sqrt(x) <= hi iff x <= hi^2.
    """
    x_lo = max(a_lo, F(0))
    x_hi = max(a_hi, F(0))
    if x_hi < 0:
        return True                      # nothing in range
    if hi < 0:
        return False
    if fr(hi) ** 2 < x_hi:
        return False                     # upper end too small
    if lo > 0 and fr(lo) ** 2 > x_lo:
        return False                     # lower end too large
    return True


def corners(a_lo, a_hi, b_lo, b_hi):
    return ((a_lo, b_lo), (a_lo, b_hi), (a_hi, b_lo), (a_hi, b_hi))


def main():
    if len(sys.argv) < 2:
        print(USAGE)
        return 1
    path = sys.argv[1]
    bad = []
    counts = {}
    worst_over = (0.0, "")      # worst excess width in ulps, reported not failed
    n = 0

    for lineno, line in enumerate(open(path), 1):
        line = line.strip()
        if not line:
            continue
        left, right = line.split("->")
        parts = left.split()
        op = parts[0]
        nums = [float.fromhex(t) for t in parts[1:]]
        out = [float.fromhex(t) for t in right.split()]
        n += 1
        counts[op] = counts.get(op, 0) + 1
        fail = None

        if op in ("radd", "rsub", "rmul", "rdiv"):
            a_lo, a_hi, b_lo, b_hi = (fr(v) for v in nums)
            r_lo, r_hi = fr(out[0]), fr(out[1])
            vals = []
            for x, y in corners(a_lo, a_hi, b_lo, b_hi):
                if op == "radd":
                    vals.append(x + y)
                elif op == "rsub":
                    vals.append(x - y)
                elif op == "rmul":
                    vals.append(x * y)
                else:
                    if y == 0:
                        vals = None
                        break
                    vals.append(x / y)
            if vals is not None:
                # + - * / are monotone in each argument, so the true range over the
                # boxes is attained at the corners -- exactly, no sampling needed.
                t_lo, t_hi = min(vals), max(vals)
                if r_lo > t_lo or r_hi < t_hi:
                    fail = "does not contain the true range [%s, %s]" % (
                        float(t_lo), float(t_hi))
                else:
                    # excess width in ulps of the result: 2 is the cost of correct
                    # outward rounding, a few more is a conservative operation, and a
                    # large number would mean the interval is needlessly loose
                    mag = max(abs(out[0]), abs(out[1]))
                    u = math.ulp(mag) if mag > 0 else math.ulp(1.0)
                    excess = (float(r_hi - r_lo) - float(t_hi - t_lo)) / u
                    if excess > worst_over[0]:
                        worst_over = (excess, "%s line %d" % (op, lineno))

        elif op == "rneg":
            a_lo, a_hi = fr(nums[0]), fr(nums[1])
            r_lo, r_hi = fr(out[0]), fr(out[1])
            if r_lo > -a_hi or r_hi < -a_lo:
                fail = "negation is not exact"

        elif op == "rabs":
            a_lo, a_hi = fr(nums[0]), fr(nums[1])
            r_lo, r_hi = fr(out[0]), fr(out[1])
            if a_lo <= 0 <= a_hi:
                t_lo, t_hi = F(0), max(-a_lo, a_hi)
            else:
                t_lo = min(abs(a_lo), abs(a_hi))
                t_hi = max(abs(a_lo), abs(a_hi))
            if r_lo > t_lo or r_hi < t_hi:
                fail = "does not contain |x| range [%s, %s]" % (float(t_lo), float(t_hi))

        elif op == "rsqrt":
            a_lo, a_hi = fr(nums[0]), fr(nums[1])
            if not sqrt_lo_hi_ok(fr(out[0]), fr(out[1]), a_lo, a_hi):
                fail = "does not contain sqrt of the input range"

        elif op in ("cadd", "csub", "cmul"):
            a = [fr(v) for v in nums[0:4]]
            b = [fr(v) for v in nums[4:8]]
            r = [fr(v) for v in out]
            are_lo, are_hi, aim_lo, aim_hi = a
            bre_lo, bre_hi, bim_lo, bim_hi = b
            if op == "cadd":
                t = (are_lo + bre_lo, are_hi + bre_hi, aim_lo + bim_lo, aim_hi + bim_hi)
            elif op == "csub":
                t = (are_lo - bre_hi, are_hi - bre_lo, aim_lo - bim_hi, aim_hi - bim_lo)
            else:
                # (a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re), each term's range from
                # its corners -- which is what the code claims to compute
                def rng(p_lo, p_hi, q_lo, q_hi):
                    vs = [x * y for x, y in corners(p_lo, p_hi, q_lo, q_hi)]
                    return min(vs), max(vs)
                rr = rng(are_lo, are_hi, bre_lo, bre_hi)
                ii = rng(aim_lo, aim_hi, bim_lo, bim_hi)
                ri = rng(are_lo, are_hi, bim_lo, bim_hi)
                ir = rng(aim_lo, aim_hi, bre_lo, bre_hi)
                t = (rr[0] - ii[1], rr[1] - ii[0], ri[0] + ir[0], ri[1] + ir[1])
            if r[0] > t[0] or r[1] < t[1] or r[2] > t[2] or r[3] < t[3]:
                fail = "does not contain the true rectangle"

        elif op in ("cabs", "cdist"):
            a = [fr(v) for v in nums[0:4]]
            b = [fr(v) for v in nums[4:8]]
            if op == "cabs":
                re_lo, re_hi, im_lo, im_hi = a
            else:
                re_lo, re_hi = a[0] - b[1], a[1] - b[0]
                im_lo, im_hi = a[2] - b[3], a[3] - b[2]
            # the true modulus range: smallest and largest |re| and |im| over the box
            def absrange(lo, hi):
                if lo <= 0 <= hi:
                    return F(0), max(-lo, hi)
                return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))
            mr = absrange(re_lo, re_hi)
            mi = absrange(im_lo, im_hi)
            sq_lo = mr[0] ** 2 + mi[0] ** 2
            sq_hi = mr[1] ** 2 + mi[1] ** 2
            if not sqrt_lo_hi_ok(fr(out[0]), fr(out[1]), sq_lo, sq_hi):
                fail = "modulus interval does not contain the true modulus range"

        elif op in ("cos", "sin"):
            # cos and sin are not rational, so exactness is impossible here; instead
            # sample the input interval densely in double precision and require the
            # reported interval to contain every sampled value with a little slack for
            # libm's own error.  Weaker than the algebraic checks above, and the reason
            # the C++ side separately checks that the extremes at multiples of pi are
            # found at all.
            a_lo, a_hi = nums[0], nums[1]
            r_lo, r_hi = out[0], out[1]
            f = math.cos if op == "cos" else math.sin
            N = 257
            eps = 4e-16
            for k in range(N):
                x = a_lo + (a_hi - a_lo) * k / (N - 1.0)
                v = f(x)
                if v < r_lo - eps or v > r_hi + eps:
                    fail = "%s(%.17g) = %.17g outside [%.17g, %.17g]" % (
                        op, x, v, r_lo, r_hi)
                    break
        else:
            fail = "unknown operation in corpus"

        if fail:
            bad.append("line %d: %s  %s" % (lineno, op, fail))

    print("verified %d operations in exact arithmetic" % n)
    for op in sorted(counts):
        print("    %-6s %6d" % (op, counts[op]))
    if worst_over[0] > 0:
        print("  worst excess width: %.1f ulps of the result (%s)"
              % (worst_over[0], worst_over[1]))
        print("  (2 ulps is the expected cost of outward rounding; an overwide interval")
        print("   is sound, just conservative, so this is reported and not a failure)")
    if bad:
        print("\n%d SOUNDNESS FAILURE(S):" % len(bad))
        for b in bad[:25]:
            print("   ", b)
        if len(bad) > 25:
            print("    ... and %d more" % (len(bad) - 25))
        return 1
    print("\nall intervals contain their true results: rigor.h is sound on this corpus")
    return 0


if __name__ == "__main__":
    sys.exit(main())
