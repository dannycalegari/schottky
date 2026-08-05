#!/usr/bin/env python3
"""spiral.py -- follow the renormalization spiral at a renormalization point.

If sigma is a renormalization point of period b -- so that the combinatorics of the limit set
repeat under z -> sigma^b z in parameter space -- then the parameters

    sigma + C sigma^{bn},   n = 0, 1, 2, ...

converge to sigma along a logarithmic spiral, and the local picture of M near sigma is
asymptotically invariant under C -> sigma^b C.  Two things can happen along such a spiral, and
they are the two sides of CKW Theorem 9.2.2:

  * C admits a LIMIT TRAP: the words u s^n x, v t^n y trap sigma + C sigma^{bn} for all large n,
    so the whole tail of the spiral lies in the interior of M.  (Use `funddom` to decide this
    for every C at once, over a fundamental domain of C^*/<sigma^b>.)
  * C lies outside T_sigma and the coincidence at sigma is unique: then sigma + C sigma^{bn} is
    DISCONNECTED for all large n, so the spiral is a chain of holes accumulating at sigma,
    which forces sigma into the boundary of M.

This script tests the second alternative directly, with the rigorous hole certificate of
holes.py: take an offset C, walk inward along C sigma^{bn}, and report at each step whether the
parameter is a certified hole.  A spiral that stays holes as far as one can compute is evidence
that sigma is on the boundary; a spiral that stops being holes says the offset is not in the
disconnected regime (or that the depth is no longer sufficient -- the certificate needs depth
growing like n, which is why the depth is scaled automatically).

usage
    spiral.py <re_sigma> <im_sigma> <b> <re_C> <im_C> [nmax] [depth0]

    spiral.py 0.371858680074136 0.519411153747943 1 0.05 0     the CKW hexahole, on dM
    spiral.py 0.25 0.6614378277661477 4 0.05 0                 the tame twindragon core
    spiral.py 0.5 0.5 1 0.05 0                                 the twindragon core

An easy way to pick C: run `holes.py scan <re> <im> <r>` to find a certified hole at distance r
from sigma, and use its offset (hole - sigma) as C -- then n = 0 is a hole by construction and
the question is whether the inward images stay holes.

WHAT TO EXPECT.  This is a probe, not a proof of boundary membership.  Each individual verdict
of "HOLE" is a rigorous certificate, but the chains reachable in practice are short: sweeping
48 angles at three radii around the CKW hexahole -- a parameter known to lie in dM -- the
longest unbroken chain found was 4.  Chains break for two quite different reasons, and the
script cannot tell them apart: either the offset lies in T_sigma after all (so the tail of that
spiral really is in the interior), or the backward-tree certificate has run out of depth.  The
diagnostic for the second case is to re-test the first failing point at much greater depth; if
it still yields nothing, the break is genuine.  For deciding which C admit limit traps -- the
complementary and much more informative question -- use `funddom`, which settles it for every C
at once over a fundamental domain of C^*/<sigma^b>.
"""
import cmath, math, sys, time

from holes import is_hole


def main(argv):
    if len(argv) < 6:
        print(__doc__)
        return 1
    sigma = complex(float(argv[1]), float(argv[2]))
    b = int(argv[3])
    C = complex(float(argv[4]), float(argv[5]))
    nmax = int(argv[6]) if len(argv) > 6 else 12
    depth0 = int(argv[7]) if len(argv) > 7 else 60

    lam = sigma ** b
    print("sigma = %.15f%+.15fi   |sigma| = %.12f   b = %d   |sigma^b| = %.12f"
          % (sigma.real, sigma.imag, abs(sigma), b, abs(lam)))
    print("offset C = %.9f%+.9fi   |C| = %.6g" % (C.real, C.imag, abs(C)))
    print("the depth grows with n: a hole at distance d needs depth ~ -log d, and")
    print("|C sigma^(bn)| shrinks geometrically, so depth = depth0 + n * b * log2(1/|sigma|) * 7/log2(10).")
    print()
    print("   n           |C sigma^bn|   depth      verdict")
    per = -math.log10(abs(lam))            # decades gained per step
    nh = 0
    for n in range(nmax + 1):
        off = C * lam ** n
        s = sigma + off
        if abs(s) >= 1.0:
            print("  %2d   %18.6e      --      |s| >= 1, outside the disc" % (n, abs(off)))
            continue
        depth = depth0 + int(7 * per * n)
        t = time.time()
        h = is_hole(s, depth)
        nh += h
        print("  %2d   %18.6e   %5d      %s   [%.1fs]"
              % (n, abs(off), depth,
                 "HOLE (certified not in M)" if h else "no certificate at this depth",
                 time.time() - t))
        sys.stdout.flush()
    print()
    print("%d of %d spiral points certified as holes." % (nh, nmax + 1))
    if nh == nmax + 1:
        print("An unbroken chain of holes accumulating at sigma: evidence that sigma is in dM.")
    elif nh == 0:
        print("No holes along this spiral: the offset is not in the disconnected regime.")
    else:
        print("The chain breaks; either the offset leaves the disconnected regime, or the")
        print("certificate needs more depth (raise depth0) beyond the break.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
