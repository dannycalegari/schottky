#!/usr/bin/env python3
"""holes.py -- rigorous detection of holes in the connectedness locus M.

For the similarity pair  f(z) = sz - 1,  g(z) = sz + 1  with attractor Lambda,

    s in M  <=>  f(Lambda) meets g(Lambda)  <=>  2/s  in  D_s := Lambda - Lambda,

and D_s is the attractor of the three contractions  z -> sz,  z -> sz + 2,  z -> sz - 2
(if p - q = 2/s then f(p) = g(q); differences of points of Lambda have digits 0, +-2).
Hence 2/s lies in D_s exactly when 2/s admits an infinite BACKWARD orbit

    a  ->  a/s,  (a-2)/s,  (a+2)/s

staying in the disc |z| <= R = 2/(1-|s|) that contains D_s.  Searching that ternary backward
tree decides membership one way:

    the tree DIES at a finite depth  =>  2/s not in D_s  =>  s NOT in M.

That is a finite rigorous certificate -- the subtree is exhausted, so no infinite backward
orbit exists.  The converse is not available: survival to the depth cap only means no
certificate was found there, so such an s is reported as "in M or undecided".  `is_hole` also
returns False when it exceeds `budget` nodes, which errs the same way; it never reports a hole
that is not one.

usage
    holes.py point  <re> <im> [depth] [budget]
    holes.py scan   <re> <im> <r> [ntheta] [depth]        circle of radius r about re+i*im
    holes.py arc    <deg_lo> <deg_hi> <steps> [|s|] [depth]      along a circle |s| = const
    holes.py map    <re> <im> <halfwidth> [RES] [depth] [out.png]     raster about a point
    holes.py region <x0> <y0> <x1> <y1> [RES] [depth] [out.png]       raster of a rectangle

Any argument ending in .png is taken as the output file wherever it appears, so the numeric
options in front of it really are optional.

In the rasters black = certified hole, white = in M or undecided, pale gray = outside the
annulus 1/2 < |s| < 1 where the question is settled anyway, and the circle |s| = 1/sqrt2 is
drawn in red where it crosses the window.  Depths default adaptively (a hole at distance r
from a renormalization point typically needs depth ~ -log r before it becomes visible).

examples
    holes.py point 0.0227 0.703                  a known hole
    holes.py point 0.5 0.5 200                   the twindragon parameter: survives
    holes.py scan  0.5 0.5 1e-3                  is there a hole near s_0?
    holes.py arc   0 90 400                      sweep the circle |s| = 1/sqrt2
    holes.py region 0 0 1 1 400 60 quadrant.png  M in the first quadrant
"""
import cmath, math, struct, sys, time, zlib


# ------------------------------------------------------------------ the hole test
def is_hole(s, depth_cap, budget=8000000):
    """True  => s is certified NOT in M (rigorous).
       False => 2/s survived to depth_cap (or to the node budget): in M, or undecided."""
    R = 2.0 / (1.0 - abs(s)) + 1e-9
    invs = 1.0 / s
    w = 2.0 / s
    if abs(w) > R:
        return True                        # 2/s lies outside the disc containing D_s
    cnt = 0
    stack = [(w, depth_cap)]
    while stack:
        a, dleft = stack.pop()
        if dleft == 0:
            return False                   # an orbit reached the cap: no certificate
        cnt += 1
        if cnt > budget:
            return False                   # conservative: never a false hole
        c = a * invs
        if abs(c) <= R:
            stack.append((c, dleft - 1))
        c = (a - 2.0) * invs
        if abs(c) <= R:
            stack.append((c, dleft - 1))
        c = (a + 2.0) * invs
        if abs(c) <= R:
            stack.append((c, dleft - 1))
    return True                            # tree exhausted => rigorous hole


def adaptive_depth(r):
    """a depth worth using when probing at distance r from a point of interest"""
    return 50 + int(7 * (-math.log10(max(r, 1e-18))))


# ------------------------------------------------------------------ tiny PNG writer
def write_png(path, W, H, rgb):
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    raw = bytearray()
    for y in range(H):
        raw.append(0)
        raw += rgb[y * 3 * W:(y + 1) * 3 * W]
    with open(path, "wb") as fp:
        fp.write(b"\x89PNG\r\n\x1a\n")
        fp.write(chunk(b"IHDR", struct.pack(">IIBBBBB", W, H, 8, 2, 0, 0, 0)))
        fp.write(chunk(b"IDAT", zlib.compress(bytes(raw), 9)))
        fp.write(chunk(b"IEND", b""))


def raster(x0, y0, x1, y1, RES, depth, path):
    """rasterise the rectangle [x0,x1] x [y0,y1]; black = certified hole"""
    Rc = 1.0 / math.sqrt(2.0)
    W = RES
    H = max(1, int(round(RES * (y1 - y0) / (x1 - x0))))
    buf = bytearray(W * H * 3)
    px = (x1 - x0) / W
    holes = 0
    t0 = time.time()
    for iy in range(H):
        im = y0 + (iy + 0.5) * (y1 - y0) / H
        for ix in range(W):
            re = x0 + (ix + 0.5) * px
            s = complex(re, im)
            p = ((H - 1 - iy) * W + ix) * 3
            a = abs(s)
            if a >= 1.0 or a <= 0.5:
                buf[p] = buf[p + 1] = buf[p + 2] = 245     # settled outside the annulus
            elif is_hole(s, depth):
                holes += 1                                  # leave black
            else:
                buf[p] = buf[p + 1] = buf[p + 2] = 255      # white
            if abs(a - Rc) < px:
                buf[p] = 220; buf[p + 1] = 0; buf[p + 2] = 0
        if iy % 25 == 0:
            print("  row %d/%d  holes=%d  [%.0fs]" % (iy, H, holes, time.time() - t0))
            sys.stdout.flush()
    write_png(path, W, H, buf)
    print("wrote %s   %dx%d   certified holes %d/%d = %.2f%%"
          % (path, W, H, holes, W * H, 100.0 * holes / (W * H)))


# ------------------------------------------------------------------ CLI
def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    mode, A = argv[1], argv[2:]

    def need(n):
        """every mode needs at least n positional arguments; say so rather than crashing"""
        if len(A) < n:
            print(__doc__)
            print("mode '%s' needs at least %d argument%s, got %d"
                  % (mode, n, "" if n == 1 else "s", len(A)))
            sys.exit(1)
    out = None                              # any .png argument is the output file, wherever
    rest = []                               # it appears, so the numeric options stay optional
    for a in A:
        if a.lower().endswith(".png"):
            out = a
        else:
            rest.append(a)
    A = rest

    if mode == "point":
        need(2)
        s = complex(float(A[0]), float(A[1]))
        depth = int(A[2]) if len(A) > 2 else 200
        budget = int(A[3]) if len(A) > 3 else 8000000
        t = time.time()
        h = is_hole(s, depth, budget)
        print("s = %.15f%+.15fi   |s| = %.12f   arg = %.6f deg"
              % (s.real, s.imag, abs(s), math.degrees(cmath.phase(s))))
        print("depth %d:  %s   [%.2fs]"
              % (depth,
                 "HOLE -- certified NOT in M" if h else "in M, or undecided at this depth",
                 time.time() - t))
        return 0

    if mode == "scan":
        need(3)
        c = complex(float(A[0]), float(A[1]))
        r = float(A[2])
        nth = int(A[3]) if len(A) > 3 else 360
        depth = int(A[4]) if len(A) > 4 else adaptive_depth(r)
        t = time.time(); h = 0; egs = []
        for i in range(nth):
            th = 2 * math.pi * i / nth
            if is_hole(c + r * cmath.exp(1j * th), depth):
                h += 1
                if len(egs) < 4:
                    egs.append(round(math.degrees(th), 1))
        print("circle of radius %.3g about %.9f%+.9fi, depth %d: %d/%d certified holes (%.2f%%)"
              % (r, c.real, c.imag, depth, h, nth, 100.0 * h / nth))
        if egs:
            print("  first few at angles (deg): %s" % egs)
        print("  [%.1fs]" % (time.time() - t))
        return 0

    if mode == "arc":
        need(3)
        lo, hi, steps = float(A[0]), float(A[1]), int(A[2])
        R = float(A[3]) if len(A) > 3 else 1.0 / math.sqrt(2.0)
        depth = int(A[4]) if len(A) > 4 else 200
        print("along |s| = %.12f, %d samples in [%g, %g] deg, depth %d"
              % (R, steps, lo, hi, depth))
        nh = 0
        for i in range(steps):
            deg = lo + (hi - lo) * (i + 0.5) / steps
            if is_hole(R * cmath.exp(1j * math.radians(deg)), depth):
                nh += 1
                print("  %.6f deg : HOLE" % deg)
        print("%d/%d samples certified as holes" % (nh, steps))
        return 0

    if mode == "map":
        need(3)
        c = complex(float(A[0]), float(A[1]))
        half = float(A[2])
        RES = int(A[3]) if len(A) > 3 else 360
        depth = int(A[4]) if len(A) > 4 else adaptive_depth(half)
        raster(c.real - half, c.imag - half, c.real + half, c.imag + half, RES, depth,
               out or "holes_map.png")
        return 0

    if mode == "region":
        need(4)
        x0, y0, x1, y1 = float(A[0]), float(A[1]), float(A[2]), float(A[3])
        RES = int(A[4]) if len(A) > 4 else 360
        depth = int(A[5]) if len(A) > 5 else 60
        raster(x0, y0, x1, y1, RES, depth, out or "holes_region.png")
        return 0

    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
