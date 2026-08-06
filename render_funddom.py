#!/usr/bin/env python3
"""render_funddom.py -- turn a funddom raster into a PNG.

funddom writes two bytes per pixel, (level, gamma): `level` is the least tail length c at
which that limit-trap parameter C was certified.  The two reserved values are the way
round funddom writes them, which is the opposite of what this docstring said until
2026-08-06: **255 = never certified up to cmax, 254 = outside the plotted disc** (in ann
mode, |C| > rho).  Whether C is in T_sigma at all is carried separately, by `gamma`.
This colors them:

    blue    certified early, i.e. at a shallow renormalization depth
    orange  certified only deep
    white   not certified at all up to cmax
    gray    C outside T_sigma

usage
    render_funddom.py <in.bin> <W> <H> <out.png> [flipy] [offset] [--annulus <rho> <sigma_abs> <b>]
    render_funddom.py bar <out.png> [W] [H]        just the color bar

W and H must be the RESX and RESY that funddom was run with; the file carries no header,
so a mismatch silently produces a sheared picture.  `flipy` (0 or 1) flips the vertical
axis, and `offset` shifts the raster horizontally, which is useful for lining a log-mode
picture up with a fundamental domain boundary.

--annulus <rho> <sigma_abs> <b>   (ann-mode rasters only) draws, in red, the two circles
    |C| = rho  and  |C| = rho * sigma_abs**b, both centered at C = 0, using the same
    pixel<->C mapping funddom used to build the raster:
        C_re = -rho + 2*rho*(ix+0.5)/W,   C_im = -rho + 2*rho*(iy+0.5)/H.
    `rho` must match the rho the raster was generated with.  The annulus between these
    two circles is a fundamental domain for the scaling action C -> sigma^b * C under
    which the limit-trap set is invariant, so drawing both boundaries shows how the
    picture tiles the rest of the plane.  Circles are drawn last, on top of everything
    else, with a ~1.5px anti-aliased red ring.  The flag can appear anywhere after the
    positional arguments and does not affect them.

examples
    render_funddom.py s0.bin 700 155 s0.png 0
    render_funddom.py bar bar.png
    render_funddom.py lm.bin 400 400 lm.png 0 --annulus 0.08 0.7071067811865476 1
"""
import struct, zlib, sys

CMIN, CMAX = 0, 30
# Viridis.  Kept in step with figexp::viridis() in figure_export.cc, which is what the
# interactive program's annulus export uses; the two used to be different ramps despite a
# comment saying otherwise, so if you change either, change both and compare them.
# Monotone in lightness, so a deeper tail length always reads as brighter, and it survives
# being printed in grayscale.
RAMP = [(0.0, (0.267004, 0.004874, 0.329415)),
        (0.1, (0.282623, 0.140926, 0.457517)),
        (0.2, (0.253935, 0.265254, 0.529983)),
        (0.3, (0.206756, 0.371758, 0.553117)),
        (0.4, (0.163625, 0.471133, 0.558148)),
        (0.5, (0.127568, 0.566949, 0.550556)),
        (0.6, (0.134692, 0.658636, 0.517649)),
        (0.7, (0.266941, 0.748751, 0.440573)),
        (0.8, (0.477504, 0.821444, 0.318195)),
        (0.9, (0.741388, 0.873449, 0.149561)),
        (1.0, (0.993248, 0.906157, 0.143936))]
WHITE = (255, 255, 255)
GRAY  = (206, 206, 210)
RED   = (237, 20, 20)

def cmap(t):
    if t <= 0: rgb = RAMP[0][1]
    elif t >= 1: rgb = RAMP[-1][1]
    else:
        for i in range(len(RAMP) - 1):
            a, ca = RAMP[i]; b, cb = RAMP[i + 1]
            if a <= t <= b:
                u = (t - a) / (b - a)
                rgb = tuple(ca[k] + u * (cb[k] - ca[k]) for k in range(3)); break
    return tuple(int(round(255 * v)) for v in rgb)

LUT = {c: cmap((c - CMIN) / float(CMAX - CMIN)) for c in range(0, 254)}

def png(path, w, h, rows):
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    raw = bytearray()
    for r in rows:
        raw.append(0); raw += r
    with open(path, "wb") as fp:
        fp.write(b"\x89PNG\r\n\x1a\n")
        fp.write(chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)))
        fp.write(chunk(b"IDAT", zlib.compress(bytes(raw), 9)))
        fp.write(chunk(b"IEND", b""))

def render(binfile, w, h, out, flipy, annulus=None):
    """annulus, if given, is (rho, [radius, ...]): circles |C| = radius are drawn in RED,
    on top of everything else, using funddom's ann-mode pixel<->C mapping (see module
    docstring).  Anti-aliased with a linear falloff over a ring about 1.5px wide."""
    dat = open(binfile, "rb").read()
    assert len(dat) == 2 * w * h, (len(dat), 2 * w * h)
    rows = []
    order = range(h - 1, -1, -1) if flipy else range(h)
    scale_x = scale_y = half_width = rho = None
    radii = ()
    if annulus is not None:
        rho, radii = annulus
        scale_x = 2.0 * rho / w
        scale_y = 2.0 * rho / h
        half_width = 0.75 * max(scale_x, scale_y)   # ~1.5px total ring width, in C-plane units
    for iy in order:
        row = bytearray()
        base = 2 * iy * w
        c_im = (-rho + scale_y * (iy + 0.5)) if annulus is not None else None
        for ix in range(w):
            lev = dat[base + 2 * ix]; gam = dat[base + 2 * ix + 1]
            if lev == 254:   col = WHITE
            elif lev == 255: col = WHITE if gam else GRAY
            else:            col = LUT[lev]
            if annulus is not None:
                c_re = -rho + scale_x * (ix + 0.5)
                dist = (c_re * c_re + c_im * c_im) ** 0.5
                alpha = 0.0
                for r in radii:
                    d = abs(dist - r)
                    if d < half_width:
                        a = 1.0 - d / half_width
                        if a > alpha: alpha = a
                if alpha > 0.0:
                    col = tuple(int(round(RED[k] * alpha + col[k] * (1.0 - alpha)))
                                for k in range(3))
            row += bytes(col)
        rows.append(row)
    png(out, w, h, rows)
    print("wrote", out)

def colorbar(out, W=560, H=26, cmin=CMIN, cmax=CMAX):
    rows = []
    for _ in range(H):
        row = bytearray()
        for ix in range(W):
            c = cmin + (cmax - cmin) * ix / float(W - 1)
            row += bytes(cmap((c - cmin) / float(cmax - cmin)))
        rows.append(row)
    png(out, W, H, rows)
    print("wrote", out)


USAGE = """render_funddom.py -- turn a funddom raster into a PNG.

    render_funddom.py <in.bin> <W> <H> <out.png> [flipy] [offset] [--annulus <rho> <sigma_abs> <b>]
    render_funddom.py bar <out.png>

<in.bin> is what `funddom` writes: two bytes per pixel, (level, gamma), where level is the
least tail length c that certifies the parameter C of that pixel (255 = never certified,
254 = outside the plotted disc) and gamma says whether C lies in T_sigma at all.

  flipy     1 for the `ann` (disc) rasters, whose rows run bottom-up; 0 for `log` ribbons.
            Default 1 if H == W, else 0.
  offset    subtract this from each level before coloring.  Use it when comparing rasters
            taken at different zoom depths: zooming in by k periods raises every level by
            k*b, so passing offset = k*b puts them on a common scale.
  --annulus <rho> <sigma_abs> <b>
            ann-mode only.  Draws, in red, the two circles |C| = rho and
            |C| = rho * sigma_abs**b that bound one fundamental domain of the scaling
            action C -> sigma^b * C.  `rho` must match the raster's rho.  May appear
            anywhere after the positional arguments.

Colors: blue = certified with a short tail, orange = only with a long one, white = not
certified up to the cmax of the run, gray = C not in T_sigma, so no limit trap can exist.
`bar` writes the matching color bar."""


def parse_annulus(argv):
    """Pull an optional '--annulus rho sigma_abs b' out of argv (it may appear anywhere
    after argv[0]) and return (remaining_argv, annulus_or_None), where annulus is
    (rho, [rho, rho * sigma_abs**b]) ready for render()."""
    args = list(argv)
    annulus = None
    i = 1
    out = args[:1]
    while i < len(args):
        if args[i] == "--annulus":
            if i + 3 >= len(args):
                sys.exit("--annulus needs 3 arguments: rho sigma_abs b")
            rho, sigma_abs, b = (float(args[i + 1]), float(args[i + 2]), float(args[i + 3]))
            annulus = (rho, [rho, rho * sigma_abs ** b])
            i += 4
        else:
            out.append(args[i])
            i += 1
    return out, annulus


def main(argv):
    argv, annulus = parse_annulus(argv)
    if len(argv) >= 3 and argv[1] == "bar":
        colorbar(argv[2])
        return 0
    if len(argv) < 5:
        print(USAGE)
        return 1
    binf, W, H, out = argv[1], int(argv[2]), int(argv[3]), argv[4]
    flipy = bool(int(argv[5])) if len(argv) > 5 else (W == H)
    off = int(argv[6]) if len(argv) > 6 else 0
    if off:
        global LUT
        LUT = {c: cmap((max(c - off, 0) - CMIN) / float(CMAX - CMIN)) for c in range(254)}
    render(binf, W, H, out, flipy, annulus=annulus)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
