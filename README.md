# schottky

Copyright 2014–2026 by Danny Calegari and Alden Walker.
Released under the GPL (see `LICENSE`). See `INSTALL` for installation instructions.

`schottky` is a tool for studying the parameter space of two-generator iterated function
systems of complex dilations of the plane. `doc/schottky-doc.pdf` is the manual: it documents
the interactive program, and its section "The headless tools" is a full reference for
`certify_arc` and `funddom` --- every subcommand and mode, with worked examples. This file is
the quick tour; the manual is the reference.

For a complex parameter `s` with `0 < |s| < 1` put

    f(z) = s z - 1,      g(z) = s z + 1,

let `Lambda = Lambda_s` be the attractor of `{f, g}`, and let

    M = { s : Lambda_s is connected }

be the connectedness locus. `M` is the object of Bandt's conjecture and of Calegari–Koch–Walker,
*Roots, Schottky semigroups, and a proof of Bandt's conjecture* (arXiv:1410.8542), whose
terminology (traps, trap-like balls, renormalization points, limit traps) is used throughout.

Alongside the interactive program this repository now contains **headless, scriptable tools**
for deciding membership in `M` and for certifying that regions of parameter space lie in its
interior:

| | |
|---|---|
| `certify_arc` | traps and rigorous interval certificates, at a point, along an arc, or over a box. Most subcommands take an angle on the circle `\|s\| = 1/sqrt2`, which is the arc the paper works on; `pointz` and `dumptlbmanyz` take the parameter directly and so reach the rest of `M` |
| `funddom` | limit traps over a fundamental domain of the elliptic curve `C*/<sigma^b>`, at a built-in core, any finite coincidence, or any enumerated landmark point |
| `holes.py` | rigorous certificates that a parameter is **not** in `M`; scans and rasters |
| `diffset_test.py` | is a given difference interior to `Lambda - Lambda`? |
| `spiral.py` | follow the renormalization spiral at a renormalization point |
| `prune_tlb.py` | reduce a trap-like ball list to its maximal elements |
| `render_funddom.py` | turn a `funddom` raster into a PNG |
| `diffset_selfcover.py` | prototype of the difference-set self-covering argument |
| `figure_export.h/.cc` | write a figure as PNG, EPS or PDF; raster plus true vector overlays |

Every program prints its own usage when run with no arguments, and the Python tools carry their
mathematical explanation in the module docstring.

## Building

    make                # schottky, certify_arc, funddom
    make schottky       # the interactive GUI (needs X11)
    make headless       # certify_arc and funddom, no X11 required
    make certify_arc    # headless trap driver
    make funddom        # limit traps (plain C)

The Python tools need only the standard library — no numpy, no matplotlib. They write PNGs
directly.

### Building without X11

Only `schottky` needs a display. `certify_arc` and `funddom` do not need X11 at all — neither the
headers nor `libX11` — so `make headless` works on a cluster node or in a container with nothing
installed but a C and a C++ compiler (`funddom` is C99, the rest is C++).

The trap code is shared between the two, and it contains debug-drawing calls: a handful of places
that open a window on the balls, the grid or the convex hull when `verbose > 0`. Every one of them
is a local drawing variable inside a verbose branch, never part of the computation. `make headless`
therefore compiles `ifs.cc`, `trap_grid.cc` and `movie.cc` a second time with `-DIFS_NO_GRAPHICS`
(into `ifs_nogfx.o`, `trap_grid_nogfx.o`, `movie_nogfx.o`), which

* omits `#include "graphics.h"`, and with it the X11 headers;
* drops the debug windows and the `draw_*` / `user_interface` members of the `ifs` class;
* replaces `TrapGrid::show`, `show_connected_components`, `show_distance_functions` and
  `ifs::draw_set_B_balls` — which exist only to open a window — by stubs that print a one-line
  note, so a verbose headless run says why no picture appeared instead of failing to link.

Computed results are unaffected: every `certify_arc` subcommand gives byte-identical output either
way. If you add to the trap code, keep new drawing inside `#ifndef IFS_NO_GRAPHICS`; the check is
that `make headless` still succeeds when the X11 headers are unreachable.

The two headless programs are documented command by command in the section "The headless tools" of
`doc/schottky-doc.pdf`, which also spells out which of them produce rigorous interval-arithmetic
certificates and which only evaluate a rigorous criterion in floating point.

## Three questions, and the tool for each

### Is this parameter in `M`?

`s in M` if and only if `f(Lambda)` meets `g(Lambda)`, if and only if `2/s` lies in the
difference set `D_s = Lambda - Lambda`, which is the attractor of `z -> sz`, `z -> sz ± 2`.
So `2/s` lies in `D_s` exactly when it admits an infinite *backward* orbit under
`a -> a/s, (a±2)/s` staying inside the disc `|z| <= 2/(1-|s|)` that contains `D_s`. If that
ternary backward tree can be exhausted, no such orbit exists and `s` is **certifiably not** in
`M` — a hole:

    python3 holes.py point 0.0227 0.703          # HOLE: certified not in M
    python3 holes.py point 0.5 0.5 200           # survives: in M, or undecided

The asymmetry matters and is deliberate: death of the tree is a finite proof of
non-membership, while survival to the depth cap is only the absence of a certificate. Scanning
and rasterising use the same test:

    python3 holes.py scan   0.5 0.5 1e-3         # any holes on a small circle about 1/2+i/2?
    python3 holes.py arc    0 90 400             # sweep the circle |s| = 1/sqrt2
    python3 holes.py region 0 0 1 1 120 60 quadrant.png    # M in the first quadrant
                                                           # (RES 400 costs ~2 hours)

### Is this parameter in the *interior* of `M`?

That needs a *trap*: four points in cyclic order on the outer boundary of `Lambda`, alternately
in `f(Lambda) - g(Lambda)` and `g(Lambda) - f(Lambda)`. A trap is an open condition, so
exhibiting one certifies a whole disc of parameters. `certify_arc` searches for traps and
verifies them, at a point or over an arc:

    ./certify_arc point  60 16 18 1e-4           # certified=1, word length m=14
    ./certify_arc sweep  10 80 200 16 18 1e-4    # sweep an arc of |s| = 1/sqrt2
    ./certify_arc pointz 0.6 0.3 18 16 1e-9      # anywhere in M, not just on that circle

The circle `|s| = 1/sqrt2` is where `dM` has its interesting part, and it is what most of these
subcommands take an angle on — but `M` is two-dimensional and nothing in the trap machinery needs
`s` to lie on it. `pointz` and `dumptlbmanyz` take the parameter directly. The interactive
program's `Trap` overlay works anywhere too, and draws the certifying pair.

The interval-arithmetic versions in `rigor.h` certify a closed *box* of parameters rather than a
floating-point point sample, using outward rounding, but they divide up differently: `rigorcover`
**searches**, covering an arc with certified boxes and checkpointing to `cov_<id>.jsonl` so a long
run can be resumed, whereas `rigorcert` and `rigorcertseg` **verify** a `(u,v)` you supply — the
latter against a relaxed 4-arc segment-persistence criterion rather than the same one in intervals.
`dumptlb` prints the trap-like balls at a parameter; `dumptlbmany` prints far more of them (see
below).

### What happens near a renormalization point?

At a *renormalization point* `sigma` the combinatorics repeat: for words `u, v` of length `a`
and `s, t` of length `b` with `pi(u s^inf, sigma) = pi(v t^inf, sigma)`, the parameters

    sigma + C sigma^{bn}

spiral into `sigma`, and (CKW Lemma 9.2.5) if the vector

    V(C; x, y) = sigma^{-a-c}(p_u - p_v) + sigma^{-c}(p_x - p_y) + sigma^{-a-c} C P'(sigma)

is trap-like at `sigma` for some words `x, y` of length `c`, then `u s^n x`, `v t^n y` trap
`sigma + C sigma^{bn}` for every large `n`: the parameter `C` *admits a limit trap*.

That condition is invariant under `C -> sigma^b C`, so it lives on the elliptic curve
`E_sigma = C*/<sigma^b>` — and covering a single fundamental domain of `E_sigma` gives *every*
`C`, hence a punctured neighborhood of `sigma` inside `M`, hence `sigma` in the interior of
`M`. `funddom` computes exactly that. Internally it uses the reduction

> `C` admits a certified limit trap **iff** the backward orbit tree of
> `Y(C) = -(p_u - p_v + C P'(sigma))/sigma^a` under the difference-set IFS
> `{z -> sigma z, sigma z ± 2}` reaches the trap-like set,

which makes each parameter cost `~(3/2)^cmax` tree nodes instead of a word search.

A worked example — the twindragon core `sigma = 1/2 + i/2`:

    ./certify_arc dumptlbmany 45 22 1e-9 40 14 8 > tlb45.txt
    python3 prune_tlb.py tlb45.txt tlb45_pruned.txt
    ./funddom s0 log 24 700 155 0.08 0 s0.bin 1 tlb45_pruned.txt
    python3 render_funddom.py s0.bin 700 155 s0.png 0

reports the fraction of one fundamental domain covered, and writes a raster of it colored by
the tail length `c` needed.

The **CKW hexahole** is the control worth running, since it is known to lie on the *boundary* of
`M` and so must **not** reach full coverage. It needs its own command line, not a substitution of
`hex` for `s0`: it sits at `|s| = 0.6388`, off the circle `|s| = 1/sqrt2` that `dumptlbmany`
parameterises by angle, so its balls come from `dumptlbmanyz`, which takes the parameter directly;
and `|dY/dC|` there is 322 rather than 11, so `rho = 0.08` would put the entire raster outside the
disc that contains the difference set and report a vacuous 0%.

    ./certify_arc dumptlbmanyz 0.371858680074136 0.519411153747943 20 1e-9 40 14 8 > tlbhex.txt
    python3 prune_tlb.py tlbhex.txt tlbhex_pruned.txt
    ./funddom hex log 24 400 90 0.0003 26 hex.bin 1 tlbhex_pruned.txt
    # covered 99.7389%   in T_sigma 99.8806%

The coverage stalls just below the `T_sigma` ceiling and stays there however far `cmax` is pushed:
the uncovered part is a scale-invariant set, which is CKW's hole spiral.

**A finite coincidence is not one of these points.** It is worth being explicit, because the
definition above admits it degenerately: if `u(0) = v(0)` then `u` and `v` are the *same* affine
map, so taking `s = t = f` makes `pi(u s^inf) = pi(v t^inf)` hold trivially, with `a = |u|`,
`b = 1` and `Delta = 0`. Nothing asymptotic is happening there, though — the agreement is already
exact — so there is no family `sigma + C sigma^{bn}` to trap and no fundamental domain of
`E_sigma` to cover, and `funddom` refuses such a `sigma` as a core for a coverage run. Those
parameters are *roots*, they are interesting in their own right, and they are covered under
[Landmark points](#landmark-points-where-the-limit-trap-mechanism-applies) below.

Three practical points, all learned the hard way:

* **A ball set belongs to its parameter.** Trap-like balls are only valid at the `s` they
  were computed at, so `dumptlb`/`dumptlbmany`/`dumptlbmanyz` write a `# s = <re> <im>` line,
  `prune_tlb.py` carries it through, and `funddom` refuses to run if it disagrees with its own
  `sigma`. Without that check, pairing a ball file with the wrong core gives a confident,
  plausible-looking and entirely meaningless coverage figure. If you write your own producer
  or consumer of these files, keep the line.

* **Normalisation.** `certify_arc` works in `f(x) = zx`, `g(x) = z(x-1)+1` with base point
  `1/2`, whereas the formula above is in `f = sz-1`, `g = sz+1` with base point `0`. The two
  are conjugate by `x -> (2x-1)/(1-s)`, so trap-like *vectors* differ by a factor `2/(1-s)`
  (modulus 2.83 at `1/2+i/2`). `funddom` applies this to the ball file it reads; if you write
  your own consumer of `dumptlb` output, do not forget it.
* **How many trap-like balls.** `ifs::trap_like_balls_from_balls` samples the 5 largest gaps of
  the convex hull at 3 positions each and keeps one ball per gap, so `dumptlb` returns at most
  10. That is far too few to cover anything. `dumptlbmany` emits as many admissible balls of
  CKW Def. 8.2.3 as you ask for — the `<ngaps>` largest hull gaps, `<ntrials>` positions along
  each, `<nradial>` depths down the inward normal — typically 1000–1500 after `prune_tlb.py`,
  with radii about 2.3 times larger. The difference decides
  whether the fundamental domain closes up or stalls in the nineties.

## Landmark points: where the limit-trap mechanism applies

A **landmark point** is a renormalization point: a `sigma` for which there are words `u, v` of
length `a` and `s, t` of length `b` with `pi(u s^inf, sigma) = pi(v t^inf, sigma)`. That is exactly
the hypothesis of CKW Lemma 9.2.5, so these are the parameters where the limit-trap machinery
above applies — `s_0`, `s_1` and the CKW hexahole are three of them.

They can be enumerated. Writing `d_j = (eps^u_j - eps^v_j)/2` in `{0,+-1}` for the difference of
the two coding sequences, `d` is some `A` of length `a` followed by some `B` of length `b` repeated
forever, and `sum_j d_j z^j = A(z) + z^a B(z)/(1 - z^b)`. Clearing the denominator, the coincidence
condition is

    Q(z) := A(z)(1 - z^b) + z^a B(z) = 0,     A_0 = -1,   A, B with coefficients in {0,+-1},

and everything a coverage run needs follows from `Q`:

    Delta = p_u - p_v = 2 A(sigma),     P(z) = 2 Q(z)/(1 - z^b),
    and since Q(sigma) = 0,             P'(sigma) = 2 Q'(sigma)/(1 - sigma^b).

Note this reaches the **infinite** (eventually periodic) coincidences, which `coin:` does not: all
three built-in cores are of that kind. So:

    ./funddom landmarks 8

which lists every landmark point of complexity `a + b <= 8`, one per line, with `sigma`, `a`, `b`,
the degree, `Delta`, `P'(sigma)`, and a spec that feeds straight back in. There are 1, 18, 99, 533,
2421, 10958 and 46201 of them for `Nmax = 3..9`.

`landmarks` is exhaustive, so it costs `3^Nmax` polynomials and passes 100,000 points at
`Nmax = 10`.  For a deep zoom -- where a hole spiral accumulates on one landmark of high
complexity -- ask for that neighbourhood instead:

    ./funddom lmnear 0.461537666483585 0.459125328106272 1e-5 12

which searches rather than enumerates, pruning coefficient prefixes against a geometric tail
bound (396 leaves instead of 531,441 polynomials, in that example).  It never discards a
landmark inside the disc, and agrees exactly with `landmarks` for `Nmax <= 9`.

For a **finite** coincidence `u(0) = v(0)` -- a root rather than a landmark -- use the `coin:`
selector, which solves for `s` as a root of `sum_j d_j s^j` with `d_j = (eps^u_j - eps^v_j)/2`:

    ./funddom coin:fgff:gfgg

    coincidence u = fgff , v = gfgg   (m = 4)
      d = (-1,+1,-1,-1)
      polynomial (degree 3):  -1 + z - z^2 - z^3  =  0
      admissible roots (1/2 < |sigma| < 1, Im sigma > 0): 1
       [0] sigma = 0.419643377607081+0.606290729207199i   |sigma| = 0.737352705760   arg = 55.311003 deg
           a = 4, b = 1, s = t = f;  Delta = 0;  P'(sigma) = 1.470353793-5.478273590i

If the polynomial has several admissible roots, `coin:<u>:<v>:<k>` names the `k`-th. This is a
**query only**: it reports the parameter, and is refused as a core for a coverage run, for the
reason given under [What happens near a renormalization point?](#what-happens-near-a-renormalization-point)
above. Note also that no finite coincidence can lie on `|s| = 1/sqrt2` at all, since a `{0,±1}`
polynomial with nonzero constant and leading coefficients is monic up to sign, making its roots
algebraic integers, while `|s|^2 = 1/2` would force `s conj(s) = 1/2` to be a rational algebraic
integer.

Both are reachable from the interactive program's `uv =` box, which dispatches on whether the
words carry brackets (eventually periodic, so landmarks) or not (finite, so roots).

The finite coincidences can also be enumerated wholesale, exactly as the landmarks can:

    ./funddom roots 8                                   # every root of degree <= 8 (7520)
    ./funddom rootsnear 0.4615377 0.4591253 1e-4 20     # only those near a point

`roots` costs `3^maxdeg` polynomials (degree 10 is 95,776 roots and ~16s); `rootsnear` prunes
coefficient prefixes the same way `lmnear` does, so degree 20 near a point takes 148 leaves
instead of the 3.5 billion polynomials an exhaustive pass would need.  Both are the `Roots:`
layer of the interactive program, which picks between them by how far you are zoomed in.

The three cores are marked in the output:

    0.500000000000000 0.500000000000000 ... 2 1 2 ... lm:-+:-:0          <== s0 (twindragon)
    0.250000000000000 0.661437827766148 ... 1 4 4 ... lm:-:+--+:0        <== s1 (tame twindragon)
    0.371858680074136 0.519411153747943 ... 8 1 8 ... lm:-+---+++:-:1    <== hex (CKW hexahole)

That last column is a core selector in its own right, so a landmark can be handed straight to a
coverage run — `lm:<A>:<B>[:k]`, with `A` and `B` as sign strings over `{-,0,+}` and `k` choosing
among the admissible roots when the polynomial has several (the hexahole is the second root of
`-1 + 2z - 2z^2 + 2z^5 - 2z^8`, hence the `:1`). Run `funddom lm:<A>:<B>[:k]` with no further
arguments to see the polynomial and list its roots before committing to a raster.

    ./funddom lm:-+---+++:-:1 log 24 400 90 0.0003 26 out.bin 1 tlbhex_pruned.txt

reproduces the built-in `hex` core **byte for byte** — which is the regression test that keeps the
formulas above honest, along with the requirement that `landmarks` recover all three cores at their
correct `(a, b)` with the `Delta` and `P'(sigma)` the `s0`/`s1`/`hex` branches hard-code.

### The asymptotic picture at a landmark

`funddom`'s `ann` mode rasterises `C` over the square `[-rho, rho]^2`, and the set of limit-trap
parameters is invariant under `C -> sigma^b C`, so the annulus `rho|sigma|^b <= |C| <= rho` is one
**fundamental domain** of `E_sigma`. `render_funddom.py --annulus <rho> <|sigma|> <b>` draws its two
bounding circles in red, which is what makes the self-similarity legible — the whole punctured disc
is tiled by rescaled copies of what lies between them:

    ./funddom lm:-+:-:0 ann 22 500 500 0.08 0 s0.bin 1 tlb45_pruned.txt
    python3 render_funddom.py s0.bin 500 500 s0.png 0 --annulus 0.08 0.7071067811865476 1

Blue means certified at a shallow tail length `c`, paler and then orange means only deep, white
means uncertified up to `cmax`, gray means outside `T_sigma`. Since `|sigma| = 1/sqrt2` and `b = 1`
at the twindragon core, the inner circle sits at `0.08/sqrt2 = 0.0566`. Swap in any spec from
`funddom landmarks` — remembering that `|sigma|` and `b` come from that same line — to get the
corresponding picture at any other landmark point.

## The difference set directly

`diffset_test.py` rasterises `Lambda`, erodes it, and asks whether a given difference `w` is an
interior point of `Lambda - Lambda`, reporting a robust margin:

    python3 diffset_test.py 0.5 0.5 20 300            # w = 2/s, the default
    python3 diffset_test.py 0.5 0.5 20 300 -2 -2      # w = fix(f) - fix(g) = -2/(1-s)

A robust margin at `w = 2/s` is the difference-set signature of `s` in the interior of `M`; a
frontier point is the signature of a boundary parameter. `diffset_selfcover.py` is a prototype
of the corresponding *proof* strategy (a compact `C` with `C` contained in `sC ∪ (sC±2)` is
contained in `D_s`); it is included because the argument is of general interest, but it is
honestly a prototype — the self-covering test fails on boundary pixels and would need interval
boxes to become a proof.

## Exporting figures

`figure_export.h/.cc` writes a figure as **PNG, EPS or PDF**. It exists because a picture from
this program is two things at once: a raster, computed one pixel at a time by a slow recursion,
and a set of *curves* drawn on top of it — marked points, the uv-graph, a trap, the circle
`|s| = 1/sqrt2`, axes, labels. In a paper the curves should stay curves, sharp at any
magnification and restylable by whoever is writing the paper. So a `Figure` is a raster **plus**
a list of vector overlays, and the EPS and PDF back ends keep them apart: the raster is embedded
as an image, the overlays are emitted as real paths. PNG, being a pure raster format,
rasterises the overlays (with antialiasing, and with its own 5x7 font for labels).

Everything is in **mathematical coordinates** — the complex plane, `y` upwards, the window given
by its corners. Nothing but `Raster` knows about pixels, and `Raster`'s row 0 is the bottom row,
so a caller never has to think about the y-flip that X11 imposes.

    figexp::Figure F;
    F.set_window(ll.real(), ll.imag(), ur.real(), ur.imag());
    F.raster.set_size(1000, 1000);
    for (...) F.raster.set_pixel(i, j, color);            // row 0 = bottom
    F.add_circle(0, 0, 1/sqrt(2), figexp::Style::stroke(1, 0, 0, 1.5));
    F.add_dot(s.real(), s.imag(), figexp::Style::fill(0, 0, 1));
    F.add_text(x, y, "s = 1/2 + i/2", figexp::Style::fill(0, 0, 0), 9, 1);

    figexp::Options opt;                  // defaults are sensible
    opt.width_pt = 360;                   // 5 inches across
    std::string err;
    if (!figexp::write_auto(F, opt, "figure.pdf", &err)) std::cerr << err << "\n";

`write_auto` takes the format from the extension; `write_figure` takes it from `opt.format`. The
options are the physical width in points, the raster resolution (`raster_px`, resampling if it
differs from the raster's own size), whether overlays are vector or burnt into the image
(`vector_overlays`), the background color, and whether to draw a frame and the axes. The
aspect ratio always comes from the mathematical window, never from the output size, so a figure
cannot be silently stretched.

**No new dependencies.** PNG and PDF need DEFLATE, which is implemented here (fixed Huffman with
LZ77, falling back to stored blocks when that would expand the data) rather than by linking
zlib — the program's only library dependency is still X11, and `figure_export.o` itself needs
neither X11 nor zlib, so the headless tools can use it too. `deflate_raw`, `zlib_wrap` and
`write_png` are exposed for other code in the tree that wants to write a PNG.

The GUI uses this: the **Write picture** button recomputes the parameter-space picture at the
chosen resolution and writes `schottky_mand_NNN.png`, `.eps` and `.pdf` plus a `.txt` sidecar
recording the window, the parameter and every enabled layer with its depth, so the figure can be
regenerated later. `NNN` increases each time, so nothing is overwritten. In the vector formats the
overlays — the circle `|s| = 1/sqrt2` and the current parameter — are real paths, not pixels.

`figure_export_test.cc` is a self-test: it round-trips the compressor on adversarial inputs,
writes one of every kind of figure, and checks that degenerate input is refused rather than
crashing. Since the compressor is hand-rolled, it is worth running once on a new machine:

    make test

The generated `test_out/*.png` can be checked with any PNG reader, and `test_out/*.{eps,pdf}`
with `gs -dNOPAUSE -dBATCH -sDEVICE=nullpage` if you have Ghostscript.

## Testing `rigor.h`

`rigor.h` is the one file here whose entire value is soundness. Everything else can be wrong in a
way that shows up as a bad picture; this can be wrong in a way that shows up as a theorem that is
not true. So it has its own known-answer test:

    make rigortest

Three layers, because each catches something the others cannot.

* **Hand-computed answers.** The exact range of a product of intervals, `|[-3,2]| = [0,3]`, that
  `cos_iv` finds the minimum at an *interior* multiple of pi rather than only at the endpoints, and
  that `PI_I` contains the true pi — which it does, but only just: the nearest double to pi is
  1.22e-16 away and one ulp there is 4.44e-16, so the outward widening is load-bearing.
* **The containment property, against an independent oracle.** `rigor_test.cc` dumps every interval
  operation it performs — inputs and computed output, in hex float so nothing is lost — and
  `rigor_test.py` re-derives each result in **exact rational arithmetic** and checks that the
  computed interval contains it. A verifier written in the same style as the code under test would
  share any misunderstanding of what outward rounding is for; `fractions.Fraction` does not. About
  3200 operations, over inputs chosen to be awkward: straddling zero, wildly different magnitudes,
  exact powers of two where the `nextafter` step size changes, and numbers with no finite binary
  expansion.
* **The closed form against a different implementation.** `rigor.h` computes a word's disk centre
  analytically; `ifs.cc` computes the same disk by composing the maps. Interval testing cannot
  catch a wrong closed form, because the float and interval versions would be wrong together, so
  the two are compared directly.

**The suite was validated by breaking the code on purpose**, which is the only way to know a test
is not decorative. Three mutants, each caught: removing the outward rounding from `radd` (caught by
the exact-arithmetic oracle); making `certify_box` accept a non-positive margin, which is the most
dangerous error possible here since it turns "could not verify" into "CERTIFIED" (caught by two
checks); and flipping a sign in the closed form for the word centre (caught by three, including the
cross-check against `ifs.cc`).

That second mutant is worth dwelling on: an earlier draft of this test did **not** catch it. Every
refusal it checked — a word linked with itself, two words on the same side, a 90-degree box — is
rejected by the *floating-point* pre-check, so the interval stage never ran. Covering the dangerous
path needs a case that is a genuine float trap on a box too wide for the intervals to separate,
where the only correct answer is `interval-FAILED`. Measured on the real pair, a box of width 0.002
about 83.6255 certifies and 0.005 does not, so the test uses 0.01.

## Authorship

The interactive program and the underlying IFS, trap and trap-like-ball machinery
(`ifs*.cc`, `graphics*`, `trap_grid*`, `movie*`, `schottky.cc`) are by Danny Calegari and Alden
Walker, 2014.

The headless tools added in 2026 — `certify_arc.cc`, `rigor.h`, `funddom.c`, `holes.py`,
`diffset_test.py`, `diffset_selfcover.py`, `spiral.py`, `prune_tlb.py`, `render_funddom.py` —
were written by Claude (Anthropic) in collaboration with Danny Calegari, in the course of the
computer-assisted parts of *Laminations and external angles for similarity pairs* by Danny
Calegari and Alden Walker.

Sixteen files of the 2014 program were modified by that work, and — unlike an earlier draft of
this note, which claimed otherwise — some things were **removed**. For the record:

* Three features were deleted, at Danny Calegari's direction, because they were not sound.
  **Set C** tested a resolution-dependent proxy for `f(Lambda)` and `g(Lambda)` meeting in a
  single point, so the set it drew changed as you zoomed. **Theta** plotted the coordinate
  `theta` from `compute_coordinates`, which was measured against a newer, trusted pipeline and
  found right to about 1% at the median parameter but badly wrong at roughly one parameter in
  seven, with no indication of failure — and not convergent in depth at the bad points.
  **Dirichlet** was a third unfinished experiment. The backing functions `close_to_set_C` and
  `compute_coordinates` went with them (about 306 lines of `ifs.cc`), as did the parameter-space
  layers and the "Find coords along path" tool built on them.
* The **Limit traps** checkbox was removed too. It called `check_limit_TLB`, which hard-codes the
  CKW hexahole and overwrites the caller's parameter with it, so it answered a question about
  that one point no matter where you asked. The general version of the same mathematics is
  `funddom` and its C API `funddom_core.h`, reached from the interactive program by the
  "Limit traps in annulus" button at any landmark point. `check_limit_TLB` itself is kept, and
  commented, as the reference implementation of the CKW section 9 argument.
* `graphics_old.cc` was deleted. It was dead: nothing referenced it, it no longer compiled
  against the current headers, it was in no makefile rule, and its own header read 12/17/2000.
  It remains in the git history.
* `heuristic_convex_hull` and `convex_hull_recurse` were deleted along with the hull they
  belonged to (see below).

Everything else present in 2014 is still present and still called from the same places. The
substantive additions to the 2014 files:

* `ifs_trap.cc` gains `trap_interleaves_topological` (a purely topological arc-alternation
  test, which is the criterion that remains usable on the marginal circle `|s| = 1/sqrt2`),
  `ckw_point_certificate` (CKW Def. 7.1.3 as a point certificate with a rigorous margin), and
  `find_trap_mixed` (trap search over word pairs of *unequal* length).
* `ifs_trap_like.cc` gains `check_TLB_mixed` and `check_TLB_bestfirst` — the latter a
  level-synchronous beam variant of `check_TLB` which explores the pair tree best-first while
  re-validating every hit with the identical containment test and the identical `eps` formula,
  so it is exactly as sound. It also **replaces `convex_hull`**: the 2014 routine combined a
  heuristic hull with a divide-and-conquer refinement, and it is now a monotone-chain (Andrew)
  hull with exact duplicate-point removal. This matters because `ball_convex_hull` feeds the
  hull to the trap-like ball construction, where a wrong hull silently corrupts the balls. The
  superseded `heuristic_convex_hull` and `convex_hull_recurse` have been deleted.
* `ifs.h` gains the corresponding declarations.

## References

* D. Calegari, S. Koch, A. Walker, *Roots, Schottky semigroups, and a proof of Bandt's
  conjecture*, arXiv:1410.8542. Traps (§7), trap-like balls (§8), renormalization points and
  limit traps (§9).
* C. Bandt, *On the Mandelbrot set for pairs of linear maps*, Nonlinearity 15 (2002).
* B. Solomyak, *On the 'Mandelbrot set' for pairs of linear maps: asymptotic self-similarity*,
  Nonlinearity 18 (2005).
* D. Calegari, A. Walker, *Laminations and external angles for similarity pairs*, in
  preparation — the source of the limit-trap and fundamental-domain tooling here.
