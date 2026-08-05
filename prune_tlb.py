#!/usr/bin/env python3
"""prune_tlb.py -- reduce a trap-like ball list to a non-redundant sub-list.

    certify_arc dumptlbmany  45 22 1e-9 40 14 8 > tlb_45.txt
    python3 prune_tlb.py tlb_45.txt tlb_45_pruned.txt
    funddom s0 log 30 1400 309 0.08 0 out.bin 4 tlb_45_pruned.txt

`certify_arc dumptlbmany` emits one ball for every (hull gap, position along the gap, depth
along the inward normal) it tries, so the list is highly redundant: many balls are contained
in a larger one and contribute nothing to the union.  Since `funddom` pays a per-node cost in
the number of balls, pruning to the maximal ones is worth doing once.

The prune is the obvious greedy one -- sort by radius descending, keep a ball only if it is
not contained in a ball already kept -- which is exact for "is this ball redundant against a
single other ball", and conservative (it never drops a ball that enlarges the union) though it
does not detect a ball covered by a *union* of several others.  Typical reduction is 2-3x.

Input format (as emitted by dumptlbmany; '#' lines are ignored):

    center_re center_im radius

Output is the same format, sorted by decreasing radius, with a header comment recording the
provenance.
"""
import sys, math


def read_balls(path):
    """returns (balls, s_line) where s_line is the '# s = re im' provenance line, if any.

    That line says which parameter the balls were computed at.  It has to survive
    pruning, because funddom checks it against its own sigma -- a ball set used at the
    wrong parameter produces a plausible but meaningless coverage number."""
    balls = []
    s_line = None
    for line in open(path):
        line = line.strip()
        if line.startswith('# s ='):
            s_line = line
            continue
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            x, y, r = float(parts[0]), float(parts[1]), float(parts[2])
        except ValueError:
            continue
        if r > 0:
            balls.append((x, y, r))
    return balls, s_line


def prune(balls):
    """keep the maximal balls: drop b if some kept ball contains it"""
    balls = sorted(balls, key=lambda t: -t[2])
    kept = []
    for x, y, r in balls:
        redundant = False
        for a, b, R in kept:
            if R < r:                       # kept list is sorted, so R >= r always
                break
            d = math.hypot(x - a, y - b)
            if d <= R - r:                  # closed containment
                redundant = True
                break
        if not redundant:
            kept.append((x, y, r))
    return kept


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    src = argv[1]
    dst = argv[2] if len(argv) > 2 else None
    balls, s_line = read_balls(src)
    if not balls:
        sys.stderr.write('no balls read from %s\n' % src)
        return 1
    kept = prune(balls)
    area_in = sum(math.pi * r * r for _, _, r in balls)
    area_out = sum(math.pi * r * r for _, _, r in kept)
    sys.stderr.write('%s: %d balls -> %d maximal (%.1fx), max radius %.6f, '
                     'sum of areas %.5f -> %.5f\n'
                     % (src, len(balls), len(kept), len(balls) / max(len(kept), 1),
                        kept[0][2], area_in, area_out))
    out = sys.stdout if dst is None else open(dst, 'w')
    out.write('# maximal trap-like balls, pruned from %s (%d -> %d)\n'
              % (src, len(balls), len(kept)))
    if s_line:                      # carry the parameter through; funddom checks it
        out.write(s_line + '\n')
    else:
        sys.stderr.write('%s: warning: no "# s = re im" line, so the pruned file does not\n'
                         '  record which parameter these balls belong to and funddom cannot\n'
                         '  check it.  Regenerate the dump with a current certify_arc.\n' % src)
    out.write('# center_re center_im radius\n')
    for x, y, r in kept:
        out.write('%.9f %.9f %.9f\n' % (x, y, r))
    if dst is not None:
        out.close()
        sys.stderr.write('wrote %s\n' % dst)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
