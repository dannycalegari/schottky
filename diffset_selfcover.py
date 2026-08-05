#!/usr/bin/env python3
"""diffset_selfcover.py -- PROTOTYPE of the difference-set self-covering argument.

A possible rigorous route to a CORE point s* in int(M).  Write
D = D_s for the attractor of Phi_s = {z->sz, z->sz-2, z->sz+2}, so that s in M <=> 2/s in D.

    COVERING LEMMA.  If C is compact and C is contained in Phi(C) := sC u (sC+2) u (sC-2),
    then C is contained in D.

(p lies in Phi(C) exactly when one of its backward children p/s, (p-2)/s, (p+2)/s lies in C.)
So the plan is: build C as a rasterised D -- the pixels having a deep backward orbit that
stays in B(0, R_D) -- and then verify
    (a) 2/s* is rho-interior to C,
    (b) every p in C has a backward child that is mu-interior to C (robust self-covering),
    (c) hence a persistence radius delta in s, giving s* in int(M).

HONEST STATUS: this is a prototype, not a proof, and it does not currently succeed.  The
self-covering test in (b) fails on the pixels at the outer boundary of C, whose backward
children leave the raster, so the reported margin mu comes out negative.  Turning this into
a proof needs interval boxes rather than pixels.  It is included because the argument is of
general interest; see the README.

usage
    diffset_selfcover.py <re> <im> [RES]

examples
    diffset_selfcover.py 0.5 0.5 200               s0, the twindragon core
    diffset_selfcover.py 0.25 0.6614378277661477   s1, the tame twindragon core
"""
import cmath, math, sys

if len(sys.argv) < 3:
    print(__doc__)
    sys.exit(1)
s = complex(float(sys.argv[1]), float(sys.argv[2]))
RES = int(sys.argv[3]) if len(sys.argv) > 3 else 500
R = 2.0/(1.0-abs(s))
M = R*1.01; LO,HI=-M,M; W=HI-LO; px=W/RES
def ix(x): return int((x-LO)/px)

# 1) build C = rasterized D via backward-survival DFS per pixel (depth D_MAX)
invs=1.0/s
def survives(w, D=70, budget=400000):
    if abs(w)>R+1e-9: return False
    cnt=0; st=[(w,D)]
    while st:
        a,d=st.pop()
        if d==0: return True
        cnt+=1
        if cnt>budget: return True
        c=a*invs
        if abs(c)<=R: st.append((c,d-1))
        c=(a-2.0)*invs
        if abs(c)<=R: st.append((c,d-1))
        c=(a+2.0)*invs
        if abs(c)<=R: st.append((c,d-1))
    return False
inC=bytearray(RES*RES); nC=0
for iy in range(RES):
    y=LO+(iy+0.5)*px
    for jx in range(RES):
        x=LO+(jx+0.5)*px
        if survives(complex(x,y)): inC[iy*RES+jx]=1; nC+=1
    if iy%80==0: print("  build row %d/%d  |C|=%d"%(iy,RES,nC)); sys.stdout.flush()
print("C = rasterized D: %d pixels (%.1f%%), px=%.4f, R_D=%.3f"%(nC,100*nC/(RES*RES),px,R))

# distance-to-complement (in pixels) via simple multi-pass chamfer (approx)
INF=10**9
dist=[INF]*(RES*RES)
for i in range(RES*RES):
    if not inC[i]: dist[i]=0
for _ in range(2):
    for iy in range(RES):
        for jx in range(RES):
            i=iy*RES+jx
            if inC[i]:
                best=dist[i]
                if jx>0: best=min(best,dist[i-1]+1)
                if iy>0: best=min(best,dist[i-RES]+1)
                if jx>0 and iy>0: best=min(best,dist[i-RES-1]+1.4142)
                if jx<RES-1 and iy>0: best=min(best,dist[i-RES+1]+1.4142)
                dist[i]=best
    for iy in range(RES-1,-1,-1):
        for jx in range(RES-1,-1,-1):
            i=iy*RES+jx
            if inC[i]:
                best=dist[i]
                if jx<RES-1: best=min(best,dist[i+1]+1)
                if iy<RES-1: best=min(best,dist[i+RES]+1)
                if jx<RES-1 and iy<RES-1: best=min(best,dist[i+RES+1]+1.4142)
                if jx>0 and iy<RES-1: best=min(best,dist[i+RES-1]+1.4142)
                dist[i]=best
def interior_depth(w):  # distance (world units) of point w to complement of C ; <=0 if outside C
    a=ix(w.real); b=ix(w.imag)
    if a<0 or a>=RES or b<0 or b>=RES or not inC[b*RES+a]: return -1.0
    return dist[b*RES+a]*px

# 2) 2/s interior margin
w0=2.0/s
rho=interior_depth(w0)
print("2/s = %.4f%+.4fi  interior depth in C = %.4f  (%s)"%(w0.real,w0.imag,rho,"INTERIOR" if rho>0 else "NOT in C"))

# 3) robust self-covering margin: for each p in C, best backward-child interior depth
mu=1e9; worst=None
for iy in range(RES):
    y=LO+(iy+0.5)*px
    for jx in range(RES):
        if not inC[iy*RES+jx]: continue
        p=complex(LO+(jx+0.5)*px,y)
        best=-1.0
        for q in (p*invs,(p-2.0)*invs,(p+2.0)*invs):
            d=interior_depth(q)
            if d>best: best=d
        if best<mu: mu=best; worst=p
print("robust self-covering margin mu = %.4f  (min over C of max-child interior depth)"%mu)
print("  worst pixel p=%.4f%+.4fi"%(worst.real,worst.imag) if worst else "")
if mu>0:
    dsp = mu*abs(s)/ (abs(w0)+2.0)     # rough persistence radius in s (child moves ~|child||ds|/|s|)
    print("  => C is ROBUSTLY self-covering (mu>0) => C subset D. 2/s interior (rho=%.3f>0) => s in int(M)."%rho)
    print("  => rough persistence radius in s ~ %.4f (interval-hardening needed for a proof)."%dsp)
else:
    print("  => self-covering margin <=0 at this resolution (refine / boundary pixels).")
