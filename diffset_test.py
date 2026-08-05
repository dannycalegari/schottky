#!/usr/bin/env python3
"""Is a given difference w an INTERIOR point of the difference set Lambda_s - Lambda_s?

For f(z)=sz-1, g(z)=sz+1 with attractor Lambda:  s in M  <=>  2/s in Lambda-Lambda,
because f(p)=g(q) exactly when p-q = 2/s.  Whether the difference is merely IN the set or
DEEP INSIDE it is the useful distinction: a robust interior margin at w = 2/s is the
difference-set signature of s in interior(M), whereas a frontier point is the signature of a
parameter on the boundary (at the CKW hexahole, for instance, f(Lambda) meets g(Lambda) in a
single point and 2/omega is not interior; at the twindragon parameter they meet in an arc and
2/s is deep interior).

Other differences are worth testing too -- e.g. w = fix(f) - fix(g) = -2/(1-s), the difference
realised by the two fixed points -- so w is a parameter, defaulting to 2/s.

Method: rasterise Lambda from its depth-N point cloud, erode by one pixel to approximate the
interior, then ask whether the interior meets its own translate by w, and if so how far the
difference can be perturbed with the overlap surviving (the robust margin).

usage: diffset_test.py <re_s> <im_s> [N] [RES] [re_w] [im_w]
  N    depth of the point cloud (2^N points; 22 is a good default)
  RES  raster resolution
  w    difference to test, default 2/s
"""
import sys

if len(sys.argv) < 3:
    print(__doc__)
    sys.exit(1)

re_s=float(sys.argv[1]); im_s=float(sys.argv[2])
N  =int(sys.argv[3]) if len(sys.argv)>3 else 22
RES=int(sys.argv[4]) if len(sys.argv)>4 else 380
s=complex(re_s,im_s)
if len(sys.argv)>6: w=complex(float(sys.argv[5]),float(sys.argv[6]))
else:               w=2.0/s               # the difference to test for interiority
R0=1.0/(1.0-abs(s))
# raster window: big enough to hold Lambda (subset of B(0,R0))
M=R0*1.05; LO,HI=-M,M; Wd=HI-LO; px=Wd/RES
occ=bytearray(RES*RES)
sys.setrecursionlimit(100000)
cnt=[0]
def rec(d,P):
    if d==N:
        ix=int((P.real-LO)/px); iy=int((P.imag-LO)/px)
        if 0<=ix<RES and 0<=iy<RES: occ[iy*RES+ix]=1
        cnt[0]+=1; return
    rec(d+1,s*P-1.0); rec(d+1,s*P+1.0)
rec(0,0j)
# erode -> interior
interior=bytearray(RES*RES); nint=0
for iy in range(1,RES-1):
    r=iy*RES
    for ix in range(1,RES-1):
        if not occ[r+ix]: continue
        ok=True
        for dy in(-1,0,1):
            b=(iy+dy)*RES+ix
            if not(occ[b-1] and occ[b] and occ[b+1]): ok=False;break
        if ok: interior[r+ix]=1; nint+=1
print("s=%.5f%+.5fi |s|=%.4f   w=%.4f%+.4fi |w|=%.4f%s"
      %(s.real,s.imag,abs(s),w.real,w.imag,abs(w),
        "  (= 2/s)" if abs(w-2.0/s)<1e-12 else ""))
print("  N=%d cloud=%d RES=%d px=%.4f  occ=%.1f%% interior=%.1f%%"
      %(N,cnt[0],RES,px,100.*sum(occ)/(RES*RES),100.*nint/(RES*RES)))
# test: exist interior p, p' with p-p' = w  (i.e. p' = p - w)
sx=int(round((-w.real)/px)); sy=int(round((-w.imag)/px))   # p' index = p index + (sx,sy)
def ov(ssx,ssy):
    c=0
    for iy in range(RES):
        jy=iy+ssy
        if jy<0 or jy>=RES: continue
        b=iy*RES
        for ix in range(RES):
            if interior[b+ix]:
                jx=ix+ssx
                if 0<=jx<RES and interior[jy*RES+jx]: c+=1
    return c
base=ov(sx,sy)
# raw (non-eroded) overlap: do Lambda and Lambda+2/s overlap at all (touch region size)?
def ovraw(ssx,ssy):
    c=0
    for iy in range(RES):
        jy=iy+ssy
        if jy<0 or jy>=RES: continue
        b=iy*RES
        for ix in range(RES):
            if occ[b+ix]:
                jx=ix+ssx
                if 0<=jx<RES and occ[jy*RES+jx]: c+=1
    return c
print("  RAW-overlap (occ) for difference w: %d cells"%ovraw(sx,sy))
print("  interior-overlap for difference w: %d cells"%base)
if base==0:
    print("  => w is NOT interior at this resolution (boundary/frontier point of Lambda-Lambda)")
else:
    Rmax=0
    for R in range(1,25):
        allok=True
        for ex in range(-R,R+1):
            for ey in range(-R,R+1):
                if max(abs(ex),abs(ey))!=R: continue
                if ov(sx+ex,sy+ey)==0: allok=False;break
            if not allok: break
        if allok: Rmax=R
        else: break
    print("  => w INTERIOR to Lambda-Lambda, robust margin >= %d px = %.4f world units"%(Rmax,Rmax*px))
    print("     (Rmax hit scan cap)" if Rmax>=24 else "     (margin resolved)")
