// rigor.h -- interval-arithmetic trap certificates.  Header-only, uses
// std::complex<double> for the float boundary combinatorics and a small sound
// interval type (outward rounding via std::nextafter) for the certificate.
//
// Certificate (CKW Def 7.1.3): words u,v; the outer boundary of the union of
// the two ball-clouds covering u.Lambda, v.Lambda alternates u/v arcs >=2 each
// (linking).  certify_box verifies, in interval arithmetic over an s-BOX, that
// every outer-loop arc's midpoint lies strictly outside the whole OPPOSITE
// cloud -- so the u/v tag of each boundary piece, hence the alternation, is
// invariant over the box => the whole box is in int(M).
#ifndef __RIGOR_H__
#define __RIGOR_H__
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <map>

namespace rig {

typedef std::complex<double> cpxf;
static const double PIf = 3.14159265358979323846;

// ---------------- sound real interval ----------------
static inline double dn(double x){ return std::nextafter(x, -INFINITY); }
static inline double up(double x){ return std::nextafter(x,  INFINITY); }

struct R {
  double lo, hi;
  R(): lo(0), hi(0) {}
  R(double x): lo(x), hi(x) {}
  R(double a, double b){ if(a<=b){lo=a;hi=b;} else {lo=b;hi=a;} }
};
static inline R radd(const R&a,const R&b){ return R(dn(a.lo+b.lo), up(a.hi+b.hi)); }
static inline R rsub(const R&a,const R&b){ return R(dn(a.lo-b.hi), up(a.hi-b.lo)); }
static inline R rneg(const R&a){ return R(-a.hi,-a.lo); }
static inline R rmul(const R&a,const R&b){
  double p1=a.lo*b.lo,p2=a.lo*b.hi,p3=a.hi*b.lo,p4=a.hi*b.hi;
  double mn=std::min(std::min(p1,p2),std::min(p3,p4));
  double mx=std::max(std::max(p1,p2),std::max(p3,p4));
  return R(dn(mn),up(mx));
}
static inline R rdiv(const R&a,const R&b){
  double p1=a.lo/b.lo,p2=a.lo/b.hi,p3=a.hi/b.lo,p4=a.hi/b.hi;
  double mn=std::min(std::min(p1,p2),std::min(p3,p4));
  double mx=std::max(std::max(p1,p2),std::max(p3,p4));
  return R(dn(mn),up(mx));
}
static inline R rsqrt(const R&a){
  double l = a.lo>0? std::sqrt(a.lo):0.0;
  double h = a.hi>0? std::sqrt(a.hi):0.0;
  return R(a.lo>0? dn(l):0.0, up(h));
}
static inline R rabsval(const R&x){ // |x| for a real interval
  if(x.lo<=0 && 0<=x.hi) return R(0.0, up(std::max(-x.lo,x.hi)));
  double m=std::min(std::fabs(x.lo),std::fabs(x.hi));
  double M=std::max(std::fabs(x.lo),std::fabs(x.hi));
  return R(m,M);
}

// ---------------- complex interval (rectangle) ----------------
struct C {
  R re, im;
  C(): re(0), im(0) {}
  C(const R&a,const R&b): re(a), im(b) {}
  C(double a,double b): re(a), im(b) {}
};
static inline C cadd(const C&a,const C&b){ return C(radd(a.re,b.re), radd(a.im,b.im)); }
static inline C csub(const C&a,const C&b){ return C(rsub(a.re,b.re), rsub(a.im,b.im)); }
static inline C cmul(const C&a,const C&b){
  return C(rsub(rmul(a.re,b.re),rmul(a.im,b.im)), radd(rmul(a.re,b.im),rmul(a.im,b.re)));
}
static inline R cabsval(const C&a){ // modulus interval
  R mre=rabsval(a.re), mim=rabsval(a.im);
  return rsqrt(radd(rmul(mre,mre),rmul(mim,mim)));
}
static inline R cdist(const C&a,const C&b){ return cabsval(csub(a,b)); }

// pi and 1/sqrt2 as intervals
static const R PI_I = R(dn(PIf), up(PIf));
static inline R inv_sqrt2(){ return rsqrt(R(0.5)); }

// interval cos/sin over theta (radians)
static inline R cos_iv(const R&th){
  std::vector<double> cand; cand.push_back(th.lo); cand.push_back(th.hi);
  long k0=(long)std::floor(th.lo/PIf)-1, k1=(long)std::ceil(th.hi/PIf)+1;
  for(long k=k0;k<=k1;++k){ double x=k*PIf; if(th.lo<=x && x<=th.hi) cand.push_back(x); }
  double mn=1e300,mx=-1e300;
  for(double x:cand){ double c=std::cos(x); mn=std::min(mn,c); mx=std::max(mx,c); }
  return R(dn(mn),up(mx));
}
static inline R sin_iv(const R&th){
  std::vector<double> cand; cand.push_back(th.lo); cand.push_back(th.hi);
  long k0=(long)std::floor((th.lo-PIf/2)/PIf)-1, k1=(long)std::ceil((th.hi-PIf/2)/PIf)+1;
  for(long k=k0;k<=k1;++k){ double x=PIf/2+k*PIf; if(th.lo<=x && x<=th.hi) cand.push_back(x); }
  double mn=1e300,mx=-1e300;
  for(double x:cand){ double s=std::sin(x); mn=std::min(mn,s); mx=std::max(mx,s); }
  return R(dn(mn),up(mx));
}
static inline C s_of_deg_box(double lo,double hi){
  R th = rmul(rmul(R(lo,hi), PI_I), R(1.0/180.0));
  R inv=inv_sqrt2();
  return C(rmul(inv,cos_iv(th)), rmul(inv,sin_iv(th)));
}

// ---------------- word-disk center (validated formula) ----------------
// center(u) = c0 + (s-1)/2 * sum_{i=1}^{n-1} eps_i s^i,  c0=s/2 (u0=f) or 1-s/2,
// eps_i=+1 if u[i]=='0' else -1.
static inline C word_center_iv(const std::string& u, const C& s){
  C half(0.5,0.0), one(1.0,0.0);
  C c0 = (u[0]=='0') ? cmul(half,s) : csub(one, cmul(half,s));
  int n=(int)u.size();
  if(n==1) return c0;
  C coeff = cmul(csub(s,one), half);
  C acc(0.0,0.0), spow=s;
  for(int i=1;i<n;++i){
    acc = cadd(acc, (u[i]=='0')? spow : C(rneg(spow.re),rneg(spow.im)));
    spow = cmul(spow,s);
  }
  return cadd(c0, cmul(coeff,acc));
}
static inline cpxf word_center_f(const std::string& u, cpxf s){
  cpxf c0 = (u[0]=='0') ? 0.5*s : (1.0 - 0.5*s);
  int n=(int)u.size(); if(n==1) return c0;
  cpxf coeff=(s-1.0)*0.5, acc(0,0), spow=s;
  for(int i=1;i<n;++i){ acc += (u[i]=='0')? spow : -spow; spow*=s; }
  return c0 + coeff*acc;
}
static inline cpxf s_of_deg_f(double deg){ double t=deg*PIf/180.0; return cpxf(std::cos(t),std::sin(t))/std::sqrt(2.0); }
static inline double mir_f(cpxf s){ return std::abs(s-1.0)*0.5/(1.0-std::abs(s)); }
static inline R mir_iv(const C& s){
  C one(1.0,0.0);
  R num=cabsval(csub(s,one));
  R den=rsub(R(1.0), cabsval(s));
  return rdiv(rmul(num,R(0.5)), den);
}

// ---------------- float disk-arrangement boundary ----------------
struct Disk { cpxf c; double r; int tag; };  // tag 0=u,1=v
struct Arc { int i; int tag; double a0,a1; cpxf P0,P1; };

static inline double normang(double a){ const double TP=2*PIf; a=std::fmod(a,TP); if(a<0)a+=TP; return a; }

static inline std::vector<cpxf> circle_ix(cpxf c1,double r1,cpxf c2,double r2){
  std::vector<cpxf> out; double d=std::abs(c2-c1);
  if(d>r1+r2 || d<std::fabs(r1-r2) || d<1e-18) return out;
  double a=(r1*r1-r2*r2+d*d)/(2*d); double h2=r1*r1-a*a; if(h2<0) return out;
  double h=std::sqrt(h2); cpxf mid=c1+a*(c2-c1)/d; cpxf perp=cpxf(0,1)*(c2-c1)/d;
  out.push_back(mid+h*perp); out.push_back(mid-h*perp); return out;
}

// build clouds (uniform radius) for word extended by all nd-bit suffixes
static inline void cloud(const std::string& w,int nd,cpxf s,double mir,int tag,std::vector<Disk>& out){
  int n=(int)w.size()+nd; double r=mir*std::pow(std::abs(s),n);
  int N=1<<nd;
  for(int k=0;k<N;++k){
    std::string suff; suff.reserve(nd);
    for(int b=nd-1;b>=0;--b) suff.push_back(((k>>b)&1)?'1':'0');
    Disk D; D.c=word_center_f(w+suff,s); D.r=r; D.tag=tag; out.push_back(D);
  }
}

struct Grid {
  std::vector<Disk>* D; double bs; std::map<long,std::vector<int> > buck;
  long key(cpxf c) const { long ix=(long)std::floor(c.real()/bs); long iy=(long)std::floor(c.imag()/bs); return ix*1000003L+iy; }
  void build(std::vector<Disk>& d){ D=&d; double mr=0; for(auto&x:d) mr=std::max(mr,x.r); bs=2.2*mr; if(bs<=0)bs=1;
    for(int i=0;i<(int)d.size();++i) buck[key(d[i].c)].push_back(i); }
  void neigh(cpxf c, std::vector<int>& out) const { out.clear(); long ix=(long)std::floor(c.real()/bs), iy=(long)std::floor(c.imag()/bs);
    for(int dx=-1;dx<=1;++dx) for(int dy=-1;dy<=1;++dy){ auto it=buck.find((ix+dx)*1000003L+(iy+dy)); if(it!=buck.end()) for(int j:it->second) out.push_back(j); } }
};

static inline bool inside_any(cpxf P, const std::vector<Disk>& D, const Grid& G, int ex1, int ex2, double eps){
  std::vector<int> nb; G.neigh(P,nb);
  for(int k:nb){ if(k==ex1||k==ex2) continue; if(std::abs(P-D[k].c) < D[k].r-eps) return true; }
  return false;
}

// boundary_arcs3: midpoint-test.  Returns arcs (all boundary arcs).
static inline void boundary_arcs3(std::vector<Disk>& D, const Grid& G, std::vector<Arc>& arcs, double eps=1e-11){
  arcs.clear();
  std::vector<int> nb;
  for(int i=0;i<(int)D.size();++i){
    cpxf ci=D[i].c; double ri=D[i].r;
    // boundary vertices on circle i
    std::vector<std::pair<double,cpxf> > Vs;
    G.neigh(ci,nb);
    for(int j:nb){ if(j==i) continue;
      std::vector<cpxf> pts=circle_ix(ci,ri,D[j].c,D[j].r);
      for(cpxf P:pts){ if(!inside_any(P,D,G,i,j,eps)){ double ang=normang(std::atan2((P-ci).imag(),(P-ci).real())); Vs.push_back(std::make_pair(ang,P)); } }
    }
    if(Vs.empty()) continue;
    std::sort(Vs.begin(),Vs.end(),[](const std::pair<double,cpxf>&a,const std::pair<double,cpxf>&b){return a.first<b.first;});
    int m=(int)Vs.size();
    for(int k=0;k<m;++k){
      double a0=Vs[k].first; cpxf P0=Vs[k].second;
      double a1=Vs[(k+1)%m].first; cpxf P1=Vs[(k+1)%m].second;
      double span=normang(a1-a0); if(span<1e-12) continue;
      double am=normang(a0+span*0.5); cpxf M=ci+ri*cpxf(std::cos(am),std::sin(am));
      if(!inside_any(M,D,G,i,-1,eps)){ Arc A; A.i=i; A.tag=D[i].tag; A.a0=a0; A.a1=normang(a0)+span; A.P0=P0; A.P1=P1; arcs.push_back(A); }
    }
  }
}

// trace the outer loop; return arc indices in cyclic order
static inline std::vector<int> trace_outer_idx(const std::vector<Arc>& arcs){
  std::vector<int> res; if(arcs.empty()) return res;
  auto vkey=[](cpxf P){ long a=(long)std::llround(P.real()*1e8); long b=(long)std::llround(P.imag()*1e8); return a*100000007L+b; };
  std::map<long,std::vector<int> > vadj;
  for(int idx=0;idx<(int)arcs.size();++idx){ vadj[vkey(arcs[idx].P0)].push_back(idx); vadj[vkey(arcs[idx].P1)].push_back(idx); }
  // rightmost vertex
  long bestv=0; double bestx=-1e300; bool has=false;
  // need actual coords for rightmost: track separately
  std::map<long,cpxf> vpt;
  for(int idx=0;idx<(int)arcs.size();++idx){ vpt[vkey(arcs[idx].P0)]=arcs[idx].P0; vpt[vkey(arcs[idx].P1)]=arcs[idx].P1; }
  for(auto&kv:vpt){ if(kv.second.real()>bestx){ bestx=kv.second.real(); bestv=kv.first; has=true; } }
  if(!has) return res;
  if(vadj[bestv].empty()) return res;
  int start=vadj[bestv][0]; int cur=start; long curv=bestv; int steps=0, maxs=(int)arcs.size()*2+5;
  while(steps<maxs){
    res.push_back(cur);
    long v0=vkey(arcs[cur].P0), v1=vkey(arcs[cur].P1);
    long nxtv = (v0==curv)? v1 : v0;
    int nxt=-1; for(int ai:vadj[nxtv]){ if(ai!=cur){ nxt=ai; break; } }
    if(nxt<0) break;
    cur=nxt; curv=nxtv; ++steps;
    if(cur==start) break;
  }
  return res;
}

static inline void alt_stats(const std::vector<int>& seq_tags, int& runs,int& ur,int& vr){
  runs=0;ur=0;vr=0; int m=(int)seq_tags.size(); if(m==0) return;
  int s=0; while(s<m && seq_tags[s]==seq_tags[0]) ++s;
  if(s==m){ runs=1; if(seq_tags[0]==0)ur=m; else vr=m; return; }
  std::vector<int> q; for(int k=0;k<m;++k) q.push_back(seq_tags[(s+k)%m]);
  runs=1; int prev=-1;
  for(int k=0;k<m;++k){ if(k>0 && q[k]!=q[k-1]) runs++; if(q[k]!=prev){ if(q[k]==0)ur++; else vr++; prev=q[k]; } }
}

// ---------------- the interval certificate over an s-box ----------------
// returns 1 if certified, 0 if not-a-trap (float), -1 if float trap but interval
// verification failed.  min_margin set on success.
static inline int certify_box(const std::string& u,const std::string& v,int nd,
                              double deg_lo,double deg_hi,double& min_margin){
  min_margin=-1.0;
  double degm=0.5*(deg_lo+deg_hi);
  cpxf sm=s_of_deg_f(degm); double mirm=mir_f(sm);
  std::vector<Disk> D; cloud(u,nd,sm,mirm,0,D); cloud(v,nd,sm,mirm,1,D);
  Grid G; G.build(D);
  std::vector<Arc> arcs; boundary_arcs3(D,G,arcs);
  std::vector<int> outer=trace_outer_idx(arcs);
  std::vector<int> tags; for(int ai:outer) tags.push_back(arcs[ai].tag);
  int runs,ur,vr; alt_stats(tags,runs,ur,vr);
  if(!(ur>=2 && vr>=2)) return 0;

  // interval disks over the box
  C S=s_of_deg_box(deg_lo,deg_hi); R mirS=mir_iv(S);
  int n_u=(int)u.size()+nd, n_v=(int)v.size()+nd;
  R aS=cabsval(S);
  R ruI=mirS, rvI=mirS;
  for(int k=0;k<n_u;++k) ruI=rmul(ruI,aS);
  for(int k=0;k<n_v;++k) rvI=rmul(rvI,aS);
  // interval centers, only need: for each outer u-arc's disk -> its center; and
  // ALL opposite-cloud centers/radii.  Compute interval centers for all disks
  // (same order as D).
  std::vector<C> Civ(D.size());
  {
    int idx=0; int Nu=1<<nd;
    for(int k=0;k<Nu;++k){ std::string suff; for(int b=nd-1;b>=0;--b) suff.push_back(((k>>b)&1)?'1':'0'); Civ[idx++]=word_center_iv(u+suff,S); }
    for(int k=0;k<Nu;++k){ std::string suff; for(int b=nd-1;b>=0;--b) suff.push_back(((k>>b)&1)?'1':'0'); Civ[idx++]=word_center_iv(v+suff,S); }
  }
  double mm=1e300;
  std::vector<int> nb;
  for(int ai:outer){
    const Arc& A=arcs[ai]; int i=A.i; int ti=A.tag;
    double span=normang(A.a1-A.a0); double am=normang(A.a0+span*0.5);
    // interval midpoint on circle i: M = Civ[i] + r_i * (cos am, sin am)
    R ri = (ti==0)? ruI : rvI;
    C M = cadd(Civ[i], C(rmul(ri,R(std::cos(am))), rmul(ri,R(std::sin(am)))));
    // certify M certainly outside every OPPOSITE-cloud disk (Z/W identity)
    G.neigh(D[i].c,nb);
    for(int k:nb){ if(k==i) continue; if(D[k].tag==ti) continue; // opposite only
      R rk=(D[k].tag==0)? ruI : rvI;
      R dist=cdist(M,Civ[k]);
      double marg = dist.lo - rk.hi;   // certainly outside iff >0
      if(marg<=0) { return -1; }
      if(marg<mm) mm=marg;
    }
  }
  min_margin=mm; return 1;
}

// ---------------- SEGMENT-PERSISTENCE certificate over an s-box ----------------
// Sound realization of the relaxed (CKW 4-arc) criterion that survives balls popping
// in/out.  Fix a center c in the overlap uD∩vD.  For a ray from c in direction θ, let
// E_u(θ)=max over u-disks (hit by the ray) of the far-intersection distance, likewise
// E_v(θ).  The outermost point of the union in direction θ is u-tagged iff E_u>E_v, and
// (since it is beyond the opposite cloud's extent on that ray) it lies OUTSIDE the
// opposite cloud.  If at 4 directions θ1<θ2<θ3<θ4 the tag alternates u,v,u,v with a
// STRICT extent gap over the whole box, the four outermost points are 4 interleaved
// points on the outer Jordan boundary (c∈W interior) -> the u-arc and v-arc link ->
// trap, for ALL s in the box.  E is a MAX over disks, so which specific ball realizes a
// segment may change with s; only a segment being fully overtaken by the opposite cloud
// (extent ordering flipping) can break it, and that is exactly what the gap test catches.
// returns 1 certified / 0 not-a-trap(float) / -1 interval-failed.  min_margin=min gap.
static inline double ray_extent_f(cpxf c,double ct,double st,cpxf dc,double r){
  cpxf w=dc-c; double p=w.real()*ct+w.imag()*st;
  double perp=-w.real()*st+w.imag()*ct; double perp2=perp*perp; // no |w|^2-p^2 cancellation
  if(perp2> r*r) return -1e300; double h=std::sqrt(r*r-perp2); double e=p+h; return e>0?e:-1e300;
}
static inline int certify_box_seg(const std::string& u,const std::string& v,int nd,
                                  double deg_lo,double deg_hi,double& min_margin){
  min_margin=-1.0;
  double degm=0.5*(deg_lo+deg_hi);
  cpxf sm=s_of_deg_f(degm); double mirm=mir_f(sm);
  std::vector<Disk> D; cloud(u,nd,sm,mirm,0,D); int nU=(int)D.size(); cloud(v,nd,sm,mirm,1,D);
  double ruf=mirm*std::pow(std::abs(sm),(int)u.size()+nd);
  double rvf=mirm*std::pow(std::abs(sm),(int)v.size()+nd);
  // center c = midpoint of the u-disk/v-disk pair with the DEEPEST overlap.  c is defined
  // as a FUNCTION OF s (midpoint of two word-centers) so that dist(c, disk-center) is a
  // difference of correlated word-centers -- far less interval blow-up than a fixed c.
  cpxf c(0,0); double bestdepth=-1e300;
  for(int j=0;j<nU;++j) for(int k=nU;k<(int)D.size();++k){
    double dd=std::abs(D[j].c-D[k].c); double depth=std::min(ruf,rvf)-0.5*dd;
    if(depth>bestdepth){ bestdepth=depth; c=0.5*(D[j].c+D[k].c); } }
  // float extents around c
  const int K=720;
  std::vector<double> Eu(K),Ev(K),gap(K); std::vector<int> tag(K);
  for(int k=0;k<K;++k){ double th=2*PIf*k/K, ct=std::cos(th), st=std::sin(th);
    double eu=-1e300,ev=-1e300;
    for(int j=0;j<nU;++j){ double e=ray_extent_f(c,ct,st,D[j].c,ruf); if(e>eu)eu=e; }
    for(int j=nU;j<(int)D.size();++j){ double e=ray_extent_f(c,ct,st,D[j].c,rvf); if(e>ev)ev=e; }
    Eu[k]=eu;Ev[k]=ev; gap[k]=std::fabs(eu-ev); tag[k]=(eu>ev)?0:1;
  }
  // runs of constant tag (cyclic); representative = argmax gap in the run
  std::vector<int> rtag; std::vector<double> rth, rgap;
  { int k0=0; while(k0<K && tag[k0]==tag[0]) ++k0; // start at a boundary (or all-same)
    if(k0==K) return 0; // single tag -> not a trap
    int start=k0; int k=start;
    for(int cnt=0;cnt<K;){ int t=tag[k%K]; double bg=-1; int bk=k;
      while(cnt<K && tag[k%K]==t){ if(gap[k%K]>bg){bg=gap[k%K];bk=k;} ++k;++cnt; }
      rtag.push_back(t); rth.push_back(2*PIf*(bk%K)/K); rgap.push_back(bg);
    }
  }
  int nr=(int)rtag.size(); int nur=0,nvr=0; for(int t:rtag){ if(t==0)nur++; else nvr++; }
  if(!(nur>=2 && nvr>=2)) return 0;
  // choose 4 consecutive runs (they alternate u,v,u,v) maximizing the min gap
  int besti=-1; double bestmin=-1;
  for(int i=0;i<nr;++i){
    double mn=1e300; for(int d=0;d<4;++d) mn=std::min(mn,rgap[(i+d)%nr]);
    // must alternate over the 4
    bool alt=true; for(int d=0;d<4;++d) if(rtag[(i+d)%nr]==rtag[(i+d+1)%nr]) alt=false;
    if(alt && mn>bestmin){ bestmin=mn; besti=i; }
  }
  if(besti<0) return 0;
  double th[4]; int tg[4];
  for(int d=0;d<4;++d){ th[d]=rth[(besti+d)%nr]; tg[d]=rtag[(besti+d)%nr]; }

  // ---- interval verification over the box ----
  C S=s_of_deg_box(deg_lo,deg_hi);
  R mirS=mir_iv(S); R aS=cabsval(S);
  R ruI=mirS,rvI=mirS; for(int k=0;k<(int)u.size()+nd;++k) ruI=rmul(ruI,aS);
  for(int k=0;k<(int)v.size()+nd;++k) rvI=rmul(rvI,aS);
  std::vector<C> Civ(D.size());
  { int idx=0,Nu=1<<nd;
    for(int k=0;k<Nu;++k){ std::string s2; for(int b=nd-1;b>=0;--b) s2.push_back(((k>>b)&1)?'1':'0'); Civ[idx++]=word_center_iv(u+s2,S); }
    for(int k=0;k<Nu;++k){ std::string s2; for(int b=nd-1;b>=0;--b) s2.push_back(((k>>b)&1)?'1':'0'); Civ[idx++]=word_center_iv(v+s2,S); }
  }
  C cC(c.real(),c.imag());  // fixed center (from float midpoint)
  // (1) c certainly in uD∩vD over the box: some u-disk and some v-disk certainly contain c.
  bool cin_u=false,cin_v=false;
  for(int k=0;k<nU;++k){ R d=cdist(cC,Civ[k]); if(d.hi<ruI.lo){ cin_u=true; break; } }
  for(int k=nU;k<(int)D.size();++k){ R d=cdist(cC,Civ[k]); if(d.hi<rvI.lo){ cin_v=true; break; } }
  if(std::getenv("RIGOR_DIAG")) std::fprintf(stderr,"[seg] cin_u=%d cin_v=%d nruns=%d bestmin_gap=%.3e\n",(int)cin_u,(int)cin_v,nr,bestmin);
  if(!(cin_u&&cin_v)) return -1;
  // (2) for each of the 4 directions, verify the strict extent ordering over the box.
  //   E_u.lo = max over u-disks that CERTAINLY hit of ext.lo ; E_v.hi = max over v-disks
  //   that MIGHT hit of ext.hi ; u-dir needs E_u.lo > E_v.hi ; v-dir the reverse.
  double gmin=1e300;
  for(int d=0;d<4;++d){
    double ct=std::cos(th[d]), st=std::sin(th[d]);
    R Eu_lo(-1e300), Ev_lo(-1e300), Eu_hi(-1e300), Ev_hi(-1e300); // track needed bounds
    // helper via lambda-like inline
    for(int side=0; side<2; ++side){
      int lo=(side==0)?0:nU, hi=(side==0)?nU:(int)D.size(); R rr=(side==0)?ruI:rvI;
      double bestcert=-1e300, bestmaybe=-1e300;
      for(int j=lo;j<hi;++j){
        C w=csub(Civ[j],cC);
        R p=radd(rmul(w.re,R(ct)),rmul(w.im,R(st)));
        R perp=radd(rmul(w.re,R(-st)),rmul(w.im,R(ct))); // perpendicular component
        R perp2=rmul(perp,perp);                          // no |w|^2-p^2 cancellation
        R r2=rmul(rr,rr);
        R diff=rsub(r2,perp2);
        if(diff.lo>0){ R h=rsqrt(diff); R ext=radd(p,h);
          if(ext.lo>bestcert) bestcert=ext.lo;
          if(ext.hi>bestmaybe) bestmaybe=ext.hi; }
        else if(diff.hi>0){ R h=rsqrt(R(0.0,diff.hi)); R ext=radd(p,h);
          if(ext.hi>bestmaybe) bestmaybe=ext.hi; }
      }
      if(side==0){ Eu_lo=R(bestcert); Eu_hi=R(bestmaybe); }
      else       { Ev_lo=R(bestcert); Ev_hi=R(bestmaybe); }
    }
    double g;
    if(tg[d]==0) g = Eu_lo.lo - Ev_hi.lo;   // u outermost: E_u.lo - E_v.hi
    else         g = Ev_lo.lo - Eu_hi.lo;   // v outermost: E_v.lo - E_u.hi
    if(g<=0) return -1;
    if(g<gmin) gmin=g;
  }
  min_margin=gmin; return 1;
}

// ---------------- THIN-TRAP interval certificate (Phase 3) ----------------
// A "thin trap" certifies a coincidence-CORE point (45 deg, ~69.3 deg) where the CKW and
// segment traps both fail.  The finite,
// interval-verifiable part is: two chords -- an f-chord [A,B] with A=wA.0, B=wB.0 (words
// starting f) and a g-chord [C,D] with C=wC.0, D=wD.0 (words starting g) -- that CROSS
// TRANSVERSELY, verified over an s-box.  If the transverse crossing holds over the box
// (an open condition), then by the contraction/self-reproduction lemma
// the substitution recursion converges to an exact infinite coincidence u'.0=v'.0 in
// fLambda ∩ gLambda for every s in the box => the box ⊂ interior(M).  THIS function
// verifies the finite crossing part rigorously; the contraction lemma is the (hand-proved)
// bridge to interior(M).
//
// seg_cross_iv: do segments A->B and C->D CERTAINLY cross (params t,u certainly in (0,1))
// with acute crossing angle CERTAINLY >= theta0, over the whole s-box?  Returns 1/-1 and
// sets margins (t,u interior margins; trans = |r x sv| - sin(theta0)|r||sv|).
static inline R rcross(const C&r,const C&sv){ // 2D cross product r x sv (real interval)
  return rsub(rmul(r.re,sv.im), rmul(r.im,sv.re));
}
static inline int seg_cross_iv(const C&A,const C&B,const C&Cc,const C&Dd,double sin_theta0,
                               double&t_marg,double&u_marg,double&trans_marg){
  C r=csub(B,A), sv=csub(Dd,Cc), qp=csub(Cc,A);
  R rxs=rcross(r,sv);
  // (1) transversality: |rxs| >= sin(theta0)*|r|*|sv| (also => rxs certainly nonzero)
  R absrxs=rabsval(rxs);
  R lensprod=rmul(cabsval(r),cabsval(sv));
  trans_marg = absrxs.lo - sin_theta0*lensprod.hi;
  t_marg=-1; u_marg=-1;
  bool rxs_nz = (rxs.lo>0.0 || rxs.hi<0.0);
  R t,u;
  if(rxs_nz){ t=rdiv(rcross(qp,sv),rxs); u=rdiv(rcross(qp,r),rxs);
    t_marg=std::min(t.lo,1.0-t.hi); u_marg=std::min(u.lo,1.0-u.hi); }
  if(std::getenv("RIGOR_DIAG"))
    std::fprintf(stderr,"[thin] rxs=[%.3e,%.3e] nz=%d trans=%.3e t=[%.4f,%.4f] u=[%.4f,%.4f]\n",
      rxs.lo,rxs.hi,(int)rxs_nz,trans_marg,t.lo,t.hi,u.lo,u.hi);
  if(!rxs_nz) return -1;                 // divisor might be 0
  if(trans_marg<=0.0) return -1;         // not certainly transverse
  if(t_marg<=0.0 || u_marg<=0.0) return -1;
  return 1;
}
// certify one thin-trap crossing (4 chord-endpoint words) over the arg-box [deg_lo,deg_hi].
static inline int certify_thin_trap(const std::string&wA,const std::string&wB,
                                    const std::string&wC,const std::string&wD,
                                    double deg_lo,double deg_hi,double theta0_deg,
                                    double&min_margin){
  C S=s_of_deg_box(deg_lo,deg_hi);
  C A=word_center_iv(wA,S), B=word_center_iv(wB,S);
  C Cc=word_center_iv(wC,S), Dd=word_center_iv(wD,S);
  double st=std::sin(theta0_deg*PIf/180.0);
  double tm,um,trm; int rc=seg_cross_iv(A,B,Cc,Dd,st,tm,um,trm);
  min_margin = (rc==1)? std::min(std::min(tm,um),trm) : -1.0;
  return rc;
}

} // namespace rig
#endif
