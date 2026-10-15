// review/exact/clipdrv.cpp — driver minimo: legge interi (come li riceve Clipper da polyclip) e calcola A ∩ B, fill nonzero,
// come interface.cpp di polyclip (Clipper c; AddPaths(A, ptSubject, true); AddPaths(B, ptClip, true); Execute(ctIntersection, ...)).
// Area esatta: 2*area in __int128 (shoelace intero). Formato input: "C <case> <variant> <nringsA> <nringsB>" poi per anello "R <n>" e n righe "x y".
#include <cstdio>
#include <cstring>
#include <string>
#include CLIPPER_HEADER
using namespace ClipperLib;
static std::string i128(__int128 v){ if(v==0) return "0"; bool neg=v<0; if(neg) v=-v; std::string s; while(v){ s.insert(s.begin(), char('0'+int(v%10))); v/=10;} return (neg?"-":"")+s; }
int main(){
  char tag[8]; int cs; char var[8]; int na, nb;
  while (scanf("%7s %d %7s %d %d", tag, &cs, var, &na, &nb) == 5) {
    Paths A(na), B(nb);
    for (int q = 0; q < na + nb; q++) { char r[4]; int n; scanf("%3s %d", r, &n); Path &p = (q < na) ? A[q] : B[q - na];
      for (int k = 0; k < n; k++) { long long x, y; scanf("%lld %lld", &x, &y); p.push_back(IntPoint(x, y)); } }
    Clipper c; c.AddPaths(A, ptSubject, true); c.AddPaths(B, ptClip, true); Paths S; c.Execute(ctIntersection, S, pftNonZero, pftNonZero);
    __int128 a2 = 0; int nv = 0;
    for (auto &p : S) { int n = p.size(); nv += n; for (int k = 0; k < n; k++) { const IntPoint &u = p[k], &v = p[(k + 1) % n]; a2 += (__int128)u.X * v.Y - (__int128)v.X * u.Y; } }
    printf("%d %s %s %d %d\n", cs, var, i128(a2).c_str(), (int)S.size(), nv);
  }
  return 0;
}
