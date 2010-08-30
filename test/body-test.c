#include"checks.h"
#include"../src/rings.h"
#include<stdio.h>

static const double eps = 1e-4;

static double
eg(const double r[3], const double v[3]) {
  return 0.5*dot(v,v) - 1.0/norm(r);
}

static int
energy_test() {
  const int n = 7;
  const double Es[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double egs[n];
  const body b = {1e-6, 0.237, 
                  {0.683811, 0.68616, 0.301055},
                  {-0.22306316859388583158, 0.19747183082696724358, 0.05658692577476481217}};
  double r[3], v[3];
  int i;

  if (!check_close(eps,eps,dot(b.L,b.A), 0.0)) return 0;

  for(i = 0; i < n; i++) {
    EtoRv(&b, Es[i], r, v);
    egs[i] = eg(r,v);
  }

  for(i = 1; i < n; i++) {
    fprintf(stderr, "E = %g, energy = %g (last energy = %g)\n", Es[i], egs[i], egs[i-1]);
    if (!check_close(eps,eps,egs[i],egs[i-1])) return 0;
  }

  if (!check_close(eps,eps,egs[0],1.0/2.0/b.a)) return 0;

  return 1;
}

int
main() {

  if (!energy_test()) return 1;
  
  return 0;
}
