#include"checks.h"
#include"../src/rings.h"
#include<stdio.h>

static const double eps = 1e-4;

static double
eg(const double r[3], const double v[3]) {
  return 0.5*dot(v,v) - 1.0/norm(r);
}

static int
mean_motion_test() {
  body b;
  double mm;

  b.m = 1e-3;
  b.a = 0.237; 
  mm = mean_motion(&b);

  return check_close(eps, eps, mm, 8.67151);
}

static int
magnitude_check(const body *b) {
  double e = norm(b->A);
  double sqrt1me2 = norm(b->L);
  
  return check_close(eps,eps,sqrt(1.0-e*e), sqrt1me2);
}

static int
energy_test() {
  const int n = 7;
  const double Es[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double egs[n];
  const body b = {1e-6, 0.237, 0.0, 0.0, 0.0,
                  {0.048121633812668, 0.346123635235211, 0.239632435658364},
                  {0.807135311203869, 0.152173745618516, -0.381883193952636},
                  {0.0, 0.0, 0.0}};
  double r[3], v[3];
  int i;

  if (!check_close(eps,eps,dot(b.L,b.A), 0.0)) return 0;
  if (!magnitude_check(&b)) return 0;

  for(i = 0; i < n; i++) {
    E_to_rv(&b, Es[i], r, v);
    egs[i] = eg(r,v);
  }

  for(i = 1; i < n; i++) {
    if (!check_close(eps,eps,egs[i],egs[i-1])) return 0;
  }

  if (!check_close(eps,eps,egs[0],-1.0/2.0/b.a)) return 0;

  return 1;
}

static int
zero_ecc_energy_test() {
  const int n = 7;
  const double Es[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double egs[n];
  const body b = {1e-6, 0.237, 0.0, 0.0, 0.0,
                  {0.113569, 0.816863, 0.565541},
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0}};
  double r[3], v[3];
  int i;

  if (!check_close(eps,eps,dot(b.L,b.A), 0.0)) return 0;
  if (!magnitude_check(&b)) return 0;

  for(i = 0; i < n; i++) {
    E_to_rv(&b, Es[i], r, v);
    egs[i] = eg(r,v);
  }

  for(i = 1; i < n; i++) {
    if (!check_close(eps,eps,egs[i],egs[i-1])) return 0;
  }

  if (!check_close(eps,eps,egs[0],-1.0/2.0/b.a)) return 0;

  return 1;
}

int
main() {

  if (!mean_motion_test()) return 2;

  if (!energy_test()) return 1;

  if (!zero_ecc_energy_test()) return 3;
  
  return 0;
}
