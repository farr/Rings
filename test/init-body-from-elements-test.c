#include"../src/rings.h"
#include"checks.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

int
main() {
  double a, e, I, Omega, omega, m;
  double myE, myI, myOmega, myomega;
  body b;
  double Qp = 0.0, inertia = 0.0, spin[3] = {0.0, 0.0, 0.0};
  gsl_rng *rng;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);

  assert(rng != 0);
  seed_random(rng);

  a = random_between(rng, 1.0, 2.0);
  m = random_between(rng, 1.0, 2.0);
  e = random_between(rng, 0.0, 1.0);
  I = random_between(rng, 0.0, 180.0);
  Omega = random_between(rng, 0.0, 360.0);
  omega = random_between(rng, 0.0, 360.0);

  init_body_from_elements(&b, m, a, e, I, Omega, omega, spin, Qp, inertia);
  elements_from_body(&b, &myE, &myI, &myOmega, &myomega);

  if (b.m != m) {
    gsl_rng_free(rng);
    return 1;
  }

  if (b.a != a) {
    gsl_rng_free(rng);
    return 2;
  }

  if (!check_close(1e-8, 1e-8, e, myE)) {
    gsl_rng_free(rng);
    return 3;
  }

  if (!check_close(1e-8, 1e-8, I, myI)) {
    gsl_rng_free(rng);
    return 4;
  }
  
  if (!check_close(1e-8, 1e-8, Omega, myOmega)) {
    fprintf(stderr, "alloc-body-from-elements: bad omega (got %g, expected %g)\n", myOmega, Omega);
    gsl_rng_free(rng);
    return 5;
  }

  if (!check_close(1e-8, 1e-8, omega, myomega)) {
      gsl_rng_free(rng);
      return 6;
  }

  gsl_rng_free(rng);
  return 0;
}
