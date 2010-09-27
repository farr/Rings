#include"../src/rings.h"
#include"checks.h"
#include<stdlib.h>
#include<stdio.h>

static double
random_between(const double a, const double b) {
  return (b-a)*(a + rand()/(RAND_MAX + 1.0));
}

int
main() {
  double a, e, I, Omega, omega, m;
  double myE, myI, myOmega, myomega;
  body b;

  sranddev();

  a = random_between(1.0, 2.0);
  m = random_between(1.0, 2.0);
  e = random_between(0.0, 1.0);
  I = random_between(0.0, 180.0);
  Omega = random_between(0.0, 360.0);
  omega = random_between(0.0, 360.0);

  init_body_from_elements(&b, m, a, e, I, Omega, omega);
  elements_from_body(&b, &myE, &myI, &myOmega, &myomega);

  if (b.m != m) {
    return 1;
  }

  if (b.a != a) {
    return 2;
  }

  if (!check_close(1e-8, 1e-8, e, myE)) {
    return 3;
  }

  if (!check_close(1e-8, 1e-8, I, myI)) {
    return 4;
  }
  
  if (!check_close(1e-8, 1e-8, Omega, myOmega)) {
    fprintf(stderr, "alloc-body-from-elements: bad omega (got %g, expected %g)\n", myOmega, Omega);
    return 5;
  }

  if (!check_close(1e-8, 1e-8, omega, myomega)) {
    return 6;
  }

  return 0;
}
