#include"../src/rings.h"
#include"checks.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

int
main() {
  double x[3], rxx[3], rzx[3];
  double theta;
  double temp;
  int i;
  gsl_rng *rng;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  theta = 2.0*(rand() / (RAND_MAX + 1.0) - 0.5);
  for (i = 0; i < 3; i++) {
    x[i] = 2.0*(rand() / (RAND_MAX + 1.0) - 0.5);
  }

  rotate_x(x, theta, rxx);
  rotate_z(x, theta, rzx);

  if (!check_close(1e-8, 1e-8, x[0], rxx[0])) {
    gsl_rng_free(rng);
    return 1;
  }

  if (!check_close(1e-8, 1e-8, x[2], rzx[2])) {
    gsl_rng_free(rng);
    return 2;
  }

  if (!check_close(1e-8, 1e-8, norm(rxx), norm(x))) {
    gsl_rng_free(rng);
    return 3;
  }

  if (!check_close(1e-8, 1e-8, norm(rzx), norm(x))) {
    gsl_rng_free(rng);
    return 4;
  }

  temp = x[0];
  x[0] = 0.0;
  rxx[0] = 0.0;
  if (!check_close(1e-8, 1e-8, dot(x, rxx) / norm(x) / norm(rxx), cos(theta))) {
    gsl_rng_free(rng);
    return 5;
  }

  x[0] = temp;
  temp = x[2];
  x[2] = 0.0;
  rzx[2] = 0.0;
  if (!check_close(1e-8, 1e-8, dot(x, rzx) / norm(x) / norm(rzx), cos(theta))) {
    gsl_rng_free(rng);
    return 6;
  }

  gsl_rng_free(rng);
  return 0;
}
