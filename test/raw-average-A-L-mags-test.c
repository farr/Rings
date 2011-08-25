#include"../src/rings.h"
#include"checks.h"

#include<stdio.h>
#include<stdlib.h>

#define EPS 1e-10

int
main() {
  body b1, b2;
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  int status = 0;
  gsl_rng *rng;
  double rhs[BODY_VECTOR_SIZE];

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != 0);
  assert(ws2 != 0);

  init_random_body(rng, &b1, random_between(rng, 1e-3, 2e-3), random_between(rng, 1.0, 1.5), 0.0, 0.0, 0.0);
  init_random_body(rng, &b2, random_between(rng, 1e-3, 2e-3), random_between(rng, 10.0, 15.0), 0.0, 0.0, 0.0);

  raw_average_rhs(0.0, &b1, &b2, ws1, ws_size, ws2, ws_size, EPS, EPS, rhs);

  if (fabs(dot(b1.A, rhs+BODY_A_INDEX) + dot(b1.L, rhs+BODY_L_INDEX)) > 1e-8) {
    fprintf(stderr, "raw-average-A-L-mags-test: non-zero change in A^2 + L^2\n");
    status = 1;
  }

  gsl_rng_free(rng);
  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
