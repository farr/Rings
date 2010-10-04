#include"../src/rings.h"
#include"checks.h"

#include<stdlib.h>

#define EPS 1e-10

int
main() {
  body b1, b2;
  gsl_rng *rng;
  const double m1 = 1e-3;
  const double m2 = 2e-3;
  const double a1 = 1.2;
  const double a2 = 2.3;
  double rhs[BODY_VECTOR_SIZE];
  double eps;
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  int status = 0;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);

  assert(ws1 != 0);
  assert(ws2 != 0);

  init_random_body(rng, &b1, m1, a1);
  init_random_body(rng, &b2, m2, a2);
  eps = random_between(rng, 1e-4, 1e-3);

  raw_average_rhs(eps, &b1, &b2, ws1, ws_size, ws2, ws_size, EPS, EPS, rhs);

  if (rhs[0] != 0.0) {
    fprintf(stderr, "raw-average-zero-adot-test: rhs[0] != 0.0\n");
    status = 1;
  }

  if (fabs(rhs[1]) > 1e-8) {
    fprintf(stderr, "raw-average-zero-adot-test: |rhs[1]| > 1e-8\n");
    status = 2;
  }

  gsl_rng_free(rng);
  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
