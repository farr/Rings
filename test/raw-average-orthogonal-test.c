#include"checks.h"
#include"../src/rings.h"

#include<gsl/gsl_rng.h>

int
main() {
  body b1, b2;
  double rhs[BODY_VECTOR_SIZE];
  gsl_rng *rng;
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  int status = 0;
  const double m1 = 1.1e-3;
  const double m2 = 2.05e-3;
  const double a1 = 1.03;
  const double a2 = 1.99;
  
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);

  assert(ws1 != 0);
  assert(ws2 != 0);

  init_random_body(rng, &b1, m1, a1);
  init_random_body(rng, &b2, m2, a2);

  raw_average_rhs(0.0, &b1, &b2, ws1, ws_size, ws2, ws_size, rhs);

  if (norm(rhs+BODY_A_INDEX) < 1e-8) {
    fprintf(stderr, "raw-average-orthogonal-test: zero-norm Adot");
    status = 2;
  }

  if (norm(rhs+BODY_L_INDEX) < 1e-8) {
    fprintf(stderr, "raw-average-orthogonal-test: zero norm Ldot");
    status = 3;
  }
  
  if (fabs(dot(b1.L, rhs+BODY_A_INDEX) + dot(b1.A, rhs+BODY_L_INDEX)) > 1e-8) {
    fprintf(stderr, "raw-average-orthogonal-test: expected orthogonal changes to A and L\n");
    status = 1;
  }

  gsl_rng_free(rng);
  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
