#include"checks.h"
#include"../src/rings.h"

#include<gsl/gsl_rng.h>

int
main() {
  body b1, b2;
  double rhs1[BODY_VECTOR_SIZE], rhs2[BODY_VECTOR_SIZE];
  gsl_rng *rng;
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  const double m1 = 1.2e-3;
  const double m2 = 0.7e-3;
  const double a1 = 1.15e-3;
  const double a2 = 1.98e-3;
  double dLtot[3];
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

  raw_average_rhs(0.0, &b1, &b2, ws1, ws_size, ws2, ws_size, rhs1);
  raw_average_rhs(0.0, &b2, &b1, ws1, ws_size, ws2, ws_size, rhs2);

  vscale(b1.m*mean_motion(&b1)*a1*a1, rhs1+BODY_L_INDEX, rhs1+BODY_L_INDEX);
  vscale(b2.m*mean_motion(&b2)*a2*a2, rhs2+BODY_L_INDEX, rhs2+BODY_L_INDEX);

  vadd(rhs1+BODY_L_INDEX, rhs2+BODY_L_INDEX, dLtot);

  if (norm(dLtot) > 1e-8) {
    fprintf(stderr, "raw-average-Ltot-test: expected dLtot/dt = 0, got %g", norm(dLtot));
    status = 1;
  }

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  gsl_rng_free(rng);
  return status;
}
