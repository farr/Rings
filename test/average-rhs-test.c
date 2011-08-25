#include"../src/rings.h"
#include"checks.h"

#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_rng.h>

#define EPS 1e-8

int
main() {
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  gsl_rng *rng;
  int status = 0;
  body b1, b2;
  const double m1 = 1.01e-3;
  const double m2 = 1.998e-3;
  const double a1 = 1.02;
  const double a2 = 10.3;
  const double eps = 0.001;
  double rhs[BODY_VECTOR_SIZE];
  double numerical_rhs[BODY_VECTOR_SIZE];
  int avg_status;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != 0);
  assert(ws2 != 0);

  init_random_body(rng, &b1, m1, a1, 0.0, 0.0, 0.0);
  init_random_body(rng, &b2, m2, a2, 0.0, 0.0, 0.0);

  raw_average_rhs(eps, &b1, &b2, ws1, ws_size, ws2, ws_size, EPS, EPS, numerical_rhs);
  avg_status = average_rhs(eps, &b1, &b2, EPS, rhs);

  if (avg_status != GSL_SUCCESS) {
    fprintf(stderr, "average-rhs-test: quad reported error %d.\n", avg_status);
    status = 2;
  }

  if (!check_vector_close(10.0*EPS, 10.0*EPS, BODY_VECTOR_SIZE, rhs, numerical_rhs)) {
    fprintf(stderr, "average-rhs-test: numerical average doesn't agree with analytical average\n");
    fprintf(stderr, "Analytic: %8g %8g %8g %8g %8g %8g %8g %8g\n",
            rhs[0], rhs[1], rhs[2], rhs[3], rhs[4], rhs[5], rhs[6],
            rhs[7]);
    fprintf(stderr, "Numerical: %8g %8g %8g %8g %8g %8g %8g %8g\n",
            numerical_rhs[0], numerical_rhs[1], numerical_rhs[2], numerical_rhs[3], numerical_rhs[4], numerical_rhs[5], numerical_rhs[6],
            numerical_rhs[7]);

    status = 1;
  }

  gsl_rng_free(rng);
  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);

  return status;
}
