#include "../src/rings.h"
#include<assert.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include"checks.h"

int main() {
  body b1, b2;
  const size_t ws_size = 100000;
  const double eps = 5e-2;
  const double epsabs = 1e-11, epsrel = 1e-11;
  gsl_integration_workspace *ws1, *ws2;
  double analytic_rhs[BODY_VECTOR_SIZE], numerical_rhs[BODY_VECTOR_SIZE];
  const double m1 = 1e-3, m2 = 2.73e-2;
  const double a = 1.0;
  const double e = 0.3;
  const double Omega = 0.0;
  const double omega = 0.0;
  const double I = 90.0;
  const double spin[] = {0.0, 0.0, 0,0};

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != NULL);
  assert(ws2 != NULL);

  init_body_from_elements(&b1, m1, a, e, 0.0, Omega, omega, spin, 0.0, 0.0, 0.0);
  init_body_from_elements(&b2, m2, a, e, I, Omega, omega, spin, 0.0, 0.0, 0.0);

  average_rhs(eps, &b1, &b2, epsabs, analytic_rhs);
  raw_average_rhs(eps, &b1, &b2, ws1, ws_size, ws2, ws_size, epsabs, epsrel, numerical_rhs);

  if (!check_vector_close(1e7*epsabs, 1e7*epsrel, BODY_VECTOR_SIZE, analytic_rhs, numerical_rhs)) {
    int i;
    fprintf(stderr, "  Numerical and analytic RHS's disagree.\n");
    for (i = 0; i < BODY_VECTOR_SIZE; i++) {
      fprintf(stderr, "    Component %d: %g (analytic) vs %g (numerical)\n", 
              i, analytic_rhs[i], numerical_rhs[i]);
    }
    return 1;
  }

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);

  return 0;
}
