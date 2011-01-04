#include"../src/rings.h"
#include<assert.h>
#include<math.h>
#include<gsl/gsl_integration.h>

static double 
inner_inclination(const double m1, const double a1, const double e1, 
                  const double m2, const double a2, const double e2, const double I2) {
  double L1 = m1*a1*sqrt(1.0-e1*e1);
  double L2 = m2*a2*sqrt(1.0-e2*e2);
  double I1 = asin(L2*sin(I2)/L1);

  assert(!isnan(I1));

  return I1;
}

int main() {
  body b1, b2;
  const size_t ws_size = 10000;
  const double eps = 0.02;
  const double epsabs = 1e-8, epsrel = 1e-8;
  gsl_integration_workspace *ws1, *ws2;
  double analytic_rhs[BODY_VECTOR_SIZE], numerical_rhs[BODY_VECTOR_SIZE];
  const double m1 = 1e-3, m2 = 2.73e-4;
  const double a1 = 1.0, a2 = 6.32;
  const double elow = 0.0, emid = 0.3, ehigh = 0.999;
  const double Omega1 = 0.0, Omega2 = M_PI;
  const double omega1 = M_PI/3.0, omega2 = sqrt(2.0)*M_PI;
  const double I2 = 1.0*M_PI/180.0; /* Outer inclination 1 degree. */
  double I1; /* Adjusts to have Ltot = Lz_tot */

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != NULL);
  assert(ws2 != NULL);  

  /* Low eccentricity */
  I1 = inner_inclination(m1, a1, elow, m2, a2, ehigh, I2);
  init_body_from_elements(&b1, m1, a1, elow, I1, Omega1, omega1);
  init_body_from_elements(&b2, m2, a2, emid, I2, Omega2, omega2);

  /* RHS on b1. */
  average_rhs(eps, &b1, &b2, epsabs, analytic_rhs);
  raw_average_rhs(eps, &b1, &b2, ws1, ws_size, ws2, ws_size, epsabs/10.0, epsrel/10.0, numerical_rhs);

  if (!check_vector_close(10*epsabs, 10*epsrel, BODY_VECTOR_SIZE, analytic_rhs, numerical_rhs)) {
    fprintf(stderr, "  At low eccentricity, analytic and numerical RHSs disagree.\n");
    return 1;
  }

  /* RHS on b2. */
  average_rhs(eps, &b2, &b1, epsabs, analytic_rhs);
  raw_average_rhs(eps, &b2, &b1, ws1, ws_size, ws2, ws_size, epsabs/10.0, epsrel/10.0, numerical_rhs);

  if (!check_vector_close(10*epsabs, 10*epsrel, BODY_VECTOR_SIZE, analytic_rhs, numerical_rhs)) {
    fprintf(stderr, "  At low eccentricity, analytic and numerical RHSs disagree.\n");
    return 1;
  }


  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  
  return 0;
}
