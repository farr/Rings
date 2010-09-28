#include"../src/rings.h"
#include"checks.h"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_integration.h>

static double
rtod(const double r) { return r*180.0/M_PI; }

static double
dtor(const double d) { return d*M_PI/180.0; }

int
main() {
  body earth, jupiter;
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  double rhs_earth[BODY_VECTOR_SIZE], rhs_jup[BODY_VECTOR_SIZE];
  int status = 0;
  double edot_earth, edot_jup;

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);

  assert(ws1 != 0);
  assert(ws2 != 0);

  init_body_from_elements(&earth, 
                          3e-6, 9.994971186791173E-01, 
                          1.644474019475725E-02,
                          1.644474019475725E-02,
                          6.745908819622363E+01,
                          3.734578297291855E+01);

  init_body_from_elements(&jupiter,
                          0.00095, 5.202363717344561E+00,
                          4.880540780991029E-02,
                          1.303833554481467E+00,
                          1.005116386574437E+02,
                          2.739290400833941E+02);

  raw_average_rhs(0.0, &earth, &jupiter, ws1, ws_size, ws2, ws_size, rhs_earth);
  raw_average_rhs(0.0, &jupiter, &earth, ws1, ws_size, ws2, ws_size, rhs_jup);

  edot_earth = vdot_to_vmagdot(earth.A, rhs_earth + BODY_A_INDEX);
  edot_jup = vdot_to_vmagdot(jupiter.A, rhs_jup + BODY_A_INDEX);

  fprintf(stderr, "edot (Earth) = %g\n", edot_earth);
  fprintf(stderr, "edot (Jupiter) = %g\n", edot_jup);

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
