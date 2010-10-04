/* Murray, C. D. and Dermott, S. F.  Solar System Dynamics.  Cambridge
   University Press, Cambridge, UK (2004).  p 282. */

#include"../src/rings.h"
#include"checks.h"

#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_integration.h>

#define EPS 1e-10

int
main() {
  body saturn, jupiter;
  const size_t ws_size = 10000;
  gsl_integration_workspace *ws1, *ws2;
  int status = 0;
  double sat_rhs[BODY_VECTOR_SIZE], jup_rhs[BODY_VECTOR_SIZE];
  double sat_deltdt[NELEMENTS], jup_deltdt[NELEMENTS];

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != 0);
  assert(ws2 != 0);
  
  init_body_from_elements(&jupiter, 9.54786e-4, 5.202545, 0.0474622, 1.30667, 100.0381, 13.983865);
  init_body_from_elements(&saturn, 2.85837e-4, 9.554841, 0.0575481, 2.48795, 113.1334, 88.719425);

  raw_average_rhs(0.0, &jupiter, &saturn, ws1, ws_size, ws2, ws_size, EPS, EPS, jup_rhs);
  raw_average_rhs(0.0, &saturn, &jupiter, ws1, ws_size, ws2, ws_size, EPS, EPS, sat_rhs);

  body_derivs_to_orbital_elements(&jupiter, jup_rhs, jup_deltdt);
  body_derivs_to_orbital_elements(&saturn, sat_rhs, sat_deltdt);

  /* fprintf(stderr, "Jupiter edot (yr^-1): %g, Idot (rad/yr): %g\n", */
  /*         jup_deltdt[E_INDEX]*2.0*M_PI, jup_deltdt[I_INDEX]*2.0*M_PI*M_PI/180.0); */
  /* fprintf(stderr, "Saturn  edot (yr^-1): %g, Idot (rad/yr): %g\n", */
  /*         sat_deltdt[E_INDEX]*2.0*M_PI, sat_deltdt[I_INDEX]*2.0*M_PI*M_PI/180.0); */

  if (fabs(jup_deltdt[E_INDEX]*2.0*M_PI - 1.2920470269576777e-6) >= 5e-7) {
    fprintf(stderr, "raw-average-jupiter-saturn: bad Jupiter edot\n");
    status = 1;
  }

  if (fabs(jup_deltdt[I_INDEX]*2.0*M_PI*M_PI/180.0 - -3.501018965845696e-7) >= 5e-8) {
    fprintf(stderr, "raw-average-jupiter-saturn: bad Jupiter Idot\n");
    status = 2;
  }

  if (fabs(sat_deltdt[E_INDEX]*2.0*M_PI - -2.620812070178865e-6) >= 5e-7) {
    fprintf(stderr, "raw-average-jupiter-saturn: bad Saturn edot\n");
    status = 3;
  }

  if (fabs(sat_deltdt[I_INDEX]*2.0*M_PI*M_PI/180.0 - 4.5322043030781326e-7) >= 5e-8) {
    fprintf(stderr, "raw-average-jupiter-saturn: bad Saturn Idot\n");
    status = 4;
  }

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
