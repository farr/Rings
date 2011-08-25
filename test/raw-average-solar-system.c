#include"../src/rings.h"
#include"checks.h"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_integration.h>

#define EPS 1e-10

static double
rtod(const double r) { return r*180.0/M_PI; }

static double
dtor(const double d) { return d*M_PI/180.0; }

int
main() {
  const int nplanet = 9; /* Yeah, we're still calling Pluto a planet. */
  body planets[nplanet];  
  const char *names[] = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
  gsl_integration_workspace *ws1, *ws2;
  const size_t ws_size = 10000;
  double rhs[nplanet][BODY_VECTOR_SIZE];
  const double Qp = 0.0;
  const double inertia = 0.0;
  double spin[3] = {0.0, 0.0, 0.0};
  int status = 0;
  double edot_earth;
  int i;

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);

  assert(ws1 != 0);
  assert(ws2 != 0);

  init_body_from_elements(planets+0, 1.66e-7, 0.387, 0.206, 7.005, 48.332, 77.456, spin, Qp, inertia);
  init_body_from_elements(planets+1, 2.448e-6, 0.723, 0.007, 3.395, 76.68,  131.532, spin, Qp, inertia);
  init_body_from_elements(planets+2, 3e-6, 1.000, 0.017, 5e-5, 348.739, 102.947, spin, Qp, inertia);
  init_body_from_elements(planets+3, 3.227e-7, 1.524, 0.093, 1.85, 49.58,   336.04, spin, Qp, inertia);
  init_body_from_elements(planets+4, 0.000955, 5.203, 0.0484, 1.31, 100.56, 14.75, spin, Qp, inertia);
  init_body_from_elements(planets+5, 0.000286, 9.537, 0.0542, 2.484, 113.715, 92.432, spin, Qp, inertia);
  init_body_from_elements(planets+6, 4.37e-5, 19.19, 0.047, 0.769, 74.23, 170.964, spin, Qp, inertia);
  init_body_from_elements(planets+7, 5.15e-5, 30.069, 0.0086, 1.769, 131.721, 44.971, spin, Qp, inertia);
  init_body_from_elements(planets+8, 6.39e-9, 39.482, 0.2488, 17.142, 110.303, 224.07, spin, Qp, inertia);

  for(i = 0; i < nplanet; i++) {
    int j;
    double deldt[NELEMENTS];
    
    for (j = 0; j < BODY_VECTOR_SIZE; j++) {
      rhs[i][j] = 0.0;
    }

    for (j = 0; j < nplanet; j++) {
      if (j != i) {
        double rhs_temp[BODY_VECTOR_SIZE];
        int k;

        raw_average_rhs(0.0, planets+i, planets+j, ws1, ws_size, ws2, ws_size, EPS, EPS, rhs_temp);

        for(k = 0; k < BODY_VECTOR_SIZE; k++) {
          rhs[i][k] += rhs_temp[k]; 
        }
      }
    }

    body_derivs_to_orbital_elements(planets+i, rhs[i], deldt);

    for (j = 0; j < NELEMENTS; j++) {
      deldt[j] = deldt[j] * (period(planets+2)*100.0); /* Convert from our time units to 100*years. */
    }
    /* fprintf(stderr, "Planet %10s derivatives: %9.3g %9.3g %9.3g %9.3g %9.3g\n", */
    /*         names[i], deldt[A_INDEX], deldt[E_INDEX],  */
    /*         deldt[I_INDEX]*3600.0, deldt[OMEGA_INDEX]*3600.0, deldt[oMEGA_INDEX]*3600.0); */
  }

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
