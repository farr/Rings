#include"rings.h"

#include<gsl/gsl_odeiv.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

typedef struct {
  size_t nbodies;
  gsl_integration_workspace *ws1;
  size_t ws1_size;
  gsl_integration_workspace *ws2;
  size_t ws2_size;
  double eps;
} system_params;

static int
rhs(double t, const double y[], double dydt[], void *vparams) {
  system_params *params = (system_params *)vparams;
  body b1, b2;
  int i;

  for (i = 0; i < params->nbodies; i++) {
    int j;
    int i0 = i*BODY_VECTOR_SIZE;
    double rhs[BODY_VECTOR_SIZE];

    vector_to_body(y+i0, &b1);
    for (j = 0; j < BODY_VECTOR_SIZE; j++) {
      dydt[j+i0] = 0.0;
    }

    for (j = 0; j < params->nbodies; j++) {
      if (j != i) {
        int j0 = j*BODY_VECTOR_SIZE;
        int k;

        vector_to_body(y+j0, &b2);

        raw_average_rhs(params->eps, &b1, &b2, 
                        params->ws1, params->ws1_size,
                        params->ws2, params->ws2_size,
                        rhs);

        for (k = 0; k < BODY_VECTOR_SIZE; k++) {
          dydt[i0+k] += rhs[k];
        }
      }
    }
  }

  return GSL_SUCCESS;
}
