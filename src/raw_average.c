#include"rings.h"

#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include<gsl/gsl_integration.h>

#define EPS 1e-12

typedef struct {
  double r1[3];
  double v1[3];
  double n1;
  const body *b1;
  const body *b2;
  double eps;
  int comp;
} inner_data;

static
double inner_fn(double E2, inner_data *data) {
  if (data->comp == BODY_M_INDEX || data->comp == BODY_Qp_INDEX || data->comp == BODY_I_INDEX || 
      data->comp == BODY_R_INDEX ||
      (data->comp >= BODY_SPIN_INDEX && data->comp < BODY_SPIN_INDEX + 3)) { 
    /* No change in these components. */
    return 0.0;
  } else {
    double r2[3], v2[3];
    double aspec[3], a[3];
    double e2 = get_e(data->b2);
    double fac = (1.0 - e2*cos(E2))/(2.0*M_PI);
    double n1 = data->n1;
    double a1 = data->b1->a;

    E_to_rv(data->b2, E2, r2, v2);

    softened_specific_acceleration(data->eps, data->r1, r2, aspec);
    vscale(data->b2->m, aspec, a);

    if (data->comp == BODY_a_INDEX) { 
      /* da1/dt */
      return fac * 2.0 * dot(data->v1, a) / (n1*n1*a1);
    } else if (BODY_L_INDEX <= data->comp && data->comp < BODY_L_INDEX+3) {
      /* dL/dt */
      int i = data->comp - BODY_L_INDEX;
      double rxa[3];

      cross(data->r1, a, rxa);

      return fac * rxa[i]/(n1*a1*a1);
    } else if (BODY_A_INDEX <= data->comp && data->comp < BODY_A_INDEX+3) {
      /* dA/dt */
      double adv = dot(a, data->v1);
      double rda = dot(data->r1, a);
      double rdv = dot(data->r1, data->v1);
      int i = data->comp - BODY_A_INDEX;

      return fac / (1.0 + data->b1->m) * (2.0*data->r1[i]*adv - data->v1[i]*rda - a[i]*rdv);
    } else {
      fprintf(stderr, "in inner_fn of raw_average_rhs with unknown component of derivative: %d", data->comp);
      abort();
    }
  }
}

typedef struct {
  const body *b1;
  const body *b2;
  gsl_integration_workspace *ws;
  size_t ws_limit;
  double eps;
  double epsabs;
  double epsrel;
  int comp;
} outer_data;

static
double outer_fn(double E1, outer_data *data) {
  inner_data idata;
  double result, error;
  gsl_function finner;
  double e1 = get_e(data->b1);
  double fac = (1.0 - e1*cos(E1))/(2.0*M_PI);

  E_to_rv(data->b1, E1, idata.r1, idata.v1);
  idata.b1 = data->b1;
  idata.b2 = data->b2;
  idata.n1 = mean_motion(data->b1);
  idata.eps = data->eps;
  idata.comp = data->comp;

  finner.function = (double (*)(double, void *))inner_fn;
  finner.params = &idata;

  gsl_integration_qag(&finner, 0.0, 2.0*M_PI, data->epsabs, data->epsrel, 
                      data->ws_limit, GSL_INTEG_GAUSS61, data->ws, &result, &error);

  return result*fac;
}

void
raw_average_rhs(const double eps, const body *b1, const body *b2, 
                gsl_integration_workspace *ws1, const size_t ws1_limit,
                gsl_integration_workspace *ws2, const size_t ws2_limit,
                const double epsabs, const double epsrel,
                double rhs[BODY_VECTOR_SIZE]) {
  outer_data odata;
  int i;
  gsl_function fouter;

  odata.b1 = b1;
  odata.b2 = b2;
  odata.ws = ws2;
  odata.ws_limit = ws2_limit;
  odata.eps = eps;
  odata.epsabs = epsabs;
  odata.epsrel = epsrel;

  fouter.function = (double (*)(double, void *))outer_fn;
  fouter.params = &odata;

  for (i = 0; i < BODY_VECTOR_SIZE; i++) {
    double result, error;
    
    odata.comp = i;

    gsl_integration_qag(&fouter, 0.0, 2.0*M_PI, epsabs, epsrel, ws1_limit, GSL_INTEG_GAUSS61, ws1, &result, &error);

    rhs[i] = result;
  }  
}
