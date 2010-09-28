#include"rings.h"

#include<math.h>

#include<gsl/gsl_integration.h>

#define EPS 1e-12

typedef struct {
  double r1[3];
  double v1[3];
  double n1;
  body *b1;
  body *b2;
  double eps;
  int comp;
} inner_data;

static
double inner_fn(double E2, inner_data *data) {
  if (data->comp == 0) { /* No change in the mass. */
    return 0;
  } else {
    double r2[3], v2[3];
    double aspec[3], a[3];
    double e2 = norm(data->b2->A);
    double fac = (1.0 - e2*cos(E2))/(2.0*M_PI);

    E_to_rv(data->b2, E2, r2, v2);

    softened_specific_acceleration(data->eps, data->r1, r2, aspec);
    vscale(data->b2->m, aspec, a);

    if (data->comp == 1) { 
      /* da1/dt */
      return fac * 2.0 * dot(data->v1, a) / (data->n1*data->n1*data->b1->a);
    } else if (2 <= data->comp && data->comp < 5) {
      /* dL/dt */
      int i = data->comp - 2;
      double a2 = data->b1->a*data->b1->a;
      double rxa[3];

      cross(data->r1, a, rxa);

      return fac * rxa[i]/(data->n1*a2);
    } else if (5 <= data->comp && data->comp < 8) {
      /* dA/dt */
      double adv = dot(a, data->v1);
      double rda = dot(data->r1, a);
      double rdv = dot(data->r1, data->v1);
      int i = data->comp - 5;

      return fac / (1.0 + data->b1->m) * (2.0*data->r1[i]*adv - data->v1[i]*rda - a[i]*rdv);
    }
  }
}

typedef struct {
  body *b1;
  body *b2;
  gsl_integration_workspace *ws;
  size_t ws_limit;
  double eps;
  int comp;
} outer_data;

static
double outer_fn(double E1, outer_data *data) {
  inner_data idata;
  double result, error;
  gsl_function finner;
  double e1 = norm(data->b1->A);
  double fac = (1.0 - e1*cos(E1))/(2.0*M_PI);

  E_to_rv(data->b1, E1, idata.r1, idata.v1);
  idata.b1 = data->b1;
  idata.b2 = data->b2;
  idata.n1 = mean_motion(data->b1);
  idata.eps = data->eps;
  idata.comp = data->comp;

  finner.function = (double (*)(double, void *))inner_fn;
  finner.params = &idata;

  gsl_integration_qag(&finner, 0.0, 2.0*M_PI, EPS, EPS, data->ws_limit, GSL_INTEG_GAUSS61, data->ws, &result, &error);

  return result*fac;
}

void
raw_average_rhs(const double eps, const body *b1, const body *b2, 
                gsl_integration_workspace *ws1, const size_t ws1_limit,
                gsl_integration_workspace *ws2, const size_t ws2_limit,
                double rhs[BODY_VECTOR_SIZE]) {
  outer_data odata;
  int i;
  gsl_function fouter;

  odata.b1 = b1;
  odata.b2 = b2;
  odata.ws = ws2;
  odata.ws_limit = ws2_limit;
  odata.eps = eps;

  fouter.function = (double (*)(double, void *))outer_fn;
  fouter.params = &odata;

  for (i = 0; i < BODY_VECTOR_SIZE; i++) {
    double result, error;
    
    odata.comp = i;

    gsl_integration_qag(&fouter, 0.0, 2.0*M_PI, EPS, EPS, ws1_limit, GSL_INTEG_GAUSS61, ws1, &result, &error);

    rhs[i] = result;
  }  
}
