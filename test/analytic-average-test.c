#include"../src/rings.h"
#include"checks.h"

#include<gsl/gsl_integration.h>
#include<gsl/gsl_rng.h>
#include<stdlib.h>

#define EPS 1e-12

typedef struct {
  double rp[3];
  body *b;
  int comp;
  double eps;
} params;

static double 
force(double E, void *vparams) {
  params *p = (params *)vparams;
  double r[3], v[3];
  double aspec[3], a[3];
  double e = norm(p->b->A);

  E_to_rv(p->b, E, r, v);
  softened_specific_acceleration(p->eps, p->rp, r, aspec);
  vscale(p->b->m, aspec, a);

  return a[p->comp]*(1.0 - e*cos(E))/(2.0*M_PI);
}

int
main() {
  body b1, b2;
  gsl_rng *rng;
  gsl_integration_workspace *ws;
  const size_t ws_size = 10000;
  int status = 0;
  double E1;
  double eps;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng != 0);
  seed_random(rng);

  ws = gsl_integration_workspace_alloc(ws_size);
  assert(ws != 0);

  init_random_body(rng, &b1, random_between(rng, 1e-3, 2e-3), random_between(rng, 1.0, 1.5));
  init_random_body(rng, &b2, random_between(rng, 1e-4, 2e-4), random_between(rng, 10.0, 15.0));

  eps = random_between(rng, 1e-2, 2e-2);

  for (E1 = 0.0; E1 < 2.0*M_PI && !status; E1 += 0.1) {
    double v1[3];
    params p;
    int i;
    double fanalytic[3];
    double fnumerical[3];

    E_to_rv(&b1, E1, p.rp, v1);
    p.b = &b2;
    p.eps = eps;

    force_averaged_unprimed(eps, p.rp, &b2, fanalytic);

    for (i = 0; i < 3; i++) {
      gsl_function f;
      double result, err;
      
      p.comp = i;
      f.function = (double (*)(double, void *))force;
      f.params = &p;
      
      gsl_integration_qag(&f, 0.0, 2.0*M_PI, EPS, EPS, ws_size, GSL_INTEG_GAUSS61, ws, &result, &err);

      fnumerical[i] = result;
    }

    if (!check_vector_close(EPS, EPS, 3, fanalytic, fnumerical)) {
      fprintf(stderr, "analytic-average-test: analytic and numerical average forces disagree.\n");
      status = 1;
    }
  }

  gsl_rng_free(rng);
  gsl_integration_workspace_free(ws);
  return status;
}
