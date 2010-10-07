#include"rings.h"

#include<gsl/gsl_integration.h>
#include<gsl/gsl_odeiv.h>
#include<gsl/gsl_errno.h>

size_t
body_size_to_vector_size(const size_t nbodies) { return nbodies*BODY_VECTOR_SIZE; }

void
bodies_to_vector(const body bs[], const size_t nbodies, double y[]) {
  size_t i;

  for(i = 0; i < nbodies; i++) {
    body_to_vector(bs+i, y+i*BODY_VECTOR_SIZE);
  }
}

void
vector_to_bodies(const double y[], const size_t nbodies, body bs[]) {
  size_t i;

  for (i = 0; i < nbodies; i++) {
    vector_to_body(y+i*BODY_VECTOR_SIZE, bs+i);
  }
}

typedef struct {
  size_t nbodies;
  double eps;
  gsl_integration_workspace *ws;
  size_t ws_size;
  double epsabs;
  double epsrel;
} fparams;

static int
f(double t, const double y[], double dydt[], void *vparams) {
  fparams *p = (fparams *)vparams;
  size_t i;

  for (i = 0; i < p->nbodies; i++) {
    size_t j;
    body bi;

    vector_to_body(y+i*BODY_VECTOR_SIZE, &bi);

    for (j = 0; j < BODY_VECTOR_SIZE; j++) {
      dydt[i*BODY_VECTOR_SIZE+j] = 0.0;
    }

    for (j = 0; j < p->nbodies; j++) {
      if (j != i) {
        body bj;
        double rhs[BODY_VECTOR_SIZE];
        size_t k;

        vector_to_body(y+j*BODY_VECTOR_SIZE, &bj);
        
        average_rhs(p->eps, &bi, &bj, p->ws, p->ws_size, p->epsabs, p->epsrel, rhs);

        for (k = 0; k < BODY_VECTOR_SIZE; k++) {
          dydt[i*BODY_VECTOR_SIZE+k] += rhs[k];
        }
      }
    }
  }

  return GSL_SUCCESS;
}

int
evolve_system(gsl_odeiv_evolve *e, gsl_odeiv_control *con, gsl_odeiv_step *step, 
              double *t, const double t1, double *h, body bs[], double y[], 
              size_t nbodies, 
              gsl_integration_workspace *ws, const size_t ws_size,
              const double epsabs, const double epsrel,
              const double eps) {
  fparams p;
  gsl_odeiv_system sys;
  int status;

  p.eps = eps;
  p.epsabs = epsabs;
  p.epsrel = epsrel;
  p.ws = ws;
  p.ws_size = ws_size;
  p.nbodies = nbodies;

  sys.function = (int (*)(double, const double [], double [], void *))f;
  sys.jacobian = 0;
  sys.dimension = body_size_to_vector_size(nbodies);
  sys.params = &p;

  bodies_to_vector(bs, nbodies, y);

  status = gsl_odeiv_evolve_apply(e, con, step, &sys, t, t1, h, y);

  vector_to_bodies(y, nbodies, bs);

  return status;
}
