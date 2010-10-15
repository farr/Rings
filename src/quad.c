#include"rings.h"

#include<gsl/gsl_errno.h>

#include<string.h>

#define MAXSUBDIV 10000

static int
cleanup(const int status, double result[BODY_VECTOR_SIZE]) {
  memset(result, 0, BODY_VECTOR_SIZE*sizeof(double));
  return status;
}

static int
should_continue(const double integral[BODY_VECTOR_SIZE], const double eps) {
  return integral[BODY_a_INDEX] > eps;
}

int
quad(const integrand f, void *fdata, const double a, const double b, const double eps, 
     double result[BODY_VECTOR_SIZE]) {
  double temp[BODY_VECTOR_SIZE];
  int status;
  int i;
  int nsubdiv = 1;

  memset(result, 0, BODY_VECTOR_SIZE*sizeof(double));

  status = f(a, fdata, temp);
  if (status != GSL_SUCCESS) {
    return cleanup(status, result);
  }

  for(i = 0; i < BODY_VECTOR_SIZE; i++) {
    result[i] += 0.5*(b-a)*temp[i];
  }

  status = f(b, fdata, temp);
  if (status != GSL_SUCCESS) {
    return cleanup(status, result);
  }

  for(i = 0; i < BODY_VECTOR_SIZE; i++) {
    result[i] += 0.5*(b-a)*temp[i];
  }

  do {
    double h;

    nsubdiv *= 2;
    if (nsubdiv > MAXSUBDIV) {
      cleanup(GSL_EMAXITER, result);
    }

    h = (b-a)/nsubdiv;

    for (i = 0; i < BODY_VECTOR_SIZE; i++) {
      result[i] *= 0.5;
    }

    for (i = 1; i < MAXSUBDIV; i += 2) {
      int j;

      f(a + h*i, fdata, temp);

      for (j = 0; j < BODY_VECTOR_SIZE; j++) {
        result[j] += h*temp[j];
      }
    }
  } while(should_continue(result, eps));

  return GSL_SUCCESS;
}
