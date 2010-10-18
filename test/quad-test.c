#include"../src/rings.h"

#include<math.h>

#include<gsl/gsl_errno.h>

static int
quad_test_rhs(const double x, void *data, double result[BODY_VECTOR_SIZE]) {
  memset(result, 0, BODY_VECTOR_SIZE*sizeof(double));

  result[0] = exp(x);
  result[1] = x*x*x + x*x + x + 1 - 79.31361854897525; /* Integrates to zero over 2Pi, as with f*v. */
  result[2] = sin(x);
  result[3] = log(x);
  result[4] = sqrt(x);
  result[5] = pow(x, 3.8);
  result[6] = cosh(x);
  result[7] = sinh(x);

  return GSL_SUCCESS;
}

int
main() {
  double result[BODY_VECTOR_SIZE];
  const double true_result[] = 
    {534.4916555247647, 0, 0, 5.264536872885934, 10.49973996381495,
     1412.609526431630, 267.7448940410165, 266.7467614837482};
  int status, i;

  if (BODY_VECTOR_SIZE < 8) {
    fprintf(stderr, "quad-test: BODY_VECTOR_SIZE too small!\n");
    return 1;
  }

  status = quad(&quad_test_rhs, NULL, 0.0, 2.0*M_PI, 1e-10, result);
  if (status != GSL_SUCCESS) {
    fprintf(stderr, "quad-test: quad failed with status %d\n", status);
  }

  for (i = 0; i < 8; i++) {
    if (!check_close(1e-4, 1e-4, result[i], true_result[i])) {
      fprintf(stderr, "quad-test: error in component %d: got %g, expected %g\n", i, result[i], true_result[i]);
      return 2;
    }
  }
  
  return 0;
}
