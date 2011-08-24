#include "../src/rings.h"
#include <gsl/gsl_rng.h>
#include"checks.h"

int main() {
  body b1, b2;
  double rhs1[BODY_VECTOR_SIZE], rhs2[BODY_VECTOR_SIZE];
  gsl_rng *rng;
  double eps;
  const double epsabs = 1e-8;
  double ori1, ori2;
  double magdot1, magdot2;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);

  eps = gsl_rng_uniform(rng);

  init_random_body(rng, &b1, 1e-1*gsl_rng_uniform(rng), 1.0+gsl_rng_uniform(rng));
  init_random_body(rng, &b2, 2e-2*gsl_rng_uniform(rng), 2.0+gsl_rng_uniform(rng));

  average_rhs(eps, &b1, &b2, epsabs, rhs1);
  average_rhs(eps, &b2, &b1, epsabs, rhs2);

  magdot1 = dot(b1.A, &(rhs1[BODY_A_INDEX])) + dot(b1.L, &(rhs1[BODY_L_INDEX]));
  magdot2 = dot(b2.A, &(rhs2[BODY_A_INDEX])) + dot(b2.L, &(rhs2[BODY_L_INDEX]));

  if (!check_close(10.0*epsabs, 0.0, magdot1, 0.0)) {
    fprintf(stderr, "L1*L1dot + A1*A1dot != 0\n");
    return 1;
  }

  if (!check_close(10.0*epsabs, 0.0, magdot2, 0.0)) {
    fprintf(stderr, "L2*L2dot + A2*A2dot != 0\n");
    return 2;
  }

  ori1 = dot(b1.A, &(rhs1[BODY_L_INDEX])) + dot(b1.L, &(rhs1[BODY_A_INDEX]));
  ori2 = dot(b2.A, &(rhs2[BODY_L_INDEX])) + dot(b2.L, &(rhs2[BODY_A_INDEX]));

  if (!check_close(10.0*epsabs, 0.0, ori1, 0.0)) {
    fprintf(stderr, "L1*A1dot + A1*L1dot != 0\n");
    return 1;
  }

  if (!check_close(10.0*epsabs, 0.0, ori2, 0.0)) {
    fprintf(stderr, "L2*A2dot + A2*L2dot != 0\n");
    return 1;
  }
    
  gsl_rng_free(rng);

  return 0;
}
