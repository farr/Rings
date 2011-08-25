#include"../src/rings.h"
#include<gsl/gsl_rng.h>
#include"checks.h"

int main() {
  body b1, b2;
  gsl_rng *rng;
  double eps;
  const double epsabs = 1e-8;
  double rhs1[BODY_VECTOR_SIZE], rhs2[BODY_VECTOR_SIZE];
  double Ldot[3];
  double n1, n2;
  double m1, m2;
  double a1, a2;
  int i;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);

  eps = gsl_rng_uniform(rng);

  a1 = 1.0+gsl_rng_uniform(rng);
  a2 = 2.0+gsl_rng_uniform(rng);
  m1 = 1e-3*gsl_rng_uniform(rng);
  m2 = 1e-4*gsl_rng_uniform(rng);

  init_random_body(rng, &b1, m1, a1, 0.0, 0.0);
  init_random_body(rng, &b2, m2, a2, 0.0, 0.0);

  n1 = mean_motion(&b1);
  n2 = mean_motion(&b2);

  average_rhs(eps, &b1, &b2, epsabs, rhs1);
  average_rhs(eps, &b2, &b1, epsabs, rhs2);

  for (i = 0; i < 3; i++) {
    Ldot[i] = m1*n1*a1*a1*rhs1[BODY_L_INDEX+i] + m2*n2*a2*a2*rhs2[BODY_L_INDEX+i];
  }

  if (!check_close(10.0*epsabs, 0.0, norm(Ldot), 0.0)) {
    fprintf(stderr, "Ldot != 0, instead {%g, %g, %g}\n", Ldot[0], Ldot[1], Ldot[2]);
    return 1;
  }

  gsl_rng_free(rng);
  return 0;
}
