#include"../src/rings.h"
#include<gsl/gsl_rng.h>
#include"checks.h"

int main() {
  double v[BODY_VECTOR_SIZE];
  body b;
  double vp[BODY_VECTOR_SIZE];
  double vc[CENTRAL_BODY_VECTOR_SIZE];
  central_body bc;
  double vcp[CENTRAL_BODY_VECTOR_SIZE];
  int i;
  gsl_rng *rng = NULL;
  const double eps = 1e-12;

  rng=gsl_rng_alloc(gsl_rng_ranlxd2);
  seed_random(rng);

  for (i = 0; i < BODY_VECTOR_SIZE; i++) {
    v[i] = gsl_rng_uniform(rng);
  }
  for (i = 0; i < CENTRAL_BODY_VECTOR_SIZE; i++) {
    vc[i] = gsl_rng_uniform(rng);
  }

  vector_to_body(v, &b);
  body_to_vector(&b, vp);

  vector_to_central_body(vc, &bc);
  central_body_to_vector(&bc, vcp);

  if (!check_vector_close(eps, eps, BODY_VECTOR_SIZE, v, vp)) return 1;
  if (!check_vector_close(eps, eps, CENTRAL_BODY_VECTOR_SIZE, vc, vcp)) return 2;

  gsl_rng_free(rng);

  return 0;
}
