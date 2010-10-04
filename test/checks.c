#include"checks.h"
#include<math.h>
#include<assert.h>

int
check_close(const double epsabs, const double epsrel, const double x, const double y) {
  double dx = fabs(x-y);
  double ave = 0.5*(fabs(x)+fabs(y));

  return dx <= epsabs + epsrel*ave;
}

int
check_vector_close(const double epsabs, const double epsrel, const size_t len, const double x[3], const double y[3]) {
  size_t i;
  for (i = 0; i < len; i++) {
    if (!check_close(epsabs, epsrel, x[i], y[i])) return 0;
  }

  return 1;
}

void
seed_random(gsl_rng *rng) {
  unsigned long int seed;
  FILE *devrandom;

  devrandom = fopen("/dev/random", "r");
  assert(devrandom != 0);

  if (fread(&seed, sizeof(unsigned long int), 1, devrandom) != 1) {
    fprintf(stderr, "%s at line %d: couldn't read seed from /dev/random", __FILE__, __LINE__);
    exit(1);
  }

  if (fclose(devrandom) != 0) {
    fprintf(stderr, "%s at line %d: could not close /dev/random", __FILE__, __LINE__);
    exit(1);
  }

  gsl_rng_set(rng, seed);
}

double
random_between(gsl_rng *rng, const double a, const double b) {
  return a + (b-a)*(gsl_rng_uniform(rng));
}

void
init_random_body(gsl_rng *rng, body *b, const double m, const double a) {
  b->m = m;
  init_body_from_elements(b, m, a,
                          random_between(rng, 0.0, 1.0),
                          random_between(rng, 0.0, 180.0),
                          random_between(rng, 0.0, 2.0*M_PI),
                          random_between(rng, 0.0, 2.0*M_PI));                          
}
