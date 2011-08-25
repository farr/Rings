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
random_vector(gsl_rng *rng, double v[3], const double scale) {
  double vol = random_between(rng, 0.0, 4.0/3.0*M_PI);
  double phi = random_between(rng, 0.0, 2.0*M_PI);
  double cos_theta = random_between(rng, -1.0, 1.0);
  double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  double r = pow(vol*3.0/(4.0*M_PI), 1.0/3.0);

  v[0] = scale*r*cos(phi)*sin_theta;
  v[1] = scale*r*sin(phi)*sin_theta;
  v[2] = scale*r*cos_theta;
}

void
init_random_body(gsl_rng *rng, body *b, const double m, const double a, const double Qp, const double inertia, const double R) {
  double n = sqrt((1.0+m)/(a*a*a));
  double spin[3];

  random_vector(rng, spin, n);
  
  b->m = m;
  init_body_from_elements(b, m, a,
                          random_between(rng, 0.0, 1.0),
                          random_between(rng, 0.0, 180.0),
                          random_between(rng, 0.0, 360.0),
                          random_between(rng, 0.0, 360.0),
                          spin, Qp, inertia, R);                          
}
