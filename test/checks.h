#ifndef __CHECKS_H__
#define __CHECKS_H__

#include<gsl/gsl_rng.h>

int
check_close(const double epsabs, const double epsrel, const double x, const double y);

int
check_vector_close(const double epsabs, const double epsrel, const double x[3], const double y[3]);

/* Random utilities. */
void
seed_random(gsl_rng *rng);

double
random_between(gsl_rng *rng, const double a, const double b);

#endif /* __CHECKS_H__ */
