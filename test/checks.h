#ifndef __CHECKS_H__
#define __CHECKS_H__

int
check_close(const double epsabs, const double epsrel, const double x, const double y);

int
check_vector_close(const double epsabs, const double epsrel, const double x[3], const double y[3]);

#endif /* __CHECKS_H__ */
