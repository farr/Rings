#include"../src/rings.h"
#include"checks.h"

#include<gsl/gsl_rng.h>
#include<string.h>

int main() {
  const double RScale = 0.1;
  const double AScale = 1.0;
  const double IScale = 1e-2;
  const double MScale = 1.0;
  const double QScale = 1.0;

  const double eps = 1e-12;

  gsl_rng *rng;

  double SpinSun[3];
  double QSun, ISun, MSun, RSun;

  double zerov[3] = {0.0, 0.0, 0.0};

  body b;
  central_body bc;

  double n;
  double a2;
  double adot;
  double ndot;
  double mu;

  double body_rhs[BODY_VECTOR_SIZE];
  double sun_rhs[3];

  double Ldot[3], hdot[3], bLdot[3], sLdot[3];

  int i;

  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  assert(rng);

  seed_random(rng);

  init_body_from_elements(&b, MScale*gsl_rng_uniform(rng), AScale*gsl_rng_uniform(rng), gsl_rng_uniform(rng),
                          180.0*gsl_rng_uniform(rng), 360.0*gsl_rng_uniform(rng), 360.0*gsl_rng_uniform(rng), 
                          zerov, QScale*gsl_rng_uniform(rng), IScale*gsl_rng_uniform(rng), RScale*gsl_rng_uniform(rng));

  random_vector(rng, b.spin, 1.0);
  random_vector(rng, SpinSun, 1.0);
  
  QSun = QScale*gsl_rng_uniform(rng);
  ISun = IScale*gsl_rng_uniform(rng);
  RSun = RScale*gsl_rng_uniform(rng);

  init_central_body(&bc, QSun, ISun, RSun, SpinSun);

  tidal_rhs(&b, &bc, body_rhs, sun_rhs);

  n = mean_motion(&b);
  a2 = b.a*b.a;

  adot = body_rhs[BODY_a_INDEX];
  ndot = -3.0/2.0*adot/b.a*n;

  mu = b.m/(1.0+b.m);

  for (i = 0; i < 3; i++) {
    hdot[i] = mu*n*a2*body_rhs[BODY_L_INDEX+i] + 2.0*mu*n*b.a*adot*b.L[i] + mu*ndot*a2*b.L[i];
    bLdot[i] = b.I*body_rhs[BODY_SPIN_INDEX+i];
    sLdot[i] = ISun*sun_rhs[i];
    Ldot[i] = hdot[i] + bLdot[i] + sLdot[i];
  }

  if (!check_close(eps, eps, norm(Ldot), 0.0)) {
    fprintf(stderr, "Ldot = {%g, %g, %g} not close to zero!\n",
            Ldot[0], Ldot[1], Ldot[2]);
    fprintf(stderr, "mu*hdot = {%g, %g, %g}\n", hdot[0], hdot[1], hdot[2]);
    fprintf(stderr, "Body LDot = {%g, %g, %g}\n", bLdot[0], bLdot[1], bLdot[2]);
    fprintf(stderr, "Sun  LDot = {%g, %g, %g}\n", sLdot[0], sLdot[1], sLdot[2]);
    gsl_rng_free(rng);
    exit(1);
  }

  if (!check_close(eps, eps, dot(b.L, body_rhs + BODY_L_INDEX) + dot(b.A, body_rhs + BODY_A_INDEX), 0.0)) {
    fprintf(stderr, "A^2 + L^2 = 1.0 magnitude constraint violated!\n");
    exit(1);
  }

  if (!check_close(eps, eps, dot(b.L, body_rhs + BODY_A_INDEX) + dot(b.A, body_rhs + BODY_L_INDEX), 0.0)) {
    fprintf(stderr, "A*L = 0 constraint violated!\n");
    exit(1);
  }

  gsl_rng_free(rng);
  return 0;
}
