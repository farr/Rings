#include"rings.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>

static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);

double
mean_motion(const body *b) {
  double a = b->a;
  double m = b->m;
  double mm = sqrt((1.0+m)/(a*a*a));
  return mm;
}

double
period(const body *b) {
  return 2.0*M_PI / mean_motion(b);
}

void
E_to_rv(const body *b, const double E, double x[3], double v[3]) {
  double xhat[3], yhat[3], zhat[3];
  double e = norm(b->A);
  double n = mean_motion(b);
  double cosE = cos(E), sinE = sin(E);
  double eCosE = e*cosE;
  double sqrt1me2 = norm(b->L);
  double denom = 1.0 - eCosE;
  double a = b->a;
  double xc = a*(cosE - e);
  double yc = a*sqrt1me2*sinE;
  double r = a*(1.0 - eCosE);
  double na2 = n*a*a;
  double xdot = -na2*sinE/r;
  double ydot = na2*sqrt1me2*cosE/r;

  int i;

  /* The reason for this if statement is to handle constructing a
     consistent coordinate system when e ~= 0 or e ~= 1. 

     In the former case, we have no real handle on xhat (because b->A
     ~= 0).  So, we get zhat from b->L.  Then we cross {1,0,0} or
     {0,1,0} with zhat (whichever is more orthogonal to zhat to begin
     with) to produce xhat in the orbital plane.  yhat then follows
     from yhat = zhat x xhat.

     In the latter case, we have no real handle on zhat (because b->L
     ~= 0).  So, we get xhat from b->A.  Then we cross either {0,1,0}
     or {0,0,1} with xhat (whichever is more orthogonal to begin with)
     to produce a yhat orthogonal to xhat.  zhat then follows from
     zhat = xhat x yhat. */
  if (e < 100.0*DBL_EPSILON) {
    double Ldx, Ldy;
    const double x[] = {1.0, 0.0, 0.0}, y[] = {0.0, 1.0, 0.0};
    double xhat_temp[3];

    unitize(b->L, zhat);

    Ldx = dot(zhat, x);
    Ldy = dot(zhat, y);

    if (Ldx > Ldy) {
      /* If L is more along x than y, then use y to cross with z to get xhat. */
      cross(y, zhat, xhat_temp);
      unitize(xhat_temp, xhat);
    } else {
      /* L is more along y than x, so we can use x to cross with z to get xhat. */
      cross(x, zhat, xhat_temp);
      unitize(xhat_temp, xhat);
    }

    cross(zhat, xhat, yhat);
  } else if (sqrt1me2 < 100.0*DBL_EPSILON) {
    double Ady, Adz;
    const double y[] = {0.0, 1.0, 0.0}, z[] = {0.0, 0.0, 1.0};
    double yhat_temp[3];

    unitize(b->A, xhat);

    Ady = dot(xhat, y);
    Adz = dot(xhat, z);

    if (Ady > Adz) {
      /* A is more along y, so use z to cross with xhat to get yhat. */
      cross(z, xhat, yhat_temp);
      unitize(yhat_temp, yhat);
    } else {
      /* A is more along x, so use y for yhat. */
      cross(y, xhat, yhat_temp);
      unitize(yhat_temp, yhat);
    }

    cross(xhat, yhat, zhat);
  } else {
    unitize(b->L, zhat);
    unitize(b->A, xhat);
    cross(zhat, xhat, yhat);
  }

  for (i = 0; i < 3; i++) {
    x[i] = xc*xhat[i] + yc*yhat[i];
    v[i] = xdot*xhat[i] + ydot*yhat[i];
  }
}

static void
rotate_to_orbit_frame(double v[3], const double I, const double Omega, const double omega) {
  double temp[3];

  memcpy(temp, v, 3*sizeof(double));

  rotate_z(temp, omega*M_PI/180.0, v); /* Periapse to omega in x-y plane. */
  rotate_x(v, I*M_PI/180.0, temp); /* Inclined orbit; line of nodes still x-axis. */
  rotate_z(temp, Omega*M_PI/180.0, v); /* Rotate line of nodes to correct location. */
  /* Correct vector stored in v. */
}

void
init_body_from_elements(body *b, 
                        const double m, const double a, const double e, const double I,
                        const double Omega, const double omega) {
  b->m = m;
  b->a = a;
  
  b->L[0] = 0.0;
  b->L[1] = 0.0;
  b->L[2] = sqrt(1.0 - e*e);

  b->A[0] = e;
  b->A[1] = 0.0;
  b->A[2] = 0.0;

  rotate_to_orbit_frame(b->L, I, Omega, omega);
  rotate_to_orbit_frame(b->A, I, Omega, omega);
}

static double
normalize_angle(const double ang, const double min, const double max) {
  assert(fabs(max - min - 2.0*M_PI) < 1e-8);

  if (ang < min) {
    return max - fmod(max - ang, 2.0*M_PI);
  } else if (ang > max) {
    return min + fmod(ang - min, 2.0*M_PI);
  } else {
    return ang;
  }    
}

void
elements_from_body(const body *b,
                   double *e, double *I, double *Omega, double *omega) {
  
  double Lmag = norm(b->L);
  double lang = atan2(b->L[1], b->L[0]); /* Between -Pi and Pi */
  double asc_node[3];
  double my_omega;
  double asc_node_x_A[3];

  *I = acos(b->L[2]/Lmag)*180.0/M_PI;
  *e = norm(b->A);
  *Omega = normalize_angle(lang + M_PI/2.0, 0.0, 2.0*M_PI);

  my_omega = *Omega;

  *Omega *= 180.0/M_PI;

  asc_node[0] = cos(my_omega);
  asc_node[1] = sin(my_omega);
  asc_node[2] = 0.0;

  *omega = acos(dot(asc_node, b->A)/norm(b->A));

  cross(asc_node, b->A, asc_node_x_A);
  if (dot(b->L, asc_node_x_A) < 0.0) { /* This means large argument of periapse. */
    *omega = 2.0*M_PI - *omega;
  }

  *omega = *omega*180.0/M_PI;
}

void
body_to_vector(const body *b, double *v) {
  v[0] = b->m;
  v[1] = b->a;
  memcpy(v+2, b->L, 3*sizeof(double));
  memcpy(v+5, b->A, 3*sizeof(double));
}

void
vector_to_body(const double *v, body *b) {
  b->m = v[0];
  b->a = v[1];
  memcpy(b->L, v+2, 3*sizeof(double));
  memcpy(b->A, v+5, 3*sizeof(double));
}

void
body_instantaneous_rhs(const double eps,
                       const body *b1, const double E1,
                       const body *b2, const double E2,
                       double dbdt[BODY_VECTOR_SIZE]) {
  double aspec[3], a[3];
  double r1[3], v1[3], r2[3], v2[3];
  double n1 = mean_motion(b1);
  double rxa[3], adv, rda, rdv;
  int i;

  E_to_rv(b1, E1, r1, v1);
  E_to_rv(b2, E2, r2, v2);

  softened_specific_acceleration(eps, r1, r2, aspec);

  vscale(b2->m, aspec, a);
  
  cross(r1, a, rxa);
  adv = dot(a, v1);
  rda = dot(a, r1);
  rdv = dot(r1, v1);

  dbdt[0] = 0.0; /* No change in m. */
  dbdt[1] = 2.0 / (n1*n1*b1->a) * dot(v1, a);
  for (i = 0; i < 3; i++) {
    dbdt[2+i] = rxa[i]/(n1*n1*b1->a);
    dbdt[5+i] = 1.0/(1.0+b1->m)*(2.0*adv*r1[i] - rda*v1[i] - rdv*a[i]);
  }
}

void
body_derivs_to_orbital_elements(const body *b, const double dbdt[BODY_VECTOR_SIZE], double deldt[NELEMENTS]) {
  double Lmag = norm(b->L);
  double Lx = b->L[0], Ly = b->L[1], Lz = b->L[2];
  double cosI = Lz/Lmag;
  double sinI = sqrt(1.0 - cosI*cosI);
  double Lmag_dot = vdot_to_vmagdot(b->L, dbdt + BODY_L_INDEX);

  deldt[A_INDEX] = dbdt[BODY_a_INDEX];

  deldt[E_INDEX] = vdot_to_vmagdot(b->A, dbdt + BODY_A_INDEX);

  deldt[I_INDEX] = 1.0/sinI*(Lz/Lmag * Lmag_dot/Lmag - dbdt[BODY_L_INDEX+2]/Lmag)*180.0/M_PI;

  deldt[OMEGA_INDEX] = 1.0/(1.0+Lx*Lx/(Ly*Ly))*(dbdt[BODY_L_INDEX+1]/Lx - Ly/(Lx*Lx)*dbdt[BODY_L_INDEX])*180.0/M_PI;

  deldt[oMEGA_INDEX] = 0.0/0.0;  /* Not implemented yet. */
}

void
body_coordinate_system(const body *b, double xhat[3], double yhat[3], double zhat[3]) {
  unitize(b->L, zhat);
  unitize(b->A, xhat);
  cross(zhat, xhat, yhat);
}
