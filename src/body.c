#include"rings.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>
#include<gsl/gsl_errno.h>

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

double
get_e(const body *b) {
  double eA = norm(b->A);

  if (eA < 0.5) {
    /* If e is small, use b->A because it implements reflecting
       boundary condition at e = 0. */
    return eA;
  } else {
    /* If e is large, use b->L because it implements reflecting
       boundary condition around e = 1. */
    double ln = norm(b->L);
    double eL = sqrt(1.0 - ln*ln);

    return eL;
  }
}

int
E_to_rv(const body *b, const double E, double x[3], double v[3]) {
  double xhat[3], yhat[3], zhat[3];
  double e = get_e(b);

  if (e < 0.0 || e >= 1.0 || isnan(e)) {
    return GSL_EINVAL;
  }

  double n = mean_motion(b);
  double cosE = cos(E), sinE = sin(E);
  double eCosE = e*cosE;
  double sqrt1me2 = sqrt(1.0-e*e);
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

  return GSL_SUCCESS;
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
                        const double Omega, const double omega, const double spin[3], 
                        const double Qp, const double inertia, const double R) {
  b->m = m;
  b->a = a;
  b->Qp = Qp;
  b->I = inertia;
  b->R = R;
  
  b->L[0] = 0.0;
  b->L[1] = 0.0;
  b->L[2] = sqrt(1.0 - e*e);

  b->A[0] = e;
  b->A[1] = 0.0;
  b->A[2] = 0.0;

  memcpy(b->spin, spin, 3*sizeof(double));  

  rotate_to_orbit_frame(b->L, I, Omega, omega);
  rotate_to_orbit_frame(b->A, I, Omega, omega);
  rotate_to_orbit_frame(b->spin, I, Omega, omega);
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
  *e = get_e(b);

  if (b->L[1]*b->L[1] + b->L[0]*b->L[0] < 1e4*DBL_EPSILON*DBL_EPSILON) {
    /* If we are basically in the x-y plane, then Omega is undefined; use x-axis. */
    *Omega = 0.0;
  } else {
    *Omega = normalize_angle(lang + M_PI/2.0, 0.0, 2.0*M_PI);
  }

  my_omega = *Omega;

  *Omega *= 180.0/M_PI;

  asc_node[0] = cos(my_omega);
  asc_node[1] = sin(my_omega);
  asc_node[2] = 0.0;

  if (*e < 100.0*DBL_EPSILON) {
    /* Then e is so small we don't really have omega. */
    *omega = 0.0;
  } else {
    double cos_omega = dot(asc_node, b->A)/get_e(b);
    
    if (cos_omega > 1.0) cos_omega = 1.0;
    if (cos_omega < -1.0) cos_omega = -1.0;

    *omega = acos(dot(asc_node, b->A)/get_e(b));
    
    cross(asc_node, b->A, asc_node_x_A);
    if (dot(b->L, asc_node_x_A) < 0.0) { /* This means large argument of periapse. */
      *omega = 2.0*M_PI - *omega;
    }
    
    *omega = *omega*180.0/M_PI;
  }
}

void
body_to_vector(const body *b, double *v) {
  v[BODY_M_INDEX] = b->m;
  v[BODY_a_INDEX] = b->a;
  v[BODY_Qp_INDEX] = b->Qp;
  v[BODY_I_INDEX] = b->I;
  v[BODY_R_INDEX] = b->R;
  memcpy(v+BODY_L_INDEX, b->L, 3*sizeof(double));
  memcpy(v+BODY_A_INDEX, b->A, 3*sizeof(double));
  memcpy(v+BODY_SPIN_INDEX, b->spin, 3*sizeof(double));
}

void
vector_to_body(const double *v, body *b) {
  b->m = v[BODY_M_INDEX];
  b->a = v[BODY_a_INDEX];
  b->Qp = v[BODY_Qp_INDEX];
  b->I = v[BODY_I_INDEX];
  b->R = v[BODY_R_INDEX];
  memcpy(b->L, v+BODY_L_INDEX, 3*sizeof(double));
  memcpy(b->A, v+BODY_A_INDEX, 3*sizeof(double));
  memcpy(b->spin, v+BODY_SPIN_INDEX, 3*sizeof(double));
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
  double LNorm = norm(b->L);
  double ANorm = norm(b->A);

  if (LNorm > 0.0 && ANorm > 0.0) {
    unitize(b->L, zhat);
    unitize(b->A, xhat);
    cross(zhat, xhat, yhat);
  } else if (LNorm > 0.0) {
    /* A norm is too small. */
    double xhattmp[3], yhattmp[3];
    unitize(b->L, zhat);
    xhat[0] = 1.0; xhat[1] = 0.0; xhat[2] = 0.0;
    cross(zhat, xhat, yhattmp);
    cross(yhat, zhat, xhattmp);

    unitize(xhattmp, xhat);
    unitize(yhattmp, yhat);
  } else {
    /* L is too small. */
    double xhattmp[3], yhattmp[3];
    unitize(b->A, xhat);
    zhat[0] = 0.0; zhat[1] = 0.0; zhat[2] = 1.0;
    cross(zhat, xhat, yhattmp);
    cross(yhat, zhat, xhattmp);

    unitize(xhattmp, xhat);
    unitize(yhattmp, yhat);
  }
}

void
body_set_synchronous_spin(body *b) {
  double n = mean_motion(b);
  double spinHat[3];

  unitize(b->L, spinHat);
  vscale(n, spinHat, b->spin);
}

double
body_system_amd(const body bs[], const size_t n) {
  double Ltot[3];
  double Lnrm;
  double amd = 0.0;
  size_t i;

  memset(Ltot, 0, 3*sizeof(double));
  for (i = 0; i < n; i++) {
    double nrm = bs[i].m*sqrt(bs[i].a);
    Ltot[0] += nrm*bs[i].L[0];
    Ltot[1] += nrm*bs[i].L[1];
    Ltot[2] += nrm*bs[i].L[2];
  }
  Lnrm = norm(Ltot);

  amd = 0.0;
  for (i = 0; i < n; i++) {
    amd += bs[i].m*sqrt(bs[i].a)*(1 - dot(bs[i].L, Ltot)/Lnrm);
  }

  return amd;
}
