#include"rings.h"

#include<stdio.h>

int
read_body(FILE *stream, body *b) {
  int nscan;

  nscan = fscanf(stream, " %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg ",
                 &(b->m), &(b->a),
                 &(b->Qp), &(b->I),
                 &(b->R),
                 b->L, b->L+1, b->L+2,
                 b->A, b->A+1, b->A+2,
                 b->spin, b->spin+1, b->spin+2);

  if (nscan == BODY_VECTOR_SIZE) {
    return 1;
  } else {
    return 0;
  }
}

int
read_central_body(FILE *stream, central_body *bc) {
  double Qp, I, R, Omega[3];
  int nscan;

  nscan = fscanf(stream, " %lg %lg %lg %lg %lg %lg ",
                 &Qp, &I, &R, &(Omega[0]), &(Omega[1]), &(Omega[2]));

  if (nscan != CENTRAL_BODY_VECTOR_SIZE) return 0;

  init_central_body(bc, Qp, I, R, Omega);

  return 1;
}

int
read_body_from_elements(FILE *stream, body *b) {
  double m, a, e, I, Omega, omega, Qp, inertia, R, spin[3];
  int nscan;

  nscan = fscanf(stream, " %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg ",
                 &m, &a,
                 &e, &I,
                 &Omega, &omega,
                 &Qp, &inertia, &R,
                 spin, spin+1, spin+2);

  if (nscan != BODY_VECTOR_SIZE - 2) return 0;

  init_body_from_elements(b, m, a, e, I, Omega, omega, spin, Qp, inertia, R);

  return 1;
}

int
read_body_bin(FILE *stream, body *b) {
  size_t nread;

  nread = fread(b, sizeof(body), 1, stream);

  return nread;
}

int
read_central_body_bin(FILE *stream, central_body *bc) {
  size_t nread;

  nread = fread(bc, sizeof(central_body), 1, stream);

  return nread;
}

int
write_body(FILE *stream, const body *b) {
  int nchar;
  
  nchar = fprintf(stream, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n",
                  b->m, b->a, b->Qp, b->I, b->R, b->L[0], b->L[1], b->L[2],
                  b->A[0], b->A[1], b->A[2],
                  b->spin[0], b->spin[1], b->spin[2]);

  if (nchar > 0) { /* Assume that any written characters == total success! */
    return 0;
  } else {
    return 1;
  }
}

int
write_central_body(FILE *stream, const central_body *bc) {
  int nchar;

  nchar = fprintf(stream, "1.0 0.0 %g %g %g 0.0 0.0 0.0 0.0 0.0 0.0 %g %g %g\n",
                  bc->Qp, bc->I, bc->R, bc->spin[0], bc->spin[1], bc->spin[2]);

  if (nchar > 0) {
    /* Any characters written == success, I guess. */
    return 1;
  } else {
    return 0;
  }
}

int write_body_elements(FILE *stream, const body *b) {
  int nchar;
  double e, I, Omega, omega;

  elements_from_body(b, &e, &I, &Omega, &omega);

  nchar = fprintf(stream, "%g %g %g %g %g %g %g %g %g %g %g %g\n",
                  b->m, b->a, e, I, Omega, omega, b->Qp, b->I, b->R, b->spin[0], b->spin[1], b->spin[2]);

  if (nchar > 0) { /* Assume that any characters written == success. */
    return 1;
  } else {
    return 0;
  }
}

int write_central_body_elements(FILE *stream, const central_body *bc) {
  int nchar;

  nchar = fprintf(stream, "1.0 0.0 0.0 0.0 0.0 0.0 %g %g %g %g %g %g\n",
                  bc->Qp, bc->I, bc->R, bc->spin[0], bc->spin[1], bc->spin[2]);

  if (nchar > 0) {
    return 1;
  } else {
    return 0;
  }
}

int
write_body_bin(FILE *stream, const body *b) {
  size_t nwrite;

  nwrite = fwrite(b, sizeof(body), 1, stream);

  return nwrite;
}

int
write_central_body_bin(FILE *stream, const central_body *bc) {
  int nwrite;

  nwrite = fwrite(bc, sizeof(central_body), 1, stream);

  return nwrite;
}
