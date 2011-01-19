#include"rings.h"

#include<stdio.h>

int
read_body(FILE *stream, body *b) {
  int nscan;

  nscan = fscanf(stream, " %lg %lg %lg %lg %lg %lg %lg %lg ",
                 &(b->m), &(b->a),
                 b->L, b->L+1, b->L+2,
                 b->A, b->A+1, b->A+2);

  if (nscan == 8) {
    return 1;
  } else {
    return 0;
  }
}

int
read_body_from_elements(FILE *stream, body *b) {
  double m, a, e, I, Omega, omega;
  int nscan;

  nscan = fscanf(stream, " %lg %lg %lg %lg %lg %lg ",
                 &m, &a,
                 &e, &I,
                 &Omega, &omega);

  if (nscan != 6) return 1;

  init_body_from_elements(b, m, a, e, I, Omega, omega);

  return 0;
}

int
read_body_bin(FILE *stream, body *b) {
  size_t nread;

  nread = fread(b, sizeof(body), 1, stream);

  return nread;
}

int
write_body(FILE *stream, const body *b) {
  int nchar;
  
  nchar = fprintf(stream, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n",
                  b->m, b->a, b->L[0], b->L[1], b->L[2],
                  b->A[0], b->A[1], b->A[2]);

  if (nchar > 0) { /* Assume that any written characters == total success! */
    return 0;
  } else {
    return 1;
  }
}

int write_body_elements(FILE *stream, const body *b) {
  int nchar;
  double e, I, Omega, omega;

  elements_from_body(b, &e, &I, &Omega, &omega);

  nchar = fprintf(stream, "%g %g %g %g %g %g\n",
                  b->m, b->a, e, I, Omega, omega);

  if (nchar > 0) { /* Assume that any characters written == success. */
    return 0;
  } else {
    return 1;
  }
}

int
write_body_bin(FILE *stream, const body *b) {
  size_t nwrite;

  nwrite = fwrite(b, sizeof(body), 1, stream);

  return nwrite;
}
