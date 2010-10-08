#include"rings.h"

#include<stdio.h>

int
main() {
  int read_count = 1;

  while (read_count == 1) {
    double t;
    double e, I, O, o;
    body b;

    read_count = fscanf(stdin, " %lg ", &t);
    if (read_count != 1) break;

    read_count = read_body(stdin, &b);
    if (read_count != 1) break;

    elements_from_body(&b, &e, &I, &O, &o);

    fprintf(stdout, "%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n",
            t, b.m, b.a, e, I, O, o);
  }

  return 0;
}
