#include"rings.h"

#include<getopt.h>

#include<stdio.h>

int 
main() {
  int read_count = 6;

  while (read_count == 6) {
    double m,a,e,I,O,o;
    body b;

    read_count = fscanf(stdin, " %lg %lg %lg %lg %lg %lg ", 
                        &m, &a, &e, &I, &O, &o);

    if (read_count == 6) {
      init_body_from_elements(&b, m, a, e, I, O, o);
      
      write_body(stdout, &b);
    }
  }

  return 0;
}
