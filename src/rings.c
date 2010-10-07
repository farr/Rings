#include<getopt.h>
#include<stdio.h>
#include<stdlib.h>

static double eps = 0.0;
static double T = 1.0;

static const char *short_opts = "e:T:";

static struct option long_opts[] = {
  { "epsilon", required_argument, NULL, 'e' },
  { "evolve-time", required_argument, NULL, 'T' },
  { NULL, 0, NULL, 0 }
};

int main (int argc, char *argv[]) {
  char ch;

  while ((ch = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
    switch(ch) {
    case 'e':
      eps = atof(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    default:
      exit(1);
      break;
    }
  }

  fprintf(stderr, "eps = %g, T = %g\n", eps, T);

  return 0;
}
