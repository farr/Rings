#include"rings.h"

#include<getopt.h>
#include<stdio.h>
#include<stdlib.h>

#include<gsl/gsl_odeiv.h>
#include<gsl/gsl_integration.h>

static double eps = 0.0;
static double T = 1.0;
static double acc = 1e-10;
static double acc_int = 1e-11;
static double h = 0.01;

static const char *short_opts = "e:T:a:h:";

/* Avoid 65--122 (A-Z, a-z, some others) */
enum LONG_ARGS {
  ACC_INT = 0
};

static struct option long_opts[] = {
  { "epsilon", required_argument, NULL, 'e' },
  { "evolve-time", required_argument, NULL, 'T' },
  { "accuracy", required_argument, NULL, 'a' },
  { "accuracy-integral", required_argument, NULL, ACC_INT },
  { "timestep", required_argument, NULL, 'h' },
  { NULL, 0, NULL, 0 }
};

int main (int argc, char *argv[]) {
  char ch;
  gsl_odeiv_control *con;
  gsl_odeiv_step *step;
  gsl_odeiv_evolve *e;
  body *bs;
  double *y;
  int bs_size;
  int bs_insert_index;
  double t = 0.0;
  
  while ((ch = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
    switch(ch) {
    case 'e':
      eps = atof(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    case 'a':
      acc = atof(optarg);
      break;
    case 'h':
      h = atof(optarg);
      break;
    case ACC_INT: 
      acc_int = atof(optarg);
      break;
    default:
      exit(1);
      break;
    }
  }

  con = gsl_odeiv_control_y_new(acc, 0.0); /* We don't really care about anything but absolute error. */

  bs = malloc(sizeof(body));
  bs_size = 1;
  bs_insert_index = 0;

  while (!feof(stdin)) {
    if (bs_insert_index == bs_size) {
      bs_size *= 2;
      bs = realloc(bs, bs_size*sizeof(body));
    }

    read_body(stdin, bs+bs_insert_index);
    bs_insert_index++;
  }
  bs = realloc(bs, bs_insert_index*sizeof(body)); /* Cut down to size. */
  bs_size = bs_insert_index;

  e = gsl_odeiv_evolve_alloc(bs_size*BODY_VECTOR_SIZE);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, bs_size*BODY_VECTOR_SIZE);

  y = malloc(bs_size*BODY_VECTOR_SIZE*sizeof(double));

  while (t != T) {
    int i;

    evolve_system(e, con, step, 
                  &t, T, &h, bs, y, bs_size,
                  acc_int,
                  eps);

    for (i = 0; i < bs_size; i++) {
      fprintf(stdout, "%8g ", t);
      write_body(stdout, bs+i);
    }
  }

  free(bs);
  free(y);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_step_free(step);
  gsl_odeiv_control_free(con);
  return 0;
}
