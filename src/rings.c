#include"rings.h"

#include<getopt.h>
#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_errno.h>

typedef enum {
  OPT_HELP,
  OPT_T,
  OPT_HINI,
  OPT_EPS,
  OPT_ACC,
  OPT_INP,
  OPT_OUT
} opt_flag;

typedef struct {
  double T;
  double hini;
  double eps;
  double epsabs;
  FILE *inp;
  FILE *out;
} configuration;

static struct option opts[] = {
  {"help", no_argument, 0, OPT_HELP},
  {"time", required_argument, 0, OPT_T},
  {"h", required_argument, 0, OPT_HINI},
  {"epsilon", required_argument, 0, OPT_EPS},
  {"accuracy", required_argument, 0, OPT_ACC},
  {"input", required_argument, 0, OPT_INP},
  {"output", required_argument, 0, OPT_OUT},
  {0}
};

static char usage[] = "rings OPTION ...\n"
"\n"
"OPTION is one of:\n"
"\n"
"--help          Display this help message\n"
"--time T        Total evolution time\n"
"--h H           Initial step size\n"
"--epsilon EPS   Softening scale\n"
"--accuracy ACC  Absolute accuracy of evolution\n"
"--input FILE    Input file\n"
"--output FILE   Output file\n";

static int parse_args(int argc, char **argv, struct option *opts, configuration *conf) {
  int flag, index;
  while ((flag = getopt_long(argc, argv, "", opts, &index))  != -1) {
    switch (flag) {
    case OPT_HELP:
      fprintf(stderr, "%s", usage);
      exit(0);
      break;
    case OPT_T:
      conf->T = atof(optarg);
      break;
    case OPT_HINI:
      conf->hini = atof(optarg);
      break;
    case OPT_EPS:
      conf->eps = atof(optarg);
      break;
    case OPT_ACC:
      conf->epsabs = atof(optarg);
      break;
    case OPT_INP:
      conf->inp = fopen(optarg, "r");
      if (conf->inp == 0) {
        fprintf(stderr, "Input file not found: %s\n", optarg);
        exit(1);
      }
      break;
    case OPT_OUT:
      conf->out = fopen(optarg, "w");
      if (conf->out == 0) {
        fprintf(stderr, "Could not open output file: %s\n", optarg);
        exit(1);
      }
      break;
    default:
      fprintf(stderr, "Unrecognized option!\n");
      exit(1);
      break;
    }
  }

  if (optind < argc) {
    fprintf(stderr, "WARNING: unknown options encountered on command line.\n");
  }

  return 0;
}

static body *read_input(FILE *inp, body *bs, int *bsize) {
  int i = 0;
  
  while (!feof(inp)) {
    if (i == *bsize) {
      bs = realloc(bs, 2*i*sizeof(body));
      *bsize = 2*i;
    }

    if (read_body_from_elements(inp, bs+i)) {
      fprintf(stderr, "Error reading body from input (%d bodies read)\n", i);
      exit(1);
    }

    i++;
  }

  if (i==0) {
    fprintf(stderr, "No bodies read from input.\n");
    exit(1);
  }

  bs=realloc(bs, i*sizeof(body));
  *bsize=i;

  return bs;
}

int main(int argc, char **argv) {
  int status;
  configuration conf = {1e9, 1.0, 0.0, 1e-8, stdin, stdout};
  body *bs = malloc(sizeof(body));
  int bsize = 1;
  double t = 0.0;
  gsl_odeiv_evolve *e;
  gsl_odeiv_control *con;
  gsl_odeiv_step *step;
  double *ys;
  int odesize;
  int i;
  double h;

  status = parse_args(argc, argv, opts, &conf);
  if (status != 0) {
    fprintf(stderr, "Error in option parsing.\n");
    exit(1);
  }

  h = conf.hini;

  bs = read_input(conf.inp, bs, &bsize);

  odesize = bsize*BODY_VECTOR_SIZE;
  ys = malloc(odesize*sizeof(double));
  e = gsl_odeiv_evolve_alloc(odesize);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, odesize);
  con = gsl_odeiv_control_secular_new(conf.epsabs);

  do {
    int status;

    for (i = 0; i < bsize; i++) {
      fprintf(conf.out, "%g ", t);
      write_body_elements(conf.out, bs+i);
    }

    status = evolve_system(e, con, step, &t, conf.T, &h, bs, ys, bsize, conf.epsabs, conf.eps);

    if (status != GSL_SUCCESS) {
      fprintf(stderr, "Error in evolution: %d (%s) at %s, line %d\n", status, gsl_strerror(status),
              __FILE__, __LINE__);
      exit(1);
    }

  } while (t != conf.T);

  for (i = 0; i < bsize; i++) {
    fprintf(conf.out, "%g ", t);
    write_body_elements(conf.out, bs+i);
  }

  if (conf.inp != stdin) fclose(conf.inp);
  if (conf.out != stdout) fclose(conf.out);
  free(bs);
  free(ys);
  gsl_odeiv_step_free(step);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(con);

  return 0;
}
