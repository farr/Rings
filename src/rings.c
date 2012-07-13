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
  OPT_EPSQUAD,
  OPT_EPSINT,
  OPT_THIN,
  OPT_INP,
  OPT_OUT,
  OPT_TIDES
} opt_flag;

typedef struct {
  double T;
  double hini;
  double eps;
  double epsquad;
  double epsint;
  int thin_factor;
  FILE *inp;
  FILE *out;
  int tides_p;
} configuration;

static struct option opts[] = {
  {"help", no_argument, 0, OPT_HELP},
  {"time", required_argument, 0, OPT_T},
  {"h", required_argument, 0, OPT_HINI},
  {"eps", required_argument, 0, OPT_EPS},
  {"epsquad", required_argument, 0, OPT_EPSQUAD},
  {"epsint", required_argument, 0, OPT_EPSINT},
  {"thin", required_argument, 0, OPT_THIN},
  {"input", required_argument, 0, OPT_INP},
  {"output", required_argument, 0, OPT_OUT},
  {"tides", no_argument, 0, OPT_TIDES},
  {0}
};

static char usage[] = "rings OPTION ...\n"
"\n"
"OPTION is one of:\n"
"\n"
"--help          Display this help message\n"
"--time T        Total evolution time\n"
"--h H           Initial step size\n"
"--eps EPS       Softening scale\n"
"--epsquad EPS   Quadrature relative accuracy\n"
"--epsint EPS    Integration relative accuracy\n"
"--thin NTHIN    Output only every NTHIN steps\n"
"--input FILE    Input file\n"
"--output FILE   Output file\n"
"--tides         Apply tidal evolution\n";

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
    case OPT_EPSQUAD:
      conf->epsquad = atof(optarg);
      break;
    case OPT_EPSINT:
      conf->epsint = atof(optarg);
      break;
    case OPT_THIN:
      conf->thin_factor = atoi(optarg);
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
    case OPT_TIDES:
      conf->tides_p = 1;
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

static body *read_input(FILE *inp, central_body *bc, body *bs, int *bsize) {
  int i = 0;

  if (read_central_body(inp, bc) != 1) {
    fprintf(stderr, "Error reading central body from input\n");
    exit(1);
  }
  
  while (!feof(inp)) {
    if (i == *bsize) {
      bs = realloc(bs, 2*i*sizeof(body));
      *bsize = 2*i;
    }

    if (read_body_from_elements(inp, bs+i) != 1) {
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

static int write_bodies(FILE *stream, const double t, const central_body *bc, const body *bs, const size_t nbs) {
  size_t i;

  fprintf(stream, "%.1f ", t);
  write_central_body_elements(stream, bc);

  for (i = 0; i < nbs; i++) {
    fprintf(stream, "%.1f ", t);
    write_body_elements(stream, &(bs[i]));
  }

  fflush(stream);

  return nbs + 1;
}

int main(int argc, char **argv) {
  int status;
  configuration conf = {1e9, 1.0, 0.0, 1e-10, 1e-8, 1, stdin, stdout, 0};
  body *bs = malloc(sizeof(body));
  central_body bc;
  int bsize = 1;
  double t = 0.0;
  gsl_odeiv_evolve *e;
  gsl_odeiv_control *con;
  gsl_odeiv_step *step;
  gsl_integration_workspace *ws;
  const size_t nws = 100000;
  double *ys;
  int odesize;
  int i;
  double h;
  const int NRETRIES_MAX = 5;
  int nstep = 0;

  status = parse_args(argc, argv, opts, &conf);
  if (status != 0) {
    fprintf(stderr, "Error in option parsing.\n");
    exit(1);
  }

  h = conf.hini;

  bs = read_input(conf.inp, &bc, bs, &bsize);

  odesize = body_size_to_vector_size(bsize);
  ys = malloc(odesize*sizeof(double));
  e = gsl_odeiv_evolve_alloc(odesize);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, odesize);
  con = gsl_odeiv_control_secular_new(conf.epsint);
  ws = gsl_integration_workspace_alloc(nws);

  do {
    int status;
    int nretries = 0;
    double told = t;

    if (nstep % conf.thin_factor == 0) {
      write_bodies(conf.out, t, &bc, bs, bsize);
    }

    do {
      status = evolve_system(e, con, step, &t, conf.T, &h, &bc, bs, ys, bsize, conf.epsquad, conf.eps);

      if (status != GSL_SUCCESS) {
        nretries++;
        h = h / 10.0;
        gsl_odeiv_evolve_reset(e);
        gsl_odeiv_step_reset(step);
        fprintf(stderr, "Retrying step at time %g, hnew = %g\n", t, h);
      }
    } while (status != GSL_SUCCESS && nretries <= NRETRIES_MAX);

    if (status != GSL_SUCCESS) {
      fprintf(stderr, "Error in evolution: %d (%s) at %s, line %d\n", status, gsl_strerror(status),
             __FILE__, __LINE__);
      exit(1);
    }

    if (told == t) {
      fprintf(stderr, "Error in evolution: driven to zero stepsize at time %g (retried %d steps), at %s, line %d\n", 
              t, nretries, __FILE__, __LINE__);
      exit(1);
    }

    nstep++;
  } while (t != conf.T);

  write_bodies(conf.out, t, &bc, bs, bsize);

  if (conf.inp != stdin) fclose(conf.inp);
  if (conf.out != stdout) fclose(conf.out);
  free(bs);
  free(ys);
  gsl_odeiv_step_free(step);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(con);
  gsl_integration_workspace_free(ws);

  return 0;
}
