#include"../src/rings.h"
#include"checks.h"

#include<gsl/gsl_integration.h>

int main() {
  body b1 = {1e-6, 0.1, {0.0, 0.0, sqrt(1.0-0.04)}, {0.2, 0.0, 0.0}};
  body b2 = {1e-6, 2.0, {0.0, 0.0, sqrt(1.0-0.09)}, {0.3, 0.0, 0.0}};
  const size_t ws_size = 10000;
  gsl_integration_workspace *ws1, *ws2;
  double rhs1[BODY_VECTOR_SIZE], rhs2[BODY_VECTOR_SIZE];
  int status = 0;

  ws1 = gsl_integration_workspace_alloc(ws_size);
  ws2 = gsl_integration_workspace_alloc(ws_size);
  assert(ws1 != 0);
  assert(ws2 != 0);

  raw_average_rhs(0.0, &b1, &b2, ws1, ws_size, ws2, ws_size, rhs1);
  raw_average_rhs(0.0, &b2, &b1, ws1, ws_size, ws2, ws_size, rhs2);
  
  /* fprintf(stderr, "%g %g %g %g %g %g %g %g\n",  */
  /*         rhs1[BODY_M_INDEX], rhs1[BODY_A_INDEX], */
  /*         rhs1[BODY_L_INDEX], rhs1[BODY_L_INDEX+1], rhs1[BODY_L_INDEX+2], */
  /*         rhs1[BODY_A_INDEX], rhs1[BODY_A_INDEX+1], rhs1[BODY_A_INDEX+2]); */
  /* fprintf(stderr, "%g %g %g %g %g %g %g %g\n",  */
  /*         rhs2[BODY_M_INDEX], rhs2[BODY_A_INDEX], */
  /*         rhs2[BODY_L_INDEX], rhs2[BODY_L_INDEX+1], rhs2[BODY_L_INDEX+2], */
  /*         rhs2[BODY_A_INDEX], rhs2[BODY_A_INDEX+1], rhs2[BODY_A_INDEX+2]); */

  gsl_integration_workspace_free(ws1);
  gsl_integration_workspace_free(ws2);
  return status;
}
