#include"../src/rings.h"
#include"checks.h"
#include<stdlib.h>

int
main() {
  double x[3], y[3], xnew[3], xmxnew[3];
  double xpy[3], xoy[3];
  int i;

  sranddev();
  
  for (i = 0; i < 3; i++) {
    x[i] = rand() / (RAND_MAX + 1.0);
    y[i] = rand() / (RAND_MAX + 1.0);
  }

  project(x,y,xpy);
  orthogonal_project(x,y,xoy);

  if (!check_close(1e-8, 1e-8, dot(xpy, y) / norm(xpy) / norm(y), 1.0)) return 1;
  if (!(dot(xpy, x) > 0.0)) return 2;
  if (!check_close(1e-8, 1e-8, dot(xoy, y), 0.0)) return 3;
  if (!(dot(xoy, x) > 0.0)) return 4;

  vadd(xpy, xoy, xnew);
  vsub(xnew, x, xmxnew);

  if (!check_close(1e-8, 1e-8, norm(xmxnew), 0.0)) return 5;

  return 0;
}
