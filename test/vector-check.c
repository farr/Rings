#include"../src/rings.h"
#include"checks.h"
#include<math.h>

static const double eps = 1e-8;

static int 
check_dot() {
  double x[] = {1.0, -2.0, 3.0};
  double y[] = {-4.0, 5.0, 6.0};

  return check_close(eps,eps,dot(x,y), 4.0);
}

static int
check_norm() {
  double x[] = {1.0, -2.0, 3.0};

  return check_close(eps,eps,norm(x),sqrt(14.0));
}

static int
check_cross() {
  double x[] = {1.0, -2.0, 3.0};
  double y[] = {-4.0, 5.0, -6.0};
  double z[3];
  double zexact[] = {-3.0, -6.0, -3.0};

  cross(x,y,z);
  
  return check_vector_close(eps,eps,z,zexact);
}

static int
check_add_sub_scale() {
  double x[] = {1.0, -2.0, 3.0};
  double y[] = {-4.0, 5.0, -6.0};
  double xpy[3];
  double xmmy[3];

  vadd(x,y,xpy);
  vscale(-1.0,y,xmmy);
  vsub(x,xmmy,xmmy);

  return check_vector_close(eps,eps,xpy,xmmy);
}

int main() {

  if (!check_dot()) return 1;

  if (!check_norm()) return 2;

  if (!check_cross()) return 3;

  if (!check_add_sub_scale()) return 4;

  return 0;
}
