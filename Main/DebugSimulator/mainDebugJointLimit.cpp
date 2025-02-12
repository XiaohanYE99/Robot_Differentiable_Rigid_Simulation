#include "Simulator/JointLimit.h"
#include "ConvexHull/Barrier.h"

using namespace PHYSICSMOTION;

int main() {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  mpfr_float::default_precision(256);
  CLogx p;
  p._x0=10;
  JointLimit::debugSimple<CLogx>(3,p);
  JointLimit::debugBall<CLogx>(3,p);
  JointLimit::debugSwing();
  JointLimit::debugExp();
  JointLimit::debugSimple();
  JointLimit::debugBall();
  return 0;
}
