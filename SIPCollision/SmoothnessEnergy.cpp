#include "SmoothnessEnergy.h"

namespace PHYSICSMOTION {
template <typename T>
SmoothnessEnergy<T>::SmoothnessEnergy(const ArticulatedBody& body,const Vec& controlPoints,const ThetaTrajectory<T>& tt,T coef,int d)
  :TrajectorySIPEnergy<T>(body,controlPoints,tt,coef) {
  T e;
  Vec g;
  MatT h;
  g.setZero(_controlPoints.size());
  h.setZero(_controlPoints.size(),_controlPoints.size());
  _thetaTrajectory.smoothnessRegularizer(_controlPoints,&e,&g,&h,d);
  _hessian=h.sparseView();
}
template <typename T>
bool SmoothnessEnergy<T>::eval(EFunc* E,GFunc* G,HFunc* H) {
  Vec g=_hessian*_controlPoints;
  if(E)
    (*E)(_controlPoints.dot(g)*_coef/2);
  if(G)
    (*G)(0,g*_coef);
  if(H)
    (*H)(0,0,_hessian*_coef);
  return true;
}
//instances
template class SmoothnessEnergy<FLOAT>;
}
