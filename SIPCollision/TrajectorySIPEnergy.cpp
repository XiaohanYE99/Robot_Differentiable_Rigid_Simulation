#include "TrajectorySIPEnergy.h"
#include <Utils/DebugGradient.h>
#include <random>

namespace PHYSICSMOTION {
//SIPEnergy
template <typename T>
bool SIPEnergy<T>::eval(T* E,Vec* G,MatT* H) {
  EFunc EF=[&](T val)->void {
    parallelAdd(*E,val);
  };
  GFunc GF=[&](int off,const Vec& val)->void {
    parallelAdd<-1>(*G,off,val);
  };
  HFunc HF=[&](int offr,int offc,const MatT& val)->void {
    parallelAdd<-1,-1>(*H,offr,offc,val);
  };
  return eval(E?&EF:NULL,G?&GF:NULL,H?&HF:NULL);
}
//TrajectorySIPEnergy
template <typename T>
TrajectorySIPEnergy<T>::TrajectorySIPEnergy(const ArticulatedBody& body,const Vec& controlPoints,const ThetaTrajectory<T>& tt,T coef)
  :_body(body),_controlPoints(controlPoints),_thetaTrajectory(tt),_coef(coef) {}
template <typename T>
void TrajectorySIPEnergy<T>::resetCoef(T newCoef) {
  _coef=newCoef;
}
template <typename T>
bool TrajectorySIPEnergy<T>::debug(T perturbRange) {
  T E0=0,E,E2;
  MatT H,H0;
  Vec G0,G,G2,dx,x0;
  DEFINE_NUMERIC_DELTA_T(T)
  E0=rand()/(T)RAND_MAX*2-1;
  E=E2=E0;
  G0.setRandom(_controlPoints.size());
  G=G2=G0;
  H0.setRandom(_controlPoints.size(),_controlPoints.size());
  H=H0;

  x0=_controlPoints;
  dx=Vec::Random(_controlPoints.size());
  const_cast<Vec&>(_controlPoints)+=dx*perturbRange;
  if(!SIPEnergy<T>::eval(&E,&G,&H))
    return false;
  const_cast<Vec&>(_controlPoints)+=dx*DELTA;
  if(!SIPEnergy<T>::eval(&E2,&G2,(MatT*)NULL))
    return false;
  //std::cout << "E= " << std::setprecision(20) << E-E0 << std::endl;
  //std::cout << "E1=" << std::setprecision(20) << E2-E0 << std::endl;
  //std::cout << "G= " << std::setprecision(20) << (G-G0).dot(dx) << std::endl;
  //std::cout << "G1=" << std::setprecision(20) << (E2-E)/DELTA << std::endl;
  //std::cout << "H= " << std::setprecision(20) << ((H-H0)*dx).norm() << std::endl;
  //std::cout << "H1=" << std::setprecision(20) << ((G2-G)/DELTA).norm() << std::endl;
  DEBUG_GRADIENT("dE",(G-G0).dot(dx),(G-G0).dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("dG",((H-H0)*dx).norm(),((H-H0)*dx-(G2-G)/DELTA).norm())
  const_cast<Vec&>(_controlPoints)=x0;
  return true;
}
//instance
template class TrajectorySIPEnergy<FLOAT>;
}
