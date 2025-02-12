#ifndef CUSTOM_PBD_ENERGY_H
#define CUSTOM_PBD_ENERGY_H

#include "CollisionGradInfo.h"

namespace PHYSICSMOTION {
template <typename T>
struct CustomPBDEnergy {
  DECL_MAT_VEC_MAP_TYPES_T
  typedef CollisionGradInfo<T> GradInfo;
  typedef std::function<void(const GradInfo&,const Vec&,MatT&,Vec&)> PBDEnergyForward;
  typedef std::function<void(const GradInfo&,const Vec&,MatT&,MatT&,MatT&)> PBDEnergyBackward;
  typedef std::function<std::tuple<MatT,Vec>(const GradInfo&,const Vec&)> PBDEnergyForwardPython;
  typedef std::function<std::tuple<MatT&,MatT&,MatT&>(const GradInfo&,const Vec&)> PBDEnergyBackwardPython;
  //The additional term is:
  //M*(grad._x-2*pos._x+lastPos._x)/dt^2+C(pos._x,(pos._x-lastPos._x)/dt)
  //The integral of this term with respect to grad._x is:
  //||grad._x-2*pos._x+lastPos._x||_{M/dt^2}^2/2 + grad._x^T*C(pos._x,(pos._x-lastPos._x)/dt)
  CustomPBDEnergy();
  CustomPBDEnergy(int nX,bool forwardOnly=false); //this version is for debug only
  T energy(GradInfo& grad,const GradInfo& pos,const GradInfo& lastPos,T dt,Vec* DE);
  void backward(GradInfo& grad,const GradInfo& pos,const GradInfo& lastPos,T dt);
  void setForward(PBDEnergyForwardPython python);
  void setBackward(PBDEnergyBackwardPython python);
 private:
  //for debug only
  std::vector<MatT> _dMdq;
  MatT _M,_dCdq,_dCddq;
  bool _forwardOnly;
  //void forward(const GradInfo& q,const Vec& dq,MatT& M,Vec& C);
  //Input:
  //    configuration (q)
  //    velocity (dq)
  //Output:
  //    mass matrix (M)
  //    centrifugal and Coriolis force (C)
  PBDEnergyForward _forward;
  PBDEnergyForwardPython _forwardPython;
  //void backward(const GradInfo& q,const Vec& RHS,MatT& dMdqRHS,MatT& dCdq,MatT& dCddq);
  //Input:
  //    configuration (q)
  //    right-hand-side vector (RHS)
  //Output:
  //    derivative of mass matrix multiplied by a given RHS vector (dM/dq * RHS)
  //    derivative of centrifugal and Coriolis force with respect to configuration (dC/dq)
  //    derivative of centrifugal and Coriolis force with respect to velocity (dC/ddq)
  PBDEnergyBackward _backward;
  PBDEnergyBackwardPython _backwardPython;
};
}
#endif
