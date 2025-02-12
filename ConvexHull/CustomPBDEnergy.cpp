#include "CustomPBDEnergy.h"

namespace PHYSICSMOTION {
//CollisionGradInfo
template <typename T>
CustomPBDEnergy<T>::CustomPBDEnergy() {
  _forward=[](const GradInfo& pos,const Vec& dotq,MatT& M,Vec& C) {
    M.setZero(pos._info._xM.size(),pos._info._xM.size());
    C.setZero(pos._info._xM.size());
  };
  _backward=[](const GradInfo& pos,const Vec& RHS,MatT& dMdqRHS,MatT& dCdq,MatT& dCddq) {
    dMdqRHS.setZero(pos._info._xM.size(),pos._info._xM.size());
    dCdq.setZero(pos._info._xM.size(),pos._info._xM.size());
    dCddq.setZero(pos._info._xM.size(),pos._info._xM.size());
  };
}
template <typename T>
CustomPBDEnergy<T>::CustomPBDEnergy(int nX,bool forwardOnly) {
  _forwardOnly=forwardOnly;
  for(int i=0; i<nX; i++) {
    MatT rand=MatT::Random(nX,nX);
    _dMdq.push_back(rand*rand.transpose());
  }
  _dCdq.setRandom(nX,nX);
  _dCddq.setRandom(nX,nX);
  _forward=[&](const GradInfo& pos,const Vec& dotq,MatT& M,Vec& C) {
    std::cout << "Calling forward!" << std::endl;
    M.setZero(pos._info._xM.size(),pos._info._xM.size());
    C=_dCdq*pos._info._xM+_dCddq*dotq;
    for(int i=0; i<(int)_dMdq.size(); i++)
      M+=_dMdq[i]*pos._info._xM[i]*pos._info._xM[i];
  };
  _backward=[&](const GradInfo& pos,const Vec& RHS,MatT& dMdqRHS,MatT& dCdq,MatT& dCddq) {
    dMdqRHS.setZero(pos._info._xM.size(),pos._info._xM.size());
    dCdq.setZero(pos._info._xM.size(),pos._info._xM.size());
    dCddq.setZero(pos._info._xM.size(),pos._info._xM.size());
    if(_forwardOnly) {
      std::cout << "Skipping backward!" << std::endl;
      return;
    } else std::cout << "Calling backward!" << std::endl;
    for(int i=0; i<(int)_dMdq.size(); i++)
      dMdqRHS.col(i)+=_dMdq[i]*RHS*2*pos._info._xM[i];
    dCdq+=_dCdq;
    dCddq+=_dCddq;
  };
}
template <typename T>
T CustomPBDEnergy<T>::energy(GradInfo& grad,const GradInfo& pos,const GradInfo& lastPos,T dt,Vec* DE) {
  //compute terms
  T coef = 1/dt/dt;
  Vec vel=(pos._info._xM-lastPos._info._xM)/dt,C;
  Vec accel=grad._info._xM-2*pos._info._xM+lastPos._info._xM;
  _forward(pos,vel,_M,C);
  ASSERT_MSG(_M.rows()==pos._info._xM.size(),"Incorrect size of M.rows()");
  ASSERT_MSG(_M.cols()==pos._info._xM.size(),"Incorrect size of M.cols()");
  ASSERT_MSG(C.size()==pos._info._xM.size(),"Incorrect size of C.size()");
  //assemble energy
  if(DE)
    *DE += _M*accel*coef+C;
  grad._HTheta += _M*coef;
  return accel.dot(_M*accel)*(coef/2) + grad._info._xM.dot(C);
}
template <typename T>
void CustomPBDEnergy<T>::backward(GradInfo& grad,const GradInfo& pos,const GradInfo& lastPos,T dt) {
  //compute terms
  T coef = 1/dt/dt;
  MatT dMdqRHS,dCdq,dCddq;
  Vec accel=grad._info._xM-2*pos._info._xM+lastPos._info._xM;
  _backward(pos,accel*coef,dMdqRHS,dCdq,dCddq);
  ASSERT_MSG(dMdqRHS.rows()==pos._info._xM.size(),"Incorrect size of dMdqRHS.rows()");
  ASSERT_MSG(dMdqRHS.cols()==pos._info._xM.size(),"Incorrect size of dMdqRHS.cols()");
  ASSERT_MSG(dCdq.rows()==pos._info._xM.size(),"Incorrect size of dCdq.rows()");
  ASSERT_MSG(dCdq.cols()==pos._info._xM.size(),"Incorrect size of dCdq.cols()");
  ASSERT_MSG(dCddq.rows()==pos._info._xM.size(),"Incorrect size of dCddq.rows()");
  ASSERT_MSG(dCddq.cols()==pos._info._xM.size(),"Incorrect size of dCddq.cols()");
  //assemble energy
  grad._HThetaL += -_M*(coef*2) + dMdqRHS + dCdq + dCddq/dt;
  grad._HThetaLL += _M*coef - dCddq/dt;
}
template <typename T>
void CustomPBDEnergy<T>::setForward(PBDEnergyForwardPython python) {
  _forwardPython=python;
  _forward=[&](const GradInfo& pos,const Vec& dotq,MatT& M,Vec& C) {
    auto ret=_forwardPython(pos,dotq);
    M=std::get<0>(ret);
    C=std::get<1>(ret);
  };
}
template <typename T>
void CustomPBDEnergy<T>::setBackward(PBDEnergyBackwardPython python) {
  _backwardPython=python;
  _backward=[&](const GradInfo& pos,const Vec& RHS,MatT& dMdqRHS,MatT& dCdq,MatT& dCddq) {
    auto ret=_backwardPython(pos,RHS);
    dMdqRHS=std::get<0>(ret);
    dCdq=std::get<1>(ret);
    dCddq=std::get<2>(ret);
  };
}
//instance
template struct CustomPBDEnergy<FLOAT>;
}
