#include "JointLimitEnergy.h"

namespace PHYSICSMOTION {
//JointLimitEnergy
template <typename T>
JointLimitEnergy<T>::JointLimitEnergy(const ArticulatedBody& body,const Vec& controlPoints,const ThetaTrajectory<T>& tt,T coef)
  :TrajectorySIPEnergy<T>(body,controlPoints,tt,coef) {
  for(int i=0; i<body.nrDOF(); i++)
    _theta.push_back(tt.getControlPointCoeffAll(i));
  _lower=_body.lowerLimit().template cast<T>();
  _upper=_body.upperLimit().template cast<T>();
}
template <typename T>
bool JointLimitEnergy<T>::eval(EFunc* E,GFunc* G,HFunc* H) {
  Vec g=Vec::Zero(_controlPoints.size());
  MatT h=MatT::Zero(_controlPoints.size(),_controlPoints.size());
  for(int i=0; i<(int)_theta.size(); i++) {
    Vec cp=_theta[i]*_controlPoints;
    Vec cpG=Vec::Zero(_theta[i].rows());
    Vec cpH=Vec::Zero(_theta[i].rows());
    for(int j=0; j<cp.size(); j++) {
      //T e=-_coef*(log(cp[j]-_lower[i])+log(_upper[i]-cp[j]));
      T eLowerPart=0;
      T eUpperPart=0;
      T gLowerPart=0;
      T gUpperPart=0;
      T hLowerPart=0;
      T hUpperPart=0;
      if(isfinite(_lower[i]))
        eLowerPart=log(cp[j]-_lower[i]);
      if(isfinite(_upper[i]))
        eUpperPart=log(_upper[i]-cp[j]);
      T e=-_coef*(eLowerPart+eUpperPart);
      if(!isfinite(e)) {
        std::cout<<"cp["<<j<<"]="<<cp[j]<<std::endl;
        std::cout<<"_lower["<<i<<"]="<<_lower[i]<<std::endl;
        std::cout<<"_upper["<<i<<"]="<<_upper[i]<<std::endl;
        std::cout<<"eLowerPart="<<eLowerPart<<std::endl;
        std::cout<<"eUpperPart="<<eUpperPart<<std::endl;
        std::cout<<"infinite joint limit energy"<<std::endl;
        return false;
      }
      if(E)
        (*E)(e);
      if(G) {
        if(isfinite(_lower[i]))
          gLowerPart=1/(cp[j]-_lower[i]);
        if(isfinite(_upper[i]))
          gUpperPart=1/(_upper[i]-cp[j]);
        //cpG[j]=-_coef*(1/(cp[j]-_lower[i])-1/(_upper[i]-cp[j]));
        cpG[j]=-_coef*(gLowerPart-gUpperPart);
      }
      if(H) {
        if(isfinite(_lower[i]))
          hLowerPart=1/(cp[j]-_lower[i])/(cp[j]-_lower[i]);
        if(isfinite(_upper[i]))
          hUpperPart=1/(_upper[i]-cp[j])/(_upper[i]-cp[j]);
        cpH[j]=_coef*(hLowerPart+hUpperPart);
        //cpH[j]=_coef*(1/(cp[j]-_lower[i])/(cp[j]-_lower[i])+1/(_upper[i]-cp[j])/(_upper[i]-cp[j]));
      }
    }
    if(G)
      g+=_theta[i].transpose()*cpG;
    if(H)
      h+=_theta[i].transpose()*cpH.asDiagonal()*_theta[i];
  }
  if(G)
    (*G)(0,g);
  if(H)
    (*H)(0,0,h);
  return true;
}
//instances
template class JointLimitEnergy<FLOAT>;
}
