#include "Articulated/PBDArticulatedGradientInfo.h"
#include "GravityEnergy.h"

namespace PHYSICSMOTION {
template <typename T>
GravityEnergy<T>::GravityEnergy(const ArticulatedBody& body,
                                const Vec& controlPoints,const ThetaTrajectory<T>& tt,
                                int JID,const Vec3T& dir,T coef,bool JTJApprox)
  : TrajectorySIPEnergy<T>(body,controlPoints,tt,coef)
  ,_JTJApporox(JTJApprox)
  ,_JID(JID) {
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body.joint(_JID)._mesh);
  if(!mesh) {
    _localPos=_localDir=Vec3T::Zero();
    _coef=0;
    return;
  }
  _localPos=_body.joint(_JID).getC().template cast<T>();
  ASSERT_MSGV(0<=JID && JID<body.nrJ(),"TargetEnergy._JID=%d out of range [0,%d)!",JID,body.nrJ())
  _localDir=dir;
}
template <typename T>
bool GravityEnergy<T>::eval(EFunc* E,GFunc* G,HFunc* H) {
  Vec GTheta;
  MatT HTheta;
  Mat3XT DTGAll=Mat3XT::Zero(3,_body.nrJ()*4);
  Vec theta=_thetaTrajectory.getPoint(_controlPoints,_thetaTrajectory.getNumSegment());
  PBDArticulatedGradientInfo<T> info(_body,theta);
  Vec3T endEffector=ROTI(info._TM,_JID)*_localPos+CTRI(info._TM,_JID);

  if(E)
    (*E)(-_localDir.dot(endEffector)*_coef);
  if(G) {
    GTheta.setZero(_body.nrDOF());
    TRANSI(DTGAll,_JID)=computeDTG(-_localDir*_coef,_localPos);
    info.DTG(_JID,_body,TRANSI(DTGAll,_JID),[&](int off,T val) {
      SIPEnergy<T>::parallelAdd(GTheta[off],val);
    });
  }
  if(H) {
    HTheta.setZero(_body.nrDOF(),_body.nrDOF());
    if(!_JTJApporox)
      info.toolB(_JID,_body,mapM(DTGAll),[&](int offr,int offc,T val) {
      SIPEnergy<T>::parallelAdd(HTheta(offr,offc),val);
    });
  }
  _thetaTrajectory.assembleEnergy(_thetaTrajectory.getNumSegment(),G,H,&GTheta,&HTheta);
  return true;
}
template <typename T> const typename GravityEnergy<T>::Vec3T& GravityEnergy<T>::localDir() const {
  return _localDir;
}
template <typename T>
int GravityEnergy<T>::JID() const {
  return _JID;
}
//instance
template class GravityEnergy<FLOAT>;
}
