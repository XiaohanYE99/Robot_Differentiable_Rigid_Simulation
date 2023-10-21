#include "TargetEnergy.h"
#include "Articulated/PBDArticulatedGradientInfo.h"

namespace PHYSICSMOTION {
template <typename T>
TargetEnergy<T>::TargetEnergy(const ArticulatedBody& body,
                              const Vec& controlPoints, const ThetaTrajectory<T>& tt,
                              int JID,const Vec3T& tar,const Vec3T& dir,const Vec3T& ori,
                              T coefP, T coefD, bool JTJApprox)
  :TrajectorySIPEnergy<T>(body,controlPoints,tt,coefP),_JTJApprox(JTJApprox),_JID(JID),_coefD(coefD) {
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body.joint(_JID)._mesh);
  if(!mesh) {
    _localDir=_globalDir=_localPos=_globalPos=Vec3T::Zero();
    _coef=_coefD=0;
    return;
  }
  ASSERT_MSGV(0<=JID && JID<body.nrJ(),"TargetEnergy._JID=%d out of range [0,%d)!",JID,body.nrJ())
  //assemble local mesh vertices
  Mat3XT vertices;
  vertices=Mat3XT::Zero(3,mesh->vss().size());
  for(int v=0; v<(int)mesh->vss().size(); ++v)
    vertices.col(v)=mesh->vss()[v].template cast<T>();
  //rest pose
  Vec theta=Vec::Zero(_body.nrDOF());
  PBDArticulatedGradientInfo<T> info(_body,theta);
  //assemble into localDir
  _localDir=dir/dir.norm();
  _localDir=ROTI(info._TM,_JID).transpose()*_localDir;
  _globalDir=tar+ori/ori.norm();
  //assemble into localPos
  Eigen::Index maxId;
  Vec(vertices.transpose()*_localDir).maxCoeff(&maxId);
  _localPos=vertices.col(maxId).dot(_localDir)*_localDir;
  Vec3T mean=vertices.rowwise().mean();
  mean-=mean.dot(_localDir)*_localDir;
  _localPos+=mean;
  _localDir+=_localPos;
  _globalPos=tar;
}
template <typename T>
TargetEnergy<T>::TargetEnergy(const TargetEnergy& ref,
                              const Vec& controlPoints,const ThetaTrajectory<T>& tt)
  :TrajectorySIPEnergy<T>(ref._body,controlPoints,tt,ref._coef),
   _localPos(ref._localPos),_globalPos(ref._globalPos),
   _localDir(ref._localDir),_globalDir(ref._globalDir),
   _JTJApprox(ref._JTJApprox),_JID(ref._JID),_coefD(ref._coefD) {}
template <typename T>
bool TargetEnergy<T>::eval(EFunc* E,GFunc* G,HFunc* H) {
  if(_coef>0)
    eval(_localPos,_globalPos,_coef,E,G,H);
  if(_coefD>0)
    eval(_localDir,_globalDir,_coefD,E,G,H);
  return true;
}
template <typename T>
const typename TargetEnergy<T>::Vec3T& TargetEnergy<T>::globalPos() const {
  return _globalPos;
}
template <typename T>
const typename TargetEnergy<T>::Vec3T& TargetEnergy<T>::localPos() const {
  return _localPos;
}
template <typename T>
const typename TargetEnergy<T>::Vec3T&TargetEnergy<T>::localDir() const {
  return _localDir;
}
template <typename T>
int TargetEnergy<T>::JID() const {
  return _JID;
}
template <typename T>
void TargetEnergy<T>::eval(const Vec3T& local,const Vec3T& global,T coef,EFunc* E,GFunc* G,HFunc* H) const {
  Vec GTheta;
  MatT HTheta;
  Mat3XT DTGAll=Mat3XT::Zero(3,_body.nrJ()*4);
  Vec theta=_thetaTrajectory.getPoint(_controlPoints,_thetaTrajectory.getNumSegment());
  PBDArticulatedGradientInfo<T> info(_body,theta);
  Vec3T endEffector=ROTI(info._TM,_JID)*local+CTRI(info._TM,_JID);
  if(E)
    (*E)((endEffector-global).squaredNorm()*coef/2.);
  if(G) {
    GTheta.setZero(_body.nrDOF());
    TRANSI(DTGAll,_JID)=computeDTG((endEffector-global)*coef,local);
    info.DTG(_JID,_body,TRANSI(DTGAll,_JID),[&](int off,T val) {
      SIPEnergy<T>::parallelAdd(GTheta[off],val);
    });
  }
  if(H) {
    HTheta.setZero(_body.nrDOF(),_body.nrDOF());
    if(!_JTJApprox)
      info.toolB(_JID,_body,mapM(DTGAll),[&](int offr,int offc,T val) {
      SIPEnergy<T>::parallelAdd(HTheta(offr,offc),val);
    });
    info.toolALR(_JID,_JID,_body,Mat3T::Identity()*coef,local,local,[&](int offr,int offc,T val) {
      SIPEnergy<T>::parallelAdd(HTheta(offr,offc),val);
    });
  }
  _thetaTrajectory.assembleEnergy(_thetaTrajectory.getNumSegment(),G,H,&GTheta,&HTheta);
}
//instance
template class TargetEnergy<FLOAT>;
}
