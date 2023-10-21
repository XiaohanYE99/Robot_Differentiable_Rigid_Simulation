#include "PBCentroidBodyDynamicsGradientInfo.h"
#include <Environment/DeformedEnvironment.h>
#include <Utils/CrossSpatialUtils.h>
#include <Utils/RotationUtils.h>

namespace PHYSICSMOTION {
Eigen::Matrix<double,3,1> XX(0.1,0.2,0.3);
//PBCentroidBodyDynamicsMapGradientInfoMap
template <typename T>
PBCentroidBodyDynamicsGradientInfoMap<T>::PBCentroidBodyDynamicsGradientInfoMap()
  :_RTM(NULL,3,4,Eigen::OuterStride<>(3)),
   _wM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTM(NULL,3,0,Eigen::OuterStride<>(3)),
   _xM(NULL,3,0,Eigen::OuterStride<>(3)),
   _fM(NULL,3,0,Eigen::OuterStride<>(3)) {}
template <typename T>
PBCentroidBodyDynamicsGradientInfoMap<T>::PBCentroidBodyDynamicsGradientInfoMap
(const PBCentroidBodyDynamicsGradientInfoMap<T>& other)
  :_RTM((Mat3X4TM&)other._RTM),
   _wM((Mat3XTM&)other._wM),
   _DTM((Mat3XTM&)other._DTM),
   _xM((Mat3XTM&)other._xM),
   _fM((Mat3XTM&)other._fM),
   _XTX(other._XTX),
   _g(other._g),
   _dt(other._dt),
   _M(other._M) {}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Vec6T
PBCentroidBodyDynamicsGradientInfoMap<T>::C
(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
 const PBCentroidBodyDynamicsGradientInfoMap<T>& IN,
 std::function<void(int,int,Vec3T)> DCDINN,
 std::function<void(int,int,Vec3T)> DCDIN,
 std::function<void(int,int,Vec3T)> DCDI,
 std::function<void(int,int,Mat3T)> DCDX,
 std::function<void(int,int,Mat3T)> DCDF) const {
  Vec6T ret=Vec6T::Zero();
  Mat3T DI12Dw1,DI12Dw2;
  Mat3T DI02Dw0,DI02Dw2;
  ret.template segment<3>(0)+=(IXX(ROT(IN._RTM),ROT(_RTM),mapM(DCDIN?&DI12Dw1:NULL),mapM(DCDI?&DI12Dw2:NULL))*2-
                               IXX(ROT(INN._RTM),ROT(_RTM),mapM(DCDINN?&DI02Dw0:NULL),mapM(DCDI?&DI02Dw2:NULL)))/_dt/_dt;
  ret.template segment<3>(3)+=(CTR(_RTM)-2*CTR(IN._RTM)+CTR(INN._RTM))*_M/_dt/_dt;
  for(int c=0; c<_xM.cols(); c++) {
    ret.template segment<3>(0)-=(_xM.col(c)-CTR(_RTM)).cross(_fM.col(c));
    ret.template segment<3>(3)-=_fM.col(c);
    if(DCDX)
      DCDX(0,c*3,cross<T>(_fM.col(c)));
    if(DCDF) {
      DCDF(0,c*3,-cross<T>(_xM.col(c)-CTR(_RTM)));
      DCDF(3,c*3,-Mat3T::Identity());
    }
    if(DCDI)
      for(int d=0; d<_DTM.cols(); d++)
        DCDI(0,d,-_fM.col(c).cross(_DTM.col(d)));
  }
  ret.template segment<3>(3)-=_g*_M;
  {
    //rotation
    if(DCDINN)
      for(int d=0; d<INN._wM.cols(); d++)
        DCDINN(0,d,-DI02Dw0*INN._wM.col(d)/_dt/_dt);
    if(DCDIN)
      for(int d=0; d<IN._wM.cols(); d++)
        DCDIN(0,d,DI12Dw1*IN._wM.col(d)*2/_dt/_dt);
    if(DCDI)
      for(int d=0; d<_wM.cols(); d++)
        DCDI(0,d,(DI12Dw2*2-DI02Dw2)*_wM.col(d)/_dt/_dt);
    //translation
    if(DCDINN)
      for(int d=0; d<INN._DTM.cols(); d++)
        DCDINN(3,d,INN._DTM.col(d)*_M/_dt/_dt);
    if(DCDIN)
      for(int d=0; d<IN._DTM.cols(); d++)
        DCDIN(3,d,-IN._DTM.col(d)*_M*2/_dt/_dt);
    if(DCDI)
      for(int d=0; d<_DTM.cols(); d++)
        DCDI(3,d,_DTM.col(d)*_M/_dt/_dt);
  }
  return ret;
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Vec6T
PBCentroidBodyDynamicsGradientInfoMap<T>::C
(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
 const PBCentroidBodyDynamicsGradientInfoMap<T>& IN,
 Mat6TM DCDINN,Mat6TM DCDIN,Mat6TM DCDI,
 Mat6XTM DCDX,Mat6XTM DCDF) const {
  if(DCDINN.data())
    DCDINN.setZero();
  if(DCDIN.data())
    DCDIN.setZero();
  if(DCDI.data())
    DCDI.setZero();
  if(DCDX.data())
    DCDX.setZero();
  if(DCDF.data())
    DCDF.setZero();
  std::function<void(int,int,Vec3T)> DCDINNFunc=[&](int r,int c,Vec3T blk) {
    DCDINN.template block<3,1>(r,c)+=blk;
  };
  std::function<void(int,int,Vec3T)> DCDINFunc=[&](int r,int c,Vec3T blk) {
    DCDIN.template block<3,1>(r,c)+=blk;
  };
  std::function<void(int,int,Vec3T)> DCDIFunc=[&](int r,int c,Vec3T blk) {
    DCDI.template block<3,1>(r,c)+=blk;
  };
  std::function<void(int,int,Mat3T)> DCDXFunc=[&](int r,int c,Mat3T blk) {
    DCDX.template block<3,3>(r,c)+=blk;
  };
  std::function<void(int,int,Mat3T)> DCDFFunc=[&](int r,int c,Mat3T blk) {
    DCDF.template block<3,3>(r,c)+=blk;
  };
  return C(INN,IN,
           DCDINN.data()?DCDINNFunc:NULL,
           DCDIN.data()?DCDINFunc:NULL,
           DCDI.data()?DCDIFunc:NULL,
           DCDX.data()?DCDXFunc:NULL,
           DCDF.data()?DCDFFunc:NULL);
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Vec6T
PBCentroidBodyDynamicsGradientInfoMap<T>::C
(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
 const PBCentroidBodyDynamicsGradientInfoMap<T>& IN) const {
  return C(INN,IN,NULL,NULL,NULL,NULL,NULL);
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfoMap<T>::setContactForce(int index,const Vec3T& pos,const Vec3T& force) {
  _xM.col(index)=pos;
  _fM.col(index)=force;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfoMap<T>::DTG(Mat3X4T GK,std::function<void(int,T)> DTG) const {
  Vec3T w=invCrossMatTrace<T>(ROT(_RTM)*ROT(GK).transpose());
  for(int d=0; d<_wM.cols(); d++)
    DTG(d,w.dot(_wM.col(d)));
  for(int d=0; d<_DTM.cols(); d++)
    DTG(d,CTR(GK).dot(_DTM.col(d)));
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Mat3X4T PBCentroidBodyDynamicsGradientInfoMap<T>::getTrans() const {
  return _RTM;
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Vec6T PBCentroidBodyDynamicsGradientInfoMap<T>::getDOF() const {
  Vec6T ret;
  ret.template segment<3>(0)=CTR(_RTM);
  if(_isEXP)
    ret.template segment<3>(3)=invExpW<T>(ROT(_RTM));
  else ret.template segment<3>(3)=invEulerX1Y3Z2<T>(ROT(_RTM));
  return ret;
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfoMap<T>::Vec3T
PBCentroidBodyDynamicsGradientInfoMap<T>::IXX(const Mat3T& Ra,const Mat3T& Rb,Mat3TM dIdwa,Mat3TM dIdwb) const {
  Vec3T I;
  T x0=Ra(1,0)*_XTX(0,0) + Ra(1,1)*_XTX(1,0) + Ra(1,2)*_XTX(2,0);
  T x1=Rb(2,0)*x0;
  T x2=Ra(1,0)*_XTX(0,1) + Ra(1,1)*_XTX(1,1) + Ra(1,2)*_XTX(2,1);
  T x3=Rb(2,1)*x2;
  T x4=Ra(1,0)*_XTX(0,2) + Ra(1,1)*_XTX(1,2) + Ra(1,2)*_XTX(2,2);
  T x5=Rb(2,2)*x4;
  T x6=x1 + x3 + x5;
  T x7=Ra(2,0)*_XTX(0,0) + Ra(2,1)*_XTX(1,0) + Ra(2,2)*_XTX(2,0);
  T x8=Rb(1,0)*x7;
  T x9=Ra(2,0)*_XTX(0,1) + Ra(2,1)*_XTX(1,1) + Ra(2,2)*_XTX(2,1);
  T x10=Rb(1,1)*x9;
  T x11=Ra(2,0)*_XTX(0,2) + Ra(2,1)*_XTX(1,2) + Ra(2,2)*_XTX(2,2);
  T x12=Rb(1,2)*x11;
  T x13=-x10 - x12 - x8;
  T x14=Rb(0,0)*x7;
  T x15=Rb(0,1)*x9;
  T x16=Rb(0,2)*x11;
  T x17=x14 + x15 + x16;
  T x18=Ra(0,0)*_XTX(0,0) + Ra(0,1)*_XTX(1,0) + Ra(0,2)*_XTX(2,0);
  T x19=Rb(2,0)*x18;
  T x20=Ra(0,0)*_XTX(0,1) + Ra(0,1)*_XTX(1,1) + Ra(0,2)*_XTX(2,1);
  T x21=Rb(2,1)*x20;
  T x22=Ra(0,0)*_XTX(0,2) + Ra(0,1)*_XTX(1,2) + Ra(0,2)*_XTX(2,2);
  T x23=Rb(2,2)*x22;
  T x24=-x19 - x21 - x23;
  T x25=Rb(1,0)*x18;
  T x26=Rb(1,1)*x20;
  T x27=Rb(1,2)*x22;
  T x28=x25 + x26 + x27;
  T x29=Rb(0,0)*x0;
  T x30=Rb(0,1)*x2;
  T x31=Rb(0,2)*x4;
  T x32=-x29 - x30 - x31;
  T x33=Rb(2,0)*x7;
  T x34=Rb(2,1)*x9;
  T x35=Rb(2,2)*x11;
  T x36=-x33 - x34 - x35;
  T x37=Rb(1,0)*x0;
  T x38=Rb(1,1)*x2;
  T x39=Rb(1,2)*x4;
  T x40=-x37 - x38 - x39;
  T x41=Rb(0,0)*x18;
  T x42=Rb(0,1)*x20;
  T x43=Rb(0,2)*x22;
  T x44=-x41 - x42 - x43;
  T x45=x33 + x34 + x35;
  T x46=x37 + x38 + x39;
  T x47=x41 + x42 + x43;
  I[0]=x13 + x6;
  I[1]=x17 + x24;
  I[2]=x28 + x32;
  if(dIdwa.data()) {
    dIdwa(0,0)=x36 + x40;
    dIdwa(1,0)=x29 + x30 + x31;
    dIdwa(2,0)=x17;
    dIdwa(0,1)=x28;
    dIdwa(1,1)=x36 + x44;
    dIdwa(2,1)=x10 + x12 + x8;
    dIdwa(0,2)=x19 + x21 + x23;
    dIdwa(1,2)=x6;
    dIdwa(2,2)=x40 + x44;
  }
  if(dIdwb.data()) {
    dIdwb(0,0)=x45 + x46;
    dIdwb(1,0)=-x25 - x26 - x27;
    dIdwb(2,0)=x24;
    dIdwb(0,1)=x32;
    dIdwb(1,1)=x45 + x47;
    dIdwb(2,1)=-x1 - x3 - x5;
    dIdwb(0,2)=-x14 - x15 - x16;
    dIdwb(1,2)=x13;
    dIdwb(2,2)=x46 + x47;
  }
  return I;
}
//PBCentroidBodyDynamicsGradientInfo
template <typename T>
PBCentroidBodyDynamicsGradientInfo<T>::PBCentroidBodyDynamicsGradientInfo() {}
template <typename T>
PBCentroidBodyDynamicsGradientInfo<T>::PBCentroidBodyDynamicsGradientInfo
(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
 std::shared_ptr<DeformedEnvironment<T>> DEnv) {
  init(body,NContact,g,dt,DEnv);
}
template <typename T>
PBCentroidBodyDynamicsGradientInfo<T>& PBCentroidBodyDynamicsGradientInfo<T>::operator=(const PBCentroidBodyDynamicsGradientInfo& other) {
  _RT=other._RT;
  _w=other._w;
  _DT=other._DT;
  _x=other._x;
  _f=other._f;
  _XTX=other._XTX;
  _g=other._g;
  _dt=other._dt;
  _M=other._M;
  resetPtr();
  return *this;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::setDeformedEnvironment(std::shared_ptr<DeformedEnvironment<T>> DEnv) {
  _DEnv=DEnv;
}
template <typename T>
typename PBCentroidBodyDynamicsGradientInfo<T>::Vec PBCentroidBodyDynamicsGradientInfo<T>::invReset(const ArticulatedBody& body,const Vec& DOF) const {
  ASSERT_MSG(DOF.size()==6,"PBCDM only supports simplified dynamics!")
  ASSERT_MSG(body.joint(0)._typeJoint==Joint::TRANS_3D,"PBCDM only supports TRANS_3D as translational joint!")
  ASSERT_MSG(body.joint(1)._typeJoint==Joint::ROT_3D_EXP ||
             body.joint(1)._typeJoint==Joint::ROT_3D_XYZ,
             "PBCDM only supports ROT_3D_EXP/ROT_3D_XYZ as rotational joint!")
  if(_DEnv) {
    Vec6T invDOF;
    invDOF.template segment<3>(0)=_DEnv->alpha(DOF.template segment<3>(0));
    //rotation
    Mat3T R,RD=_DEnv->rotation(invDOF.template segment<3>(0),NULL);
    if(body.joint(1)._typeJoint==Joint::ROT_3D_EXP) {
      R=expWGradV<T,Vec3T>(DOF.template segment<3>(3));
      invDOF.template segment<3>(3)=invExpW<T>(RD.transpose()*R);
    } else {
      R=eulerX1Y3Z2<T,Vec3T>(DOF.template segment<3>(3),NULL,NULL);
      invDOF.template segment<3>(3)=invEulerX1Y3Z2<T>(RD.transpose()*R);
    }
    return invDOF;
  } else return DOF;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::reset(const ArticulatedBody& body,const Vec& DOF) {
  ASSERT_MSG(DOF.size()==6,"PBCDM only supports simplified dynamics!")
  ASSERT_MSG(body.joint(0)._typeJoint==Joint::TRANS_3D,"PBCDM only supports TRANS_3D as translational joint!")
  ASSERT_MSG(body.joint(1)._typeJoint==Joint::ROT_3D_EXP ||
             body.joint(1)._typeJoint==Joint::ROT_3D_XYZ,
             "PBCDM only supports ROT_3D_EXP/ROT_3D_XYZ as rotational joint!")

  if(_DEnv) {
    DEnvReset(body,DOF);
  } else {
    Vec3T w[3];
    CTR(_RT)=DOF.template segment<3>(0);
    if(body.joint(1)._typeJoint==Joint::ROT_3D_EXP)
      ROT(_RT)=expWGradV<T,Vec3T>(DOF.template segment<3>(3),w);
    else ROT(_RT)=eulerX1Y3Z2<T,Vec3T>(DOF.template segment<3>(3),w,NULL);
    _w.setZero(3,6);
    for(int d=0; d<3; d++)
      _w.col(d+3)=w[d];
    _DT.setIdentity(3,3);
  }
  resetPtr();
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::DEnvReset(const ArticulatedBody& body,const Vec& DOF) {
  //R=Q(J(DOF[0:3]))*R(DOF[3:6])
  //t=forward(DOF[0:3])
  Mat3T J,w,DR;
  Vec3T pos,dw[3];
  ROT(_RT)=_DEnv->rotation(DOF.template segment<3>(0),&w,&pos,&J);
  CTR(_RT)=pos;
  if(body.joint(1)._typeJoint==Joint::ROT_3D_EXP)
    DR=expWGradV<T,Vec3T>(DOF.template segment<3>(3),dw);
  else DR=eulerX1Y3Z2<T,Vec3T>(DOF.template segment<3>(3),dw,NULL);
  _w.setZero(3,6);
  _w.template block<3,3>(0,0)=w;
  for(int d=0; d<3; d++)
    _w.col(d+3)=ROT(_RT)*dw[d];
  ROT(_RT)*=DR;
  _DT=J;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::init
(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
 std::shared_ptr<DeformedEnvironment<T>> DEnv) {
  ASSERT_MSG(body.joint(0)._typeJoint==Joint::TRANS_3D,"PBCDM only supports TRANS_3D as translational joint!")
  ASSERT_MSG(body.joint(1)._typeJoint==Joint::ROT_3D_EXP ||
             body.joint(1)._typeJoint==Joint::ROT_3D_XYZ,
             "PBCDM only supports ROT_3D_EXP/ROT_3D_XYZ as rotational joint!")
  ASSERT_MSG(body.joint(0)._trans == Joint::Mat3X4T::Identity(),"PBCDM only supports Joint._trans==Identity!")
  ASSERT_MSG(body.joint(1)._trans == Joint::Mat3X4T::Identity(),"PBCDM only supports Joint._trans==Identity!")

  _x.setZero(3,NContact);
  _f.setZero(3,NContact);
  resetPtr();

  _XTX=body.joint(body.rootJointId())._MCCT.template cast<T>();
  _g=g;
  _dt=dt;
  _M=body.joint(body.rootJointId())._M;
  _isEXP=body.joint(1)._typeJoint==Joint::ROT_3D_EXP;
  _DEnv=DEnv;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::resetPtr() {
  new (&_RTM)   Mat3X4TM(mapM(_RT));
  new (&_wM)    Mat3XTM(mapM(_w));
  new (&_DTM)   Mat3XTM(mapM(_DT));
  new (&_xM)    Mat3XTM(mapM(_x));
  new (&_fM)    Mat3XTM(mapM(_f));
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::debug
(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
 std::shared_ptr<DeformedEnvironment<T>> DEnv) {
  DEFINE_NUMERIC_DELTA_T(T)
  PBCentroidBodyDynamicsGradientInfo GNN(body,NContact,g,dt,DEnv);
  PBCentroidBodyDynamicsGradientInfo GN(body,NContact,g,dt,DEnv);
  PBCentroidBodyDynamicsGradientInfo G(body,NContact,g,dt,DEnv);
  PBCentroidBodyDynamicsGradientInfo GRef(body,NContact,g,dt,DEnv);
  Vec dI=Vec::Random(6);
  Vec INN=Vec::Random(6);
  Vec IN=Vec::Random(6);
  Vec I=Vec::Random(6);
  GNN.reset(body,INN);
  GN.reset(body,IN);
  G.reset(body,I);
  Mat3XT x=G._x.setRandom(3,NContact);
  Mat3XT f=G._f.setRandom(3,NContact);
  Mat6T DCDINN,DCDIN,DCDI;
  Mat6XT DCDX=Mat6XT::Random(6,NContact*3);
  Mat6XT DCDF=Mat6XT::Random(6,NContact*3);
  Vec C=G.C(GNN,GN,mapM(DCDINN),mapM(DCDIN),mapM(DCDI),mapM(DCDX),mapM(DCDF)),C2;

  //invDOF
  Vec IRef=GRef.invReset(body,G.getDOF());
  DEBUG_GRADIENT("invReset",IRef.norm(),(IRef-I).norm())

  //DCDINN
  GNN.reset(body,INN+dI*DELTA);
  C2=G.C(GNN,GN);
  DEBUG_GRADIENT("DCDINN",(DCDINN*dI).norm(),(DCDINN*dI-(C2-C)/DELTA).norm())
  GNN.reset(body,INN);

  //DCDIN
  GN.reset(body,IN+dI*DELTA);
  C2=G.C(GNN,GN);
  DEBUG_GRADIENT("DCDIN",(DCDIN*dI).norm(),(DCDIN*dI-(C2-C)/DELTA).norm())
  GN.reset(body,IN);

  //DCDI
  G.reset(body,I+dI*DELTA);
  C2=G.C(GNN,GN);
  DEBUG_GRADIENT("DCDI",(DCDI*dI).norm(),(DCDI*dI-(C2-C)/DELTA).norm())
  G.reset(body,I);

  //DCDX
  Vec dx=Vec::Random(NContact*3);
  G._x=x+Eigen::Map<const Mat3XT>(dx.data(),3,NContact)*DELTA;
  C2=G.C(GNN,GN);
  DEBUG_GRADIENT("DCDX",(DCDX*dx).norm(),(DCDX*dx-(C2-C)/DELTA).norm())
  G._x=x;

  //DCDF
  Vec df=Vec::Random(NContact*3);
  G._f=f+Eigen::Map<const Mat3XT>(df.data(),3,NContact)*DELTA;
  C2=G.C(GNN,GN);
  DEBUG_GRADIENT("DCDF",(DCDF*df).norm(),(DCDF*df-(C2-C)/DELTA).norm())
  G._f=f;
}
template <typename T>
void PBCentroidBodyDynamicsGradientInfo<T>::debugContinuous(const ArticulatedBody& body,int NContact) {
  DEFINE_NUMERIC_DELTA_T(float)
  PBCentroidBodyDynamicsGradientInfo GNN(body,NContact,Vec3T::Zero(),DELTA);
  PBCentroidBodyDynamicsGradientInfo GN(body,NContact,Vec3T::Zero(),DELTA);
  PBCentroidBodyDynamicsGradientInfo G(body,NContact,Vec3T::Zero(),DELTA);
  GNN.reset(body,Vec::Random(6));
  GN.reset(body,Vec::Random(6));
  G.reset(body,Vec::Random(6));
  Vec3T w0=Vec3T::Random();
  Vec3T w=Vec3T::Random();
  Vec3T dw=Vec3T::Random();
  Vec3T v=Vec3T::Random();
  Vec3T dv=Vec3T::Random();
  ROT(GNN._RT)=expWGradV<T,Vec3T>(w0);
  ROT(GN._RT)=expWGradV<T,Vec3T>(w*DELTA)*expWGradV<T,Vec3T>(w0);
  ROT(G._RT)=expWGradV<T,Vec3T>((w+dw*DELTA)*DELTA)*
             expWGradV<T,Vec3T>(w*DELTA)*
             expWGradV<T,Vec3T>(w0);
  CTR(GN._RT)=CTR(GNN._RT)+v*DELTA;
  CTR(G._RT)=CTR(GN._RT)+(v+dv*DELTA)*DELTA;
  G._x.setZero(3,NContact);
  G._f.setZero(3,NContact);
  Vec C=G.C(GNN,GN);

  Mat3T R0=ROT(GNN._RT);
  Mat3T R1=ROT(GN._RT);
  Mat3T R2=ROT(G._RT);
  Mat3T ddR=(cross(dw)+cross(w)*cross(w))*R0;
  Mat3T ddRRef=(R2-2*R1+R0)/DELTA/DELTA;
  DEBUG_GRADIENT("ddR",ddR.norm(),(ddR-ddRRef).norm())

  //CRef
  Vec6T CRef;
  Mat3T I0=Mat3T::Identity()*G._XTX.trace()-G._XTX;
  Mat3T I=R0*I0*R0.transpose();
  CRef.template segment<3>(0)=I*dw+w.cross(I*w);
  CRef.template segment<3>(3)=dv*G._M;
  DEBUG_GRADIENT("C0",C.template segment<3>(0).norm(),(C-CRef).template segment<3>(0).norm())
  DEBUG_GRADIENT("C1",C.template segment<3>(3).norm(),(C-CRef).template segment<3>(3).norm())
}
//instance
template struct PBCentroidBodyDynamicsGradientInfoMap<FLOAT>;
template struct PBCentroidBodyDynamicsGradientInfo<FLOAT>;
#ifdef FORCE_ADD_DOUBLE_PRECISION
template struct PBCentroidBodyDynamicsGradientInfoMap<double>;
template struct PBCentroidBodyDynamicsGradientInfo<double>;
#endif
}
