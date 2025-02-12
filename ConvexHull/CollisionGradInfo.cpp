#include "CollisionGradInfo.h"
#include "ConvexHullDistanceEnergy.h"
#include <Environment/MeshExact.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
//CollisionGradInfo
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo() {}
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo(const ArticulatedBody& body,const Vec& theta) {
  reset(body,theta);
}
template <typename T>
void CollisionGradInfo<T>::reset(const ArticulatedBody& body,const Vec& theta) {
  _HTheta.setZero(body.nrDOF(),body.nrDOF());
  _HThetaL.setZero(body.nrDOF(),body.nrDOF());
  _HThetaLL.setZero(body.nrDOF(),body.nrDOF());
  _HThetaPTarget.setZero(body.nrDOF(),body.nrDOF());
  _HThetaDTarget.setZero(body.nrDOF(),body.nrDOF());
  _HThetaN=MatX4T::Zero(body.nrDOF(),4);
  _HThetaU=MatX3T::Zero(body.nrDOF(),3);
  _DTG.setZero(3,4*body.nrJ());
  _info.reset(body,theta);
  //data
  _polytopes.resize(body.nrJ());
  Eigen::Matrix<int,2,1> off(0,0);
  for(int i=0; i<body.nrJ(); i++) {
    std::shared_ptr<MeshExact> m=std::dynamic_pointer_cast<MeshExact>(body.joint(i)._mesh);
    if(!m)
      _polytopes[i]=GJKPolytope<T>();
    else {
      if(_polytopes[i].jid()==-1)
        _polytopes[i]=GJKPolytope<T>(i,body,*this);
      else _polytopes[i].resetGlobalVss(this);
      off[1]=off[0]+(int)m->vss().size();
      _polytopes[i].setVertexId(off);
      off[0]+=m->vss().size();
    }
  }
  _nrVex=off[0];
  _HThetaD.setZero(body.nrDOF(),_nrVex*3);
  _HPos.setZero(body.nrDOF(),_nrVex*3);
}
template <typename T>
void CollisionGradInfo<T>::getTransformation(const ArticulatedBody& body,int jid,Mat3T& R,Vec3T& t) const {
  ASSERT_MSG(jid>=0 && jid<body.nrJ(),"Invalid joint id when calling getTransformation!")
  R=ROTI(_info._TM,jid);
  t=CTRI(_info._TM,jid);
}
template <typename T>
void CollisionGradInfo<T>::getJacobian(const ArticulatedBody& body,int jid,Mat3XT& JV,Mat3XT& JW) const {
  ASSERT_MSG(jid>=0 && jid<body.nrJ(),"Invalid joint id when calling getJacobian!")
  JV.setZero(3,body.nrDOF());
  JW.setZero(3,body.nrDOF());
  _info.JCSparse(body,jid,[&](int col,const Vec3T& v) {
    JV.col(col)=v;
  });
  _info.JRSparse(body,jid,[&](int col,const Vec3T& w) {
    JW.col(col)=w;
  });
}
template <typename T>
void CollisionGradInfo<T>::getJacobianDeriv(const ArticulatedBody& body,int jid,std::vector<Mat3XT>& dJVdq,std::vector<Mat3XT>& dJWdq) const {
  ASSERT_MSG(jid>=0 && jid<body.nrJ(),"Invalid joint id when calling getJacobianDeriv!")
  //jacobian
  Mat3XT JW=Mat3XT::Zero(3,body.nrDOF());
  _info.JRSparse(body,jid,[&](int col,const Vec3T& w) {
    JW.col(col)=w;
  });
  //derivative
  dJVdq.resize(body.nrDOF(),Mat3XT::Zero(3,body.nrDOF()));
  dJWdq.resize(body.nrDOF(),Mat3XT::Zero(3,body.nrDOF()));
  Mat3XT G=Mat3XT::Zero(3,body.nrJ()*4);
  for(int d=0; d<3; d++) {
    //translation
    G.setZero();
    CTRI(G,jid).setUnit(d);
    _info.toolB(body,mapM(G),[&](int row,int col,T val) {
      dJVdq[row].col(col)[d]+=val;
    });
    //rotation double derivative:
    //trace([w_{1,2}]R\bar{R}^T)+trace([w_1][w_2]R\bar{R}^T)
    G.setZero();
    ROTI(G,jid)=cross<T>(Vec3T::Unit(d))*0.5f*ROTI(_info._TM,jid);
    _info.toolB(jid,body,mapM(G),[&](int row,int col,T val) {
      dJWdq[row].col(col)[d]+=val;
    });
  }
  //Subtract trace([w_1][w_2]R\bar{R}^T)
  for(int row=0; row<body.nrDOF(); row++)
    for(int col=0; col<body.nrDOF(); col++)
      dJWdq[row].col(col)+=invCrossMatTrace<T>(cross<T>(JW.col(col))*cross<T>(JW.col(row)))*0.5f;
}
template <typename T>
void CollisionGradInfo<T>::debugJacobianDeriv(const ArticulatedBody& body) {
  DEFINE_NUMERIC_DELTA_T(T)
  Vec x=Vec::Random(body.nrDOF());
  Vec dx=Vec::Random(body.nrDOF());
  CollisionGradInfo<T> pos(body,x);
  CollisionGradInfo<T> pos2(body,x+dx*DELTA);
  //pos._info.debug(body);

  //DTG
  Vec DTG=Vec::Zero(body.nrDOF());
  Vec DTGRef=Vec::Zero(body.nrDOF());
  Mat3XT G=Mat3XT::Random(3,4*body.nrJ());
  for(int jid=0; jid<body.nrJ(); jid++) {
    Mat3XT JV,JW;
    pos.getJacobian(body,jid,JV,JW);
    DTG+=JV.transpose()*CTRI(G,jid);
    for(int c=0; c<JW.cols(); c++)
      DTG[c]+=(cross<T>(JW.col(c))*ROTI(pos._info._TM,jid)*ROTI(G,jid).transpose()).trace();
  }
  pos._info.DTG(body,mapM(G),mapV(DTGRef));
  DEBUG_GRADIENT("DTG",DTG.norm(),(DTG-DTGRef).norm())

  //dJdq
  for(int jid=0; jid<body.nrJ(); jid++) {
    Mat3XT JV,JW,JV2,JW2;
    pos.getJacobian(body,jid,JV,JW);
    pos2.getJacobian(body,jid,JV2,JW2);

    Mat3XT dJV,dJW;
    std::vector<Mat3XT> dJVdq,dJWdq;
    dJV.setZero(JV.rows(),JV.cols());
    dJW.setZero(JW.rows(),JW.cols());
    pos.getJacobianDeriv(body,jid,dJVdq,dJWdq);
    for(int i=0; i<JV.cols(); i++)
      dJV+=dJVdq[i]*dx[i];
    for(int i=0; i<JW.cols(); i++)
      dJW+=dJWdq[i]*dx[i];

    DEBUG_GRADIENT("JV",dJV.norm(),(dJV-(JV2-JV)/DELTA).norm())
    DEBUG_GRADIENT("JW",dJW.norm(),(dJW-(JW2-JW)/DELTA).norm())
  }
}
//instance
template struct CollisionGradInfo<FLOAT>;
}
