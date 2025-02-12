#ifndef COLLISION_GRAD_INFO_H
#define COLLISION_GRAD_INFO_H

#include "GJKPolytope.h"
#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/Utils.h>

namespace PHYSICSMOTION {
template <typename T>
struct CollisionGradInfo {
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  CollisionGradInfo();
  CollisionGradInfo(const ArticulatedBody& body,const Vec& theta);
  void reset(const ArticulatedBody& body,const Vec& theta);
  void getTransformation(const ArticulatedBody& body,int jid,Mat3T& R,Vec3T& t) const;
  void getJacobian(const ArticulatedBody& body,int jid,Mat3XT& JV,Mat3XT& JW) const;
  void getJacobianDeriv(const ArticulatedBody& body,int jid,std::vector<Mat3XT>& dJVdq,std::vector<Mat3XT>& dJWdq) const;
  static void debugJacobianDeriv(const ArticulatedBody& body);
  //data
  int _nrVex;
  std::vector<Vec3T> _centre;
  MatT _HTheta,_HThetaL,_HThetaLL;      //theta gradient
  MatT _HThetaD;                        //vertex gradient
  MatT _HThetaPTarget,_HThetaDTarget;   //PD gradient
  MatT _HPos;                           //pos gradient
  MatT _HThetaDesign;                   //design gradient
  MatX4T _HThetaN;
  MatX3T _HThetaU;
  Mat3XT _DTG;
  std::vector<GJKPolytope<T>> _polytopes;
  PBDArticulatedGradientInfo<T> _info;
};
}
#endif
