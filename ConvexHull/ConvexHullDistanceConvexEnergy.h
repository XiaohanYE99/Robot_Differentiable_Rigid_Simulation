#ifndef CONVEX_HULL_DISTANCE_CONVEX_ENERGY_H
#define CONVEX_HULL_DISTANCE_CONVEX_ENERGY_H

#include "ConvexHullDistanceEnergy.h"

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH=typename HigherPrecisionTraits<T>::TH>
class CCBarrierConvexEnergy : public CCBarrierEnergy<T,PFunc,TH> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef Eigen::Matrix<TH,3,1> Vec3TH;
  typedef Eigen::Matrix<TH,4,1> Vec4TH;
  typedef Eigen::Matrix<TH,3,3> Mat3TH;
  typedef Eigen::Matrix<TH,4,4> Mat4TH;
  using CCDistanceEnergy<T>::_p1;
  using CCDistanceEnergy<T>::_p2;
  using CCBarrierEnergy<T,PFunc,TH>::_p;
  using CCBarrierEnergy<T,PFunc,TH>::_d0;
  using CCBarrierEnergy<T,PFunc,TH>::_d0Half;
  using CCBarrierEnergy<T,PFunc,TH>::_coef;
  using CCBarrierEnergy<T,PFunc,TH>::_implicit;
  using CCBarrierEnergy<T,PFunc,TH>::_output;
  using CCBarrierEnergy<T,PFunc,TH>::_x;
  using CCBarrierEnergy<T,PFunc,TH>::initialize;
  using CCBarrierEnergy<T,PFunc,TH>::debugEnergy;
  CCBarrierConvexEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit=true);
  bool initialize(Vec4T* res,const ArticulatedBody* body) override;
  static void debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
 private:
  bool energy(const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const override;
  bool optimize(Vec4TH& x,TH& E,Vec4TH& G,Mat4TH& H) const override;
  //gradient/hessian
  void sensitivity(const Vec4T& x,const Vec4T& G,const Mat4T& H,Mat4T& invH) const override;
};
}
#endif
