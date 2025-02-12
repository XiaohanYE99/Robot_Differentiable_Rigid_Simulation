#ifndef CONVEX_HULL_MESH_DISTANCE_FRICTION_ENERGY_H
#define CONVEX_HULL_MESH_DISTANCE_FRICTION_ENERGY_H

#include "ConvexHullMeshDistanceEnergy.h"

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH=typename HigherPrecisionTraits<T>::TH>
class CCBarrierMeshFrictionEnergy : public CCBarrierMeshEnergy<T,PFunc,TH> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef Eigen::Matrix<TH,2,1> Vec2TH;
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
  typedef GJKPolytope<T> const* GJKPolytopePtr;
  using typename CCBarrierMeshEnergy<T,PFunc,TH>::MAll;
  using typename CCBarrierMeshEnergy<T,PFunc,TH>::MPair;
  struct FrictionTerm {
    FrictionTerm();
    FrictionTerm(const Eigen::Matrix<int,4,1>& pss,const Eigen::Matrix<int,4,1>& vid,const Vec3T V[4],const Vec4T& bary,const Vec3T& n,T D);
    Eigen::Matrix<int,4,1> _pss;
    Eigen::Matrix<int,4,1> _vid;
    Vec3T _V[4];
    Vec4T _bary;
    Vec3T _n;
    T _D;
  };
  CCBarrierMeshFrictionEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const GJKPolytope<T>& pl1,const GJKPolytope<T>& pl2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,T dt,const std::vector<FrictionTerm>* terms=NULL);
  static void debugGradient(const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
  bool eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta);
  bool evalbackward(T *E,const ArticulatedBody* body,CollisionGradInfo<T>* grad);
  const std::vector<FrictionTerm>& terms() const;
 private:
  void computeDTGH(const FrictionTerm& term,const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec3T& G,const Mat3T& H,MAll& m) const;
  void computeHBackward(const FrictionTerm& term,const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec3T& G,const Mat3T& H,MAll& m) const;
  //data
  std::vector<FrictionTerm> _terms;
  const GJKPolytope<T>& _pl1;
  const GJKPolytope<T>& _pl2;
  const T _dt,_eps,_fri;
};
}
#endif
