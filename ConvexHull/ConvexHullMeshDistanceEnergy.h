#ifndef CONVEX_HULL_MESH_DISTANCE_ENERGY_H
#define CONVEX_HULL_MESH_DISTANCE_ENERGY_H

#include "ConvexHullDistanceEnergy.h"

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH=typename HigherPrecisionTraits<T>::TH>
class CCBarrierMeshEnergy : public CCBarrierEnergy<T,PFunc,TH> {
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
  typedef GJKPolytope<T> const* GJKPolytopePtr;
  struct MPair {
    Mat3T _Mww=Mat3T::Zero();
    Mat3T _Mwt=Mat3T::Zero();
    Mat3T _Mtw=Mat3T::Zero();
    Mat3T _Mtt=Mat3T::Zero();
    Mat3T _DMwt=Mat3T::Zero();
    Mat3T _DMww=Mat3T::Zero();
  };
  struct MAll {
    MPair _m11,_m12,_m21,_m22;
  };
  CCBarrierMeshEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit=true);
  virtual bool eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,Vec4T* x=NULL);
  virtual bool evalbackward(T *E,const ArticulatedBody* body,CollisionGradInfo<T>* grad);
  static void debugGradient(const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
 protected:
  bool evalBF(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward=false) const;
  bool evalBvh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward=false) const;
  T evalBsh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward=false) const;
  bool evalEE(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  bool evalVT(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  //bool evalVT(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  void computeDTGH(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const;
  void computeHBackward(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const;
  void contractHAll(const ArticulatedBody& body,CollisionGradInfo<T>& grad,const MAll& m) const;
  bool _useBVH=true;
  bool _useLRI=true;
};
}
#endif
