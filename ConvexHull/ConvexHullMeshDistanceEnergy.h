#ifndef CONVEX_HULL_MESH_DISTANCE_ENERGY_H
#define CONVEX_HULL_MESH_DISTANCE_ENERGY_H

#include "ConvexHullDistanceEnergy.h"
#include "ConvexHullDistanceConvexEnergy.h"
#include <Environment/ContactGenerator.h>

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
  using CCBarrierEnergy<T,PFunc,TH>::_grad;
  typedef GJKPolytope<T> const* GJKPolytopePtr;
  typedef ContactGenerator::ContactManifold ContactManifold;
  using typename CCBarrierEnergy<T,PFunc,TH>::MAll;
  using typename CCBarrierEnergy<T,PFunc,TH>::MPair;
  using typename CCBarrierEnergy<T,PFunc,TH>::GAll;
  using typename CCBarrierEnergy<T,PFunc,TH>::GPair;
  CCBarrierMeshEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit=false);
  virtual bool eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,std::vector<ContactManifold>* m=NULL);
  virtual bool evalbackward(T *E,const ArticulatedBody* body,CollisionGradInfo<T>* grad);
  static void debugGradient(const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
 protected:
  bool evalBF(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward=false) const;
  bool evalBvh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward=false) const;
  bool evalBsh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<ContactManifold>* ml,bool backward=false) const;
  bool evalEE(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  bool evalVT(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  //bool evalVT(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward=false) const;
  void computeDTGH(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const;
  void computeHBackward(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const;
  void contractHAll(const ArticulatedBody& body,CollisionGradInfo<T>& grad,const MAll& m) const;
  void contractMAll(MAll& m,Mat3T Rxi,Mat3T Rxj,Mat3T rxi,Mat3T rxj,Mat6T H) const;
  void contractGAll(GAll& g,Mat3T Rxi,Mat3T Rxj,Vec3T L) const;
  bool ComputePotential(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,int id1, int id2, T* P,
                        Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g,CollisionGradInfo<T>* grad,
                        const ArticulatedBody& body,std::vector<ContactManifold>* ml,bool* flag) const;
  void addMAll(MAll& m1,MAll& m2,T alpha) const;
  void addGAll(GAll& g1,GAll& g2,T alpha) const;
  void clearMAll(MAll& m) const;
  void clearGAll(GAll& g) const;
  void mergeGAll(GAll& g1,GAll& g2,MAll& m) const;
  bool _useBVH=true;
  bool _useLRI=true;
  T _ep=.1;
  Vec4T _h;
};
}
#endif
