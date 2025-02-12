#ifndef CONVEX_HULL_DISTANCE_FRICTION_ENERGY_H
#define CONVEX_HULL_DISTANCE_FRICTION_ENERGY_H

#include "ConvexHullDistanceConvexEnergy.h"

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH=typename HigherPrecisionTraits<T>::TH>
class CCBarrierFrictionEnergy : public CCBarrierEnergy<T,PFunc,TH> {
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
  CCBarrierFrictionEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const GJKPolytope<T>& pl1,const GJKPolytope<T>& pl2,const Vec4T& x,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,T dt,bool implicit=true);
  bool initialize(Vec4T* res,const ArticulatedBody* body) override;
  static void debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
  bool eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta);
  bool evalBackward(const ArticulatedBody* body,CollisionGradInfo<T>* grad,CollisionGradInfo<T>* Pos,
                    std::vector<Mat3X4T>* DNDX,std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2);
  virtual void debugEnergyP(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,T perturbRange)const;
  virtual void debugEnergyN(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,T perturbRange)const;
  virtual void debugEnergy(const Vec3TH& x) const;
  T fri() const;
  T& fri();
 private:
  bool energyP(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,TH& E,Vec3TH* G,Mat3TH* H) const;
  bool energyN(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,TH& E,Vec3TH* G,Mat3TH* H) const;
  bool energy(const Vec3TH& u,TH& E,Vec3TH* G,Mat3TH* H) const ;
  bool optimize(Vec3TH& u,TH& E,Vec3TH& G,Mat3TH& H) const;
  //gradient/hessian
  void sensitivity(const Vec3T& u,const Vec3T& G,const Mat3T& H,Mat3T& invH) const ;
  void energyPDTGH(const Vec3T& v,const Vec3T& vLast,const Vec3T& vl,const Vec3T& u,const Vec4T& x,Mat3X4T* DTG,
                   const Vec3T& Rvl,const Vec3T* Rvll,Mat3T* LRH,Mat3T* wLRH,Mat3T* LRH1,Mat3T* wLRH1,Mat3X4T* LRH2,Mat3X4T* wLRH2,
                   Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,Mat3T *Hxx=NULL) const;
  void energyNDTGH(const Vec3T& v,const Vec3T& vLast,const Vec3T& vl,const Vec3T& u,const Vec4T& x,Mat3X4T* DTG,
                   const Vec3T& Rvl,const Vec3T* Rvll,Mat3T* LRH,Mat3T* wLRH,Mat3T* LRH1,Mat3T* wLRH1,Mat3X4T* LRH2,Mat3X4T* wLRH2,
                   Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,Mat3T *Hxx=NULL) const;
  void energyPDUDN(const Vec3T& v,const Vec3T& vLast,const Vec3T& u,const Vec4T& x,Mat3X4T* LRH) const;
  void energyNDUDN(const Vec3T& v,const Vec3T& vLast,const Vec3T& u,const Vec4T& x,Mat3X4T* LRH) const;
  void computeDTGH(const ArticulatedBody& body,CollisionGradInfo<T>& info,
                   const Vec3T& u,const Vec4T& x,const Vec3T* G,const Mat3T* H) const;
  void computeHLBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,CollisionGradInfo<T>& infoL,
                         const Vec3T& u,const Vec4T& x,const Vec3T* G,const Mat3T* H,const std::vector<Mat3X4T>* DNDX) const;
  void computeHBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,
                        const Vec3T& u,const Vec4T& x,const Vec3T& G,const Mat3T& H,
                        std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) const;
  //virtual void contractHTheta(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,
  //                    const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const;
  const GJKPolytope<T>& _pl1;
  const GJKPolytope<T>& _pl2;
  T _dt,_eps,_fri;
  Vec3TH _u;
};
}
#endif
