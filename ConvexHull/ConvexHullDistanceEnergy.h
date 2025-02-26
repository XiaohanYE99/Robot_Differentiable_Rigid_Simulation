#ifndef CONVEX_HULL_DISTANCE_ENERGY_H
#define CONVEX_HULL_DISTANCE_ENERGY_H

#include "Barrier.h"
#include "GJKPolytope.h"
#include "CollisionGradInfo.h"
#include <Environment/MeshExact.h>
#include <Utils/VTKWriter.h>

namespace PHYSICSMOTION {
template <typename T>
struct HigherPrecisionTraits {
#ifdef QUADMATH_SUPPORT
  typedef double TH;
#else
  typedef double TH;
#endif
};
template <>
struct HigherPrecisionTraits<float128> {
  typedef mpfr_float TH;
};
template <>
struct HigherPrecisionTraits<mpfr_float> {
  typedef mpfr_float TH;
};
template <typename T>
class CCDistanceEnergy {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  CCDistanceEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2);
  virtual bool eval(T* E,Vec12T* G,Mat12T* H);
 protected:
  const GJKPolytope<T>& _p1;
  const GJKPolytope<T>& _p2;
};
template <typename T,typename PFunc,typename TH=typename HigherPrecisionTraits<T>::TH>
class CCBarrierEnergy : public CCDistanceEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef Eigen::Matrix<TH,3,1> Vec3TH;
  typedef Eigen::Matrix<TH,4,1> Vec4TH;
  typedef Eigen::Matrix<TH,3,3> Mat3TH;
  typedef Eigen::Matrix<TH,4,4> Mat4TH;
  using CCDistanceEnergy<T>::_p1;
  using CCDistanceEnergy<T>::_p2;
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
  struct GPair {
    Vec3T _w=Vec3T::Zero();
    Vec3T _t=Vec3T::Zero();
  };
  struct GAll {
    GPair _g1,_g2;
  };
  CCBarrierEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit=true);
  virtual ~CCBarrierEnergy() {}
  bool update(Vec4T* res);
  Eigen::Matrix<T,4,1> getX() const;
  virtual bool initialize(Vec4T* res,const ArticulatedBody* body);
  void initialize(const Vec4T& x);
  virtual bool eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,Vec4T* x=NULL);
  bool evalLRI(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,
               Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g,Vec4T* x=NULL);
  virtual bool evalBackward(const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2);
  static void debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0=0.01,bool output=false);
  static void debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0=0.01,bool output=false);
  virtual void debugEnergyP(const Vec3TH& v,const Vec4TH& x,T perturbRange)const;
  virtual void debugEnergyN(const Vec3TH& v,const Vec4TH& x,T perturbRange)const;
  void setOutput(bool output);
 protected:
  bool energyP(const Vec3TH& v,const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const;
  bool energyN(const Vec3TH& v,const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const;
  virtual bool energy(const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const;
  virtual bool optimize(Vec4TH& x,TH& E,Vec4TH& G,Mat4TH& H) const;
  virtual void debugEnergy(const Vec4TH& x) const;
  //gradient/hessian
  virtual void sensitivity(const Vec4T& x,const Vec4T& G,const Mat4T& H,Mat4T& invH) const;
  virtual void energyPDTGH(const Vec3T& v,const Vec3T& vl,const Vec4T& x,Mat3X4T* DTG,
                           const Vec3T& Rvl,Mat3X4T* LRH,Mat3X4T* wLRH,
                           Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,GAll* g=NULL) const;
  virtual void energyNDTGH(const Vec3T& v,const Vec3T& vl,const Vec4T& x,Mat3X4T* DTG,
                           const Vec3T& Rvl,Mat3X4T* LRH,Mat3X4T* wLRH,
                           Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,GAll* g=NULL) const;
  void computeDTGHLRI(const ArticulatedBody& body,CollisionGradInfo<T>& info,
                      const Vec4T& x,const Vec4T* G,const Mat4T* H,std::vector<Mat3X4T>* DNDX,
                      Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g) const;
  virtual void computeDTGH(const ArticulatedBody& body,CollisionGradInfo<T>& info,
                           const Vec4T& x,const Vec4T* G,const Mat4T* H,std::vector<Mat3X4T>* DNDX) const;
  void computeHBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,
                        const Vec4T& x,const Vec4T& G,const Mat4T& H,
                        std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) const;
  virtual void contractHTheta(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,
                              const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const;
  virtual void contractHThetaL(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,CollisionGradInfo<T>& infoL,
                               const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const;
  void contractHBackward(int k,const ArticulatedBody& body,const CollisionGradInfo<T>& info,
                         MatX4T& HThetaX,const Mat3X4T& Mwx,const Mat3X4T& Mtx) const;
  void contractHBackward(int k,const ArticulatedBody& body,const CollisionGradInfo<T>& info,
                         MatX3T& HThetaX,const Mat3T& Mwt,const Mat3T& Mtt) const;
 protected:
  const PFunc& _p;
  T _d0,_d0Half,_coef;
  bool _implicit;
  bool _output;
  Vec4TH _x;
  const CollisionGradInfo<T>* _grad;
};
}
#endif
