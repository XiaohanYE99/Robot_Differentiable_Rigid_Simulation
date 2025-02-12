#ifndef DISTANCE_ENERGY_H
#define DISTANCE_ENERGY_H

#include "Barrier.h"
#include <Utils/SparseUtils.h>
#include <Utils/CrossSpatialUtils.h>
#include <math.h>

namespace PHYSICSMOTION {
template <typename T>
class VTDistanceEnergy {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  VTDistanceEnergy(const Vec3T triA,const Vec3T& triB,const Vec3T& triC,const Vec3T& point);
  VTDistanceEnergy(const Vec3T tri[3],const Vec3T& point);
  virtual bool eval(T* E,Vec12T* G,Mat12T* H);
  static const SMatT& getIF();
  static const SMatT& getIE();
  bool debug(const Vec3T tri[3],const Vec3T& point,T perturbRange,const Eigen::Matrix<char,2,1>& feat);
  void debug(T perturbRange,const Eigen::Matrix<char,2,1>& feat);
 private:
  Eigen::Matrix<char,2,1> _feat;
  Vec3T _tri[3];
  Vec3T _point;
};
template <typename T>
class EEDistanceEnergy {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  EEDistanceEnergy(const Vec3T edge1A,const Vec3T edge1B,const Vec3T edge2A,const Vec3T edge2B,bool JTJApprox=false);
  EEDistanceEnergy(const Vec3T edge1[2],const Vec3T edge2[2],bool JTJApprox=false);
  virtual bool eval(T* E,Vec12T* G,Mat12T* H);
  bool debug(const Vec3T edge1[2],const Vec3T edge2[2],T perturbRange,const Eigen::Matrix<char,2,1>& feat);
  void debug(T perturbRange,const Eigen::Matrix<char,2,1>& feat);
  static const SMatT& getIEE();
 protected:
  Eigen::Matrix<char,2,1> _feat;
  Vec3T _edge1[2],_edge2[2];
  T _mollifierCoef;
  bool _JTJApprox;
};
template <typename T,typename PFunc>
class VTBarrierEnergy : public VTDistanceEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  VTBarrierEnergy(const Vec3T triA,const Vec3T& triB,const Vec3T& triC,const Vec3T& point,const PFunc& p,T d0,T coef);
  VTBarrierEnergy(const Vec3T tri[3],const Vec3T& point,const PFunc& p,T d0,T coef);
  bool eval(T* E,Vec12T* G,Mat12T* H) override;
 private:
  const PFunc& _p;
  T _d0,_coef;
};
template <typename T,typename PFunc>
class EEBarrierEnergy : public EEDistanceEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  using EEDistanceEnergy<T>::_edge1;
  using EEDistanceEnergy<T>::_edge2;
  using EEDistanceEnergy<T>::_JTJApprox;
  EEBarrierEnergy(const Vec3T edge1A,const Vec3T edge1B,const Vec3T edge2A,const Vec3T edge2B,bool JTJApprox,const PFunc& p,T d0,T coef,T mollifierCoef=1e-2);
  EEBarrierEnergy(const Vec3T edge1[2],const Vec3T edge2[2],bool JTJApprox,const PFunc& p,T d0,T coef,T mollifierCoef=1e-2);
  bool eval(T* E,Vec12T* G,Mat12T* H) override;
 private:
  bool mollifier(T* E,Vec12T* G,Mat12T* H) const;
  const PFunc& _p;
  T _d0,_coef,_mollifierCoef;
};
}
#endif

