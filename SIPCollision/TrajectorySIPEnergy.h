#ifndef TRAJECTORY_SIP_ENERGY_H
#define TRAJECTORY_SIP_ENERGY_H

#include "ThetaTrajectory.h"
#include <Articulated/ArticulatedBody.h>

namespace PHYSICSMOTION {
template <typename T>
class SIPEnergy {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  typedef Eigen::SparseVector<T,0,int> SVec;
  typedef std::function<void(const T)> EFunc;
  typedef std::function<void(int,const Vec&)> GFunc;
  typedef std::function<void(int,int,const MatT&)> HFunc;
  virtual bool eval(EFunc* E,GFunc* G,HFunc* H)=0;
  virtual bool eval(T* E,Vec* G,MatT* H);
  static inline void parallelAdd(T& lhs,T val) {
#ifndef FORCE_ADD_DOUBLE_PRECISION
    OMP_ATOMIC_
#else
    OMP_CRITICAL_
#endif
    lhs+=val;
  }
  template <int R>
  static inline void parallelAdd(Eigen::Matrix<T,R,1>& DG,int off,const Vec& val) {
    for(int i=0; i<val.rows(); i++)
      parallelAdd(DG[off+i],val[i]);
  }
  template <int R,int C>
  static inline void parallelAdd(Eigen::Matrix<T,R,C>& DH,int offr,int offc,const MatT& val) {
    for(int r=0; r<val.rows(); r++)
      for(int c=0; c<val.cols(); c++)
        parallelAdd(DH(offr+r,offc+c),val(r,c));
  }
  template <int R,int C>
  static inline void parallelAdd(Eigen::Matrix<T,R,C>& DH,int offr,int offc,const Eigen::MatrixBase<Mat3X4T>& val) {
    for(int r=0; r<val.rows(); r++)
      for(int c=0; c<val.cols(); c++)
        parallelAdd(DH(offr+r,offc+c),val(r,c));
  }
  template <int R,int C>
  static inline void parallelAdd(Eigen::Matrix<T,R,C>& DH,int offr,int offc,const Eigen::MatrixBase<Mat3T>& val) {
    for(int r=0; r<val.rows(); r++)
      for(int c=0; c<val.cols(); c++)
        parallelAdd(DH(offr+r,offc+c),val(r,c));
  }
  template <int R>
  static inline void parallelAdd(Eigen::Matrix<T,R,1>& DG,int off,const SVec& val) {
    for(typename SVec::InnerIterator it(val); it; ++it)
      parallelAdd(DG[it.index()+off],it.value());
  }
  template <int R,int C>
  static inline void parallelAdd(Eigen::Matrix<T,R,C>& DH,int offr,int offc,const SMatT& val) {
    for(int k=0; k<val.outerSize(); ++k)
      for(typename SMatT::InnerIterator it(val,k); it; ++it)
        parallelAdd(DH(it.row()+offr,it.col()+offc),it.value());
  }
  static inline Mat3X4T computeDTG(const Vec3T& G,const Vec3T& vLocal) {
    Mat3X4T DTG;
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        DTG(r,c)=G[r]*(c<3?vLocal[c]:1);
    return DTG;
  }
};
template <typename T>
class TrajectorySIPEnergy : public SIPEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  REUSE_MAP_FUNCS_T(SIPEnergy<T>)
  using typename SIPEnergy<T>::STrip;
  using typename SIPEnergy<T>::STrips;
  using typename SIPEnergy<T>::SMatT;
  using typename SIPEnergy<T>::EFunc;
  using typename SIPEnergy<T>::GFunc;
  using typename SIPEnergy<T>::HFunc;
  TrajectorySIPEnergy(const ArticulatedBody& body,const Vec& controlPoints,const ThetaTrajectory<T>& tt,T coef=1);
  void resetCoef(T newCoef);
  virtual bool debug(T perturbRange);
 protected:
  const ArticulatedBody& _body;
  const Vec& _controlPoints;
  const ThetaTrajectory<T>& _thetaTrajectory;
  T _coef;
};
}
#endif
