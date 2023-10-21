#ifndef TARGET_ENERGY_H
#define TARGET_ENERGY_H

#include "TrajectorySIPEnergy.h"

namespace PHYSICSMOTION {
template <typename T>
class TargetEnergy : public TrajectorySIPEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  REUSE_MAP_FUNCS_T(TrajectorySIPEnergy<T>)
  using typename TrajectorySIPEnergy<T>::STrip;
  using typename TrajectorySIPEnergy<T>::STrips;
  using typename TrajectorySIPEnergy<T>::SMatT;
  using typename TrajectorySIPEnergy<T>::EFunc;
  using typename TrajectorySIPEnergy<T>::GFunc;
  using typename TrajectorySIPEnergy<T>::HFunc;
  using TrajectorySIPEnergy<T>::debug;
  using TrajectorySIPEnergy<T>::computeDTG;
  using TrajectorySIPEnergy<T>::_body;
  using TrajectorySIPEnergy<T>::_controlPoints;
  using TrajectorySIPEnergy<T>::_thetaTrajectory;
  using TrajectorySIPEnergy<T>::_coef;
  TargetEnergy(const ArticulatedBody& body,
               const Vec& controlPoints,const ThetaTrajectory<T>& tt,
               int JID,const Vec3T& tar,const Vec3T& dir,const Vec3T& ori,
               T coefP=1,T coefD=1,bool JTJApprox=false);
  TargetEnergy(const TargetEnergy& ref,
               const Vec& controlPoints,const ThetaTrajectory<T>& tt);
  bool eval(EFunc* E,GFunc* G,HFunc* H) override;
  const Vec3T& globalPos() const;
  const Vec3T& localPos() const;
  const Vec3T& localDir() const;
  int JID() const;
 private:
  void eval(const Vec3T& local,const Vec3T& global,T coef,EFunc* E,GFunc* G,HFunc* H) const;
  Vec3T _localPos,_globalPos;
  Vec3T _localDir,_globalDir;
  bool _JTJApprox;
  int _JID;
  T _coefD;
};
}
#endif
