#ifndef SMOOTHNESS_ENERGY_H
#define SMOOTHNESS_ENERGY_H

#include "TrajectorySIPEnergy.h"

namespace PHYSICSMOTION {
template <typename T>
class SmoothnessEnergy : public TrajectorySIPEnergy<T> {
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
  SmoothnessEnergy(const ArticulatedBody& body,const Vec& controlPoints,const ThetaTrajectory<T>& tt,T coef=1,int d=1);
  bool eval(EFunc* E,GFunc* G,HFunc* H) override;
 private:
  SMatT _hessian;
};
}
#endif
