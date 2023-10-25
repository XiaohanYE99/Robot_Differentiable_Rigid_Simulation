#ifndef COLLISION_BARRIER_ENERGY_H
#define COLLISION_BARRIER_ENERGY_H

#include "EntityId.h"
#include "DistanceEnergy.h"
#include "CCSeparatingPlanes.h"
#include "TrajectorySIPEnergy.h"
#include "ConvexHullDistanceEnergy.h"
#include "ConvexHullDistanceConvexEnergy.h"
#include <Articulated/PBDArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
template <typename T,typename PFunc>
class CollisionBarrierEnergy : public TrajectorySIPEnergy<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  REUSE_MAP_FUNCS_T(TrajectorySIPEnergy<T>)
  using typename TrajectorySIPEnergy<T>::STrip;
  using typename TrajectorySIPEnergy<T>::STrips;
  using typename TrajectorySIPEnergy<T>::SMatT;
  using typename TrajectorySIPEnergy<T>::EFunc;
  using typename TrajectorySIPEnergy<T>::GFunc;
  using typename TrajectorySIPEnergy<T>::HFunc;
  typedef CCBarrierConvexEnergy<T,PFunc> CCBarrierEnergyType;
  //typedef CCBarrierEnergy<T,PFunc> CCBarrierEnergyType;
  typedef CollisionGradInfo<T> GradInfo;
  using TrajectorySIPEnergy<T>::debug;
  using TrajectorySIPEnergy<T>::computeDTG;
  using TrajectorySIPEnergy<T>::_body;
  using TrajectorySIPEnergy<T>::_controlPoints;
  using TrajectorySIPEnergy<T>::_thetaTrajectory;
  using TrajectorySIPEnergy<T>::_coef;
  CollisionBarrierEnergy(const ArticulatedBody& body,
                         const Vec& controlPoints,const ThetaTrajectory<T>& tt,
                         std::vector<std::pair<EntityId<T>,EntityId<T>>>& TTPairs,CCSeparatingPlanes<T>& CCPlanes,
                         std::unordered_map<rational,GradInfo>& gradInfo,T d0,T x0,T coef=1,bool JTJApprox=false,bool implicit=true);
  const std::vector<std::pair<EntityId<T>,EntityId<T>>>& getTriPairs() const;
  const std::vector<std::pair<EntityId<T>,EntityId<T>>>& getFailedTT() const;
  const std::vector<std::pair<EntityId<T>,EntityId<T>>>& getFailedCC() const;
  const CCSeparatingPlanes<T>& getCCPlanes() const;
  bool eval(EFunc* E,GFunc* G,HFunc* H) override;
  void randomizeCCPairs(std::shared_ptr<MeshExact> obstacle);
  void randomizeTTPairs(std::shared_ptr<MeshExact> obstacle);
  void setConstructGradInfo(bool construct);
  void setOutput(bool output);
  bool implicit() const;
  bool updateCCPlanes();
 private:
  bool evalCC(const std::pair<EntityId<T>,EntityId<T>>& pair,CCSeparatingPlane<T>& plane,EFunc* E,bool DG,bool DH,GradInfo& gradInfo);
  bool evalTT(const std::pair<EntityId<T>,EntityId<T>>& pair,EFunc* E,bool DG,bool DH,GradInfo& gradInfo) const;
  std::vector<std::pair<EntityId<T>,EntityId<T>>>& _TTPairs;
  CCSeparatingPlanes<T>& _CCPlanes;
  std::unordered_map<rational,GradInfo>& _gradInfo;
  bool _JTJApprox,_constructGradInfo,_implicit,_output;
  PFunc _pTT,_pCC;
  T _d0;
  //temporary data
  std::vector<std::pair<EntityId<T>,EntityId<T>>> _CCPairs,_failedTT,_failedCC;
  std::vector<CCSeparatingPlane<T>*> _planeHandles;
};
}
#endif
