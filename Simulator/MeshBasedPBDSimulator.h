#ifndef MESHBASED_SIMULATOR_H
#define MESHBASED_SIMULATOR_H

#include "ConvHullPBDSimulator.h"
#include <ConvexHull/ConvexHullMeshDistanceFrictionEnergy.h>

namespace PHYSICSMOTION {
class MeshBasedPBDSimulator : public ConvHullPBDSimulator {
 public:
  MeshBasedPBDSimulator(T dt);
  void step() override;
 protected:
  void detectContact(const Mat3XT& t) override;
  T normalEnergy(GradInfo& grad,Vec* DE,bool backward=false) override;
  T tangentEnergy(GradInfo& grad,Vec* DE,bool backward=false) override;
  std::map<ContactGenerator::ContactManifold,std::vector<CCBarrierMeshFrictionEnergy<T,Barrier>::FrictionTerm>> _termMap;
};
}

#endif
