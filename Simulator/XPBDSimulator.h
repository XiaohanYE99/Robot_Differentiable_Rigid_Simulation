#ifndef XPBD_SIMULATOR_H
#define XPBD_SIMULATOR_H

#include "PBDSimulator.h"

namespace PHYSICSMOTION {
//This solver implements a modified version of XPBD:
//XPBD: Position-Based Simulation of Compliant Constrained Dynamics
class XPBDSimulator : public PBDSimulator {
 public:
  XPBDSimulator(T dt);
  void step() override;
  void debugEnergy(T scale);
  void setCrossTerm(bool cross,bool CRBA);
 protected:
  T energy(SolverState& state,bool updateTangentBound=false) override;
  void projectBody(SolverState& state);
  void contactNormal(const GradInfo& newPos,const ContactManifold& m,const ContactPoint& p,int& off);
  void contactTangent(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,int d,int& off);
  template <int N=3>
  void project(GradInfo& newPos,int beg,int end,T* limit,bool updateKinematics);
  std::vector<XPBDConstraint> _Css;
  T _epsDenom;
};
}

#endif
