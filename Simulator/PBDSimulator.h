#ifndef PBD_SIMULATOR_H
#define PBD_SIMULATOR_H

#include "Simulator.h"
#include <Articulated/PBDArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
//This solver implements PBAD:
//The dynamics of articulated body formulation follows: Position-Based Time-Integrator for Frictional Articulated Body Dynamics
//The normal contact force formulation follows:         Simulation of Non-penetrating Elastic Bodies Using Distance Fields
//The frictional force formulation follows:             Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics
class PBDMatrixSolver;
class PBDSimulator : public Simulator {
 public:
  DECL_MAP_FUNCS
  struct SolverState {
    SolverState();
    SolverState(const ArticulatedBody& body,const Vec& x);
    void reset(const ArticulatedBody& body,const Vec& x);
    //data
    GradInfo _pos;
    Vec _DE;
    //Temporary variable for CRBA/Native
    MatT _DDE;
    //Temporary variable for ABA does not use the assembled matrix.
    //Rather, it uses these temporary variables to assemble online.
    Mat3XT _MRR,_MRt,_MtR,_Mtt;
    MatT _diag;
  };
  PBDSimulator(T dt);
  void setArticulatedBody(std::shared_ptr<ArticulatedBody> body) override;
  void step() override;
  Vec pos() const override;
  void setPos(const Vec& pos) override;
  Vec vel() const override;
  void setVel(const Vec& vel) override;
  void detectCurrentContact() override;
  void debugEnergy(T scale);
  void setJTJ(bool JTJ);
  virtual void setCrossTerm(bool cross,bool CRBA=false);
 protected:
  void solveBody(SolverState& state,SolverState& state2);
  virtual void update(const SolverState& state,SolverState& state2,T alpha);
  void computeLocalContactPos(const Mat3XT& t) override;
  virtual T energy(SolverState& state,bool updateTangentBound=false);
  T contactEnergy(const ContactManifold& m,ContactPoint& p,SolverState& state,bool updateTangentBound);
  T normalEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H) const;
  T tangentEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H,bool updateTangentBound) const;
  //data
  GradInfo _pos,_lastPos;
  T _gTol,_epsV;
  bool _JTJ,_crossTerm;
  std::shared_ptr<PBDMatrixSolver> _sol;
  int _maxIt;
  //temporary variable for energy assembly
  Mat3XT _GB,_G;
  MatT _DDER;
};
}

#endif
