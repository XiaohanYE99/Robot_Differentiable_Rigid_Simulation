#ifndef CONVHULLPBD_SIMULATOR_H
#define CONVHULLPBD_SIMULATOR_H

#include "Simulator.h"
#include <SIPCollision/Barrier.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <SIPCollision/CollisionBarrierEnergy.h>

namespace PHYSICSMOTION {
//This solver implements PBAD:
//The dynamics of articulated body formulation follows: Position-Based Time-Integrator for Frictional Articulated Body Dynamics
//The normal contact force formulation follows:         Simulation of Non-penetrating Elastic Bodies Using Distance Fields
//The frictional force formulation follows:             Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics
class PBDMatrixSolver;
class ConvHullPBDSimulator : public Simulator {
 public:
  typedef PBDArticulatedGradientInfo<T> GradInfo;
  typedef GJKPolytope<T> GJ;
  typedef Px Barrier;
  DECL_MAP_FUNCS
  ConvHullPBDSimulator(T dt);
  void setArticulatedBody(std::shared_ptr<ArticulatedBody> body) override;
  void step() override;
  Vec pos() const override;
  void setPos(const Vec& pos) override;
  Vec vel() const override;
  void setVel(const Vec& vel) override;
  void setGravity(const Vec3T& g) override;
  void detectCurrentContact() override;
  void simpleCollisionHandle();
  void debugEnergy(T scale);
  void setOutput(bool output);
  void setJTJ(bool JTJ);
  void setCrossTerm(bool cross);
 protected:
  virtual void update(const GradInfo& newPos,GradInfo& newPos2,Vec& D,const Vec& DE,const MatT& DDE,T alpha) const;
  void computeLocalContactPos(const Mat3XT& t) override;
  void linesearch(const GradInfo& newPos,const T e,const Vec& DE,const MatT& DDE,T &alpha,const T alphaMin);
  void mask(Vec& diag,Vec& DE,MatT& DDE) const;
  virtual T energy(CollisionGradInfo<T>& grad,Vec& DE,MatT& DDE,bool updateTangentBound=false);
  T contactEnergy(const ContactManifold& m,ContactPoint& p,
                  const GradInfo& newPos,Vec& DE,MatT& DDE,Mat3XT& G,
                  Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool updateTangentBound) const;
  T normalEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H) const;
  T tangentEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H,bool updateTangentBound) const;
  //data
  GradInfo _pos,_lastPos,_subPos;

  T _gTol,_alpha,_epsV;
  bool _output,_JTJ,_crossTerm;
  std::shared_ptr<PBDMatrixSolver> _sol;
  Mat3XT _JRCF;
  int _maxIt;
  //temporary variable for ABA
  //ABA does not use the assembled matrix.
  //Rather, it uses these temporary variables to assemble online.
  Mat3XT _MRR,_MRt,_MtR,_Mtt;
  Vec _diag;
  Px _barrier;
  T _d0=1e-3;
  bool _ispenetrated;
};
}

#endif
