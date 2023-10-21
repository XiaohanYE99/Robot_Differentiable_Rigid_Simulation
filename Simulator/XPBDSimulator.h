#ifndef XPBD_SIMULATOR_H
#define XPBD_SIMULATOR_H

#include "PBDSimulator.h"

namespace PHYSICSMOTION {
//This solver implements a modified version of XPBD:
//XPBD: Position-Based Simulation of Compliant Constrained Dynamics
//The only difference is that we linearize the entire dynamic system and then use (linear) Gauss-Seidel to solve
class XPBDSimulator : public PBDSimulator {
 public:
  XPBDSimulator(T dt);
  void setArticulatedBody(std::shared_ptr<ArticulatedBody> body) override;
  void debugEnergy(T scale);
 protected:
  void update(const GradInfo& newPos,GradInfo& newPos2,const Vec& DE,const MatT& DDE,T alpha) const override;
  void schurComplement(Vec& DE,MatT& DDE) const;
  void updateLambda(const GradInfo& newPos,Vec& lambdaLimit,VecCM dx);
  void GaussSeidel(const GradInfo& newPos,Vec& dx,Vec& r,const Vec& DE,MatT DDE,T alpha,T tol,int maxIter,int debug) const;
  T energy(const GradInfo& newPos,Vec& DE,MatT& DDE,bool updateTangentBound=false) override;
  T dynamics(const GradInfo& newPos,MatTM M,VecM g);
  T jointLimitConstraint(const GradInfo& newPos,Vec& lambdaLimit,MatTM H,VecM g,VecCM dx,int& off) const;
  T normalConstraint(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,MatTM H,VecM g,VecCM dx,int& off) const;
  T tangentConstraint(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,MatTM H,VecM g,VecCM dx,int& off,bool updateTangentBound) const;
  T dragConstraint(const GradInfo& newPos,DragEnergy& drag,MatTM H,VecM g,VecCM dx,int& off) const;
  Vec _lambdaLimit;
};
}

#endif
