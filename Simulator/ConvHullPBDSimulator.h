#ifndef CONV_HULL_PBD_SIMULATOR_H
#define CONV_HULL_PBD_SIMULATOR_H

#include "Simulator.h"
#include <ConvexHull/Barrier.h>
#include <ConvexHull/CollisionGradInfo.h>
#include <ConvexHull/CustomPBDEnergy.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <ConvexHull/ConvexHullDistanceConvexEnergy.h>
#include <ConvexHull/ConvexHullDistanceFrictionEnergy.h>

namespace PHYSICSMOTION {
//This solver implements PBAD:
//The dynamics of articulated body formulation follows: Position-Based Time-Integrator for Frictional Articulated Body Dynamics
//The normal contact force formulation follows:         Simulation of Non-penetrating Elastic Bodies Using Distance Fields
//The frictional force formulation follows:             Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics
class PBDMatrixSolver;
class ConvHullPBDSimulator : public Simulator {
 public:
  typedef CollisionGradInfo<T> GradInfo;
  typedef GJKPolytope<T> GJ;
  typedef Logx Barrier;
  DECL_MAP_FUNCS
  ConvHullPBDSimulator(T dt);
  virtual void clearShape() override;
  virtual void addShape(std::shared_ptr<ShapeExact> shape) override;
  void setCustomEnergy(std::shared_ptr<CustomPBDEnergy<T>> custom);
  Vec3T getCentre(int k);
  MatT getdPTarget();
  MatT getdDTarget();
  MatT getdL();
  MatT getdLL();
  MatT getdPos();
  MatT getdD();
  Vec getConvexPoints();
  Vec getGlobalPoints(bool backward=true);
  void getJointPosGrad(GradInfo& grad,int k);
  void reset();
  void resetWithPos(Vec pos);
  void setArticulatedBody(std::shared_ptr<ArticulatedBody> body) override;
  void step() override;
  void stepwithbackward(bool debug);
  void backward(GradInfo& grad,bool debug=false);
  void backward();
  void detectCurrentContact() override;
  T gTol() const;
  void setGTol(T gTol);
  T x0() const;
  void setX0(T x0);
  bool hardLimit() const;
  void setHardLimit(bool hardLimit);
  Vec pos() const override;
  void setPos(const Vec& pos) override;
  Vec vel() const override;
  void setVel(const Vec& vel) override;
  void setclass(Eigen::Matrix<int,-1,1> c);
  void setCoefBarrier(T coefBarrier);
  void debugEnergy(T scale,const T* customDelta=NULL);
  void debugBackward(T scale,const T* customDelta=NULL);
  void debugBVHEnergy(GradInfo& grad);
  void setIsDesign(const std::vector<char>& isDesign,int startJID=0);
  void setPD(const Vec& P,const Vec& D,int st=0,int en=0,T* kp=NULL,T* kd=NULL);
  void updateConvexPoints(const Vec& X);
 protected:
  void detectContact(const Mat3XT& t) override;
  virtual bool detectLastContact();
  virtual void update(const GradInfo& newPos,GradInfo& newPos2,Vec& D,const Vec& DE,const MatT& DDE,T alpha) const;
  virtual T energy(GradInfo& grad,Vec* DE);
  virtual T normalEnergy(GradInfo& grad,Vec* DE,bool backward=false);
  virtual T tangentEnergy(GradInfo& grad,Vec* DE,bool backward=false);
  //data
  GradInfo _pos,_lastPos;
  T _gTol,_alpha,_epsV,_coefBarrier,_d0=1e-3f;
  Barrier _barrier;
  bool _hardLimit;
  std::shared_ptr<CustomPBDEnergy<T>> _custom;
  std::shared_ptr<PBDMatrixSolver> _sol;
  std::vector<GJKPolytope<T>> _obs;
  int _maxIt;
  //temporary variable
  Mat3XT _MRR,_MRt,_MtR,_Mtt;
  Mat3XT _MRRL,_MRtL,_MtRL,_MttL;
  std::vector<ContactManifold> _manifoldsLast;
};
}

#endif
