#ifndef SOFT_JOINT_H
#define SOFT_JOINT_H

#include "Simulator.h"
#include <Articulated/PBDArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
class SoftJoint {
 public:
  typedef FLOAT T;
  typedef PBDArticulatedGradientInfo<T> GradInfo;
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  SoftJoint();
  Vec3T posA(const GradInfo& pos) const;
  Vec3T posB(const GradInfo& pos) const;
  Vec3T rotA(const GradInfo& pos,int d) const;
  Vec3T rotB(const GradInfo& pos,int d) const;
  Vec3T rotAL(int d) const;
  Vec3T rotBL(int d) const;
  void setRandom(const ArticulatedBody& body);
  //energy model
  T energyLinear(const ArticulatedBody& body,const GradInfo& newPos,int d,Vec* DE,MatT& DDE,Mat3XT& G,
                 Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm) const;
  T energyBall(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,Mat3XT& G,
               Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm,T tol) const;
  T energyAngular(const ArticulatedBody& body,const GradInfo& newPos,int d,Vec* DE,MatT& DDE,Mat3XT& G,
                  Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm) const;
  T energy(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,Mat3XT& G,
           Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm,T tol=1e-6f) const;
  T energy(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,bool JTJ,bool crossTerm,T tol=1e-6f) const;
  //constraint model
  void constraintLinear(const ArticulatedBody& body,const GradInfo& newPos,int d,int nrDOF,
                        std::vector<Simulator::XPBDConstraint>& Css,int& off) const;
  void constraintBall(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,
                      std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol) const;
  void constraintAngular(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,int d,
                         std::vector<Simulator::XPBDConstraint>& Css,int& off) const;
  void constraint(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,
                  std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol=1e-6f) const;
  //debug
  void debug(const ArticulatedBody& body);
  void debugInner(const std::string& name,const ArticulatedBody& body);
 public:
  Mat3T _linearLimit;
  Mat3T _angularLimit;
  Mat3X4T _transA,_transB;
  int _jidA=-1,_jidB=-1;
};
}

#endif
