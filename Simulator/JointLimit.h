#ifndef JOINT_LIMIT_H
#define JOINT_LIMIT_H

#include "Simulator.h"

namespace PHYSICSMOTION {
class JointLimit {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  //energy model
  static T energyExp(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE,T tol=1e-6f);
  static T energySimple(VecCM x,const Joint& J,int c,int nrD,Vec* DE,MatT* DDE);
  template <typename Barrier>
  static bool energySimple(VecCM x,const Joint& J,int c,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p);
  static T energyBall(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE,T tol=1e-6f);
  template <typename Barrier>
  static bool energyBall(VecCM x,const Joint& J,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p,T tol=1e-6f);
  static T energy(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE);
  template <typename Barrier>
  static bool energy(VecCM x,const Joint& J,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p);
  //constraint model
  static void constraintExp(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol=1e-6f);
  static void constraintSimple(VecCM x,const Joint& J,int c,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off);
  static void constraintBall(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol=1e-6f);
  static void constraint(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off);
  //debug
  static void debugSwing(T tol=1e-6f);
  static void debugExp(int off=3,T tol=1e-6f);
  static void debugSimple(int off=3);
  template <typename Barrier>
  static void debugSimple(int off,const Barrier& p);
  static void debugBall(int off=3);
  template <typename Barrier>
  static void debugBall(int off,const Barrier& p);
  static void debug(const std::string& name,const Vec& x,const Joint& J,std::function<bool(const Vec&)> func);
  template <typename Barrier>
  static void debugHard(const std::string& name,const Vec& x,const Joint& J,const Barrier& p);
 private:
  static T getASinY(T y,T z,Vec2T* DY,Mat2T* DDY,T tol);
  static T getASinZ(T y,T z,Vec2T* DZ,Mat2T* DDZ,T tol);
  static T getSwingAngleY(T y, T z,Vec2T* DY,Mat2T* DDY,T tol);
  static T getSwingAngleZ(T y, T z,Vec2T* DZ,Mat2T* DDZ,T tol);
};
}

#endif
