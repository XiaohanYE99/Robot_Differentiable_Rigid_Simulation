#ifndef ARTICULATED_UTILS_H
#define ARTICULATED_UTILS_H

#include "ArticulatedBody.h"
#include <Environment/CompositeShapeExact.h>
#include <tinyxml2.h>

namespace PHYSICSMOTION {
class ArticulatedUtils {
 public:
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  struct TransInfo {
    Mat3X4T _TLast,_X0;
    Vec3T _lenX,_lenY,_lenZ,_origin;
    T _rad,_rad2,_rad3,_rho;
  };
  ArticulatedUtils(ArticulatedBody& body);
  void assembleGlobalVars(const tinyxml2::XMLElement& pt);
  void assemble(const tinyxml2::XMLElement& pt);
  void assembleJoints(const tinyxml2::XMLElement& pt,std::vector<TransInfo>& infos);
  void assembleJoint(const tinyxml2::XMLElement& pt,int parent,std::vector<TransInfo>& infos);
  //articulated body topology operations
  static void mergeMesh(Joint& joint,
                        const std::vector<std::shared_ptr<ShapeExact>>& geomsMerged,
                        const std::vector<CompositeShapeExact::Mat3X4T>& transMerged);
  void addBase(int dim,const Vec3T& planeNormal);
  void combine(const std::vector<ArticulatedBody>& bodies);
  Vec mergeChildren(int jid,const Vec& DOF);
  Vec fix(std::function<bool(int,const Joint&)> canFix,const Vec& DOF);
  Vec eliminate(std::function<bool(int,const Joint&)> canEliminate,const Vec& DOF);
  Vec simplify(std::function<bool(int,const Joint&)> canSimplify,const Vec& DOF,int nrDebug);
  Vec simplify(const Vec& DOF,int nrDebug);
  Vec replaceJoint(const Vec& DOF,int jid,Mat3XT axes);
  void convexDecompose(T rho=1);
  void addBody(ArticulatedBody& body);
  //mesh operation
  void tessellate(bool rebuildBVH=false);
  void BBApproxiate(bool rebuildBVH=false);
  void makeConvex();
  //rigid transformation
  Mat3X4T transformTorso(const Mat3X4T& trans);
  Mat3X4T transformTorso(const Mat3T& R);
  Mat3X4T transformTorso(const Vec3T& t);
  Mat3X4T transformTorsoToCOM();
  Mat3X4T scaleBody(T coef);
  //mass transformation
  void scaleMass(T coef);
  T totalMass() const;
 private:
  ArticulatedBody& _body;
  Vec _canSimplify;
  //transient data, no included in read/write
  T _globalRho;
  T _globalLimitCoef;
  T _globalControlCoef;
  T _globalDampingCoef;
  T _globalSampleDist;
  T _globalSampleDistLocal;
};
}

#endif
