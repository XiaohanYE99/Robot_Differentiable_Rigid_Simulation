#ifndef ENTITY_ID_H
#define ENTITY_ID_H

#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/Utils.h>

namespace PHYSICSMOTION {
template <typename T>
class GJKPolytope;
template <typename T>
struct CollisionGradInfo {
  DECL_MAT_VEC_MAP_TYPES_T
  CollisionGradInfo();
  CollisionGradInfo(const ArticulatedBody& body,const Vec& theta);
  void reset(const ArticulatedBody& body,const Vec& theta);
  //data
  MatT _HTheta;
  Mat3XT _DTG;
  std::vector<Mat3XT> _globalVss;
  std::vector<std::shared_ptr<GJKPolytope<T>>> _polytopes;
  PBDArticulatedGradientInfo<T> _info;
};
}
#endif
