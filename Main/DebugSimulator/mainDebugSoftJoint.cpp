#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SoftJoint.h>

using namespace PHYSICSMOTION;

int main() {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),0,10,0.5f,0.1f,D2R(30),0,0,0,3,0,0,0);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  //debug
  SoftJoint().debug(*body);
  return 0;
}
