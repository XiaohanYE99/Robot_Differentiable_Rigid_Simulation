#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Simulator/PBDSimulator.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  std::vector<ArticulatedBody> bodies(2);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_XYZ,1,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[0]);
    utils2.assemble(*(pt.RootElement()));
  }
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_XYZ,10,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[1]);
    utils2.assemble(*(pt.RootElement()));
  }
  ArticulatedUtils(*body).combine(bodies);
  CTR(body->joint(11)._trans)=Vec3T(0,0,1).template cast<ArticulatedBody::T>();

  //debug Jacobian derivative
  CollisionGradInfo<T>::debugJacobianDeriv(*body);
  return 0;
}
