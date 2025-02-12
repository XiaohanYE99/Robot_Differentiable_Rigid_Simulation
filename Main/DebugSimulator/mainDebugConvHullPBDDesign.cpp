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
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_XYZ,10,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
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

  //PBDSimulator
  ConvHullPBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.setGravity(Vec3T(0,0,-9.81f));
  sim.addShape(floor);
  sim.setGTol(1e-40f);
  sim.setX0(0.9f);
  sim.setOutput(true);

  std::shared_ptr<CustomPBDEnergy<T>> custom(new CustomPBDEnergy<T>(body->nrDOF(),false));
  sim.setCustomEnergy(custom);

  std::vector<char> isDesign(body->nrDOF(),false);
  isDesign[5]=true;
  isDesign[7]=true;
  sim.setIsDesign(isDesign);

  T customDelta=1e-8f;
  sim.debugBackward(0.1,&customDelta);
  sim.debugEnergy(0.01);
  return 0;
}
