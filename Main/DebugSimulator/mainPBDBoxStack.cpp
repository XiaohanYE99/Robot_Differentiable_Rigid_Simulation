#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/PBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  std::vector<ArticulatedBody> bodies(10);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  for(int d=0; d<10; d++) {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createBox(*(pt.RootElement()),.5f,.5f,.5f,0);
    ArticulatedUtils utils(bodies[d]);
    utils.assemble(*(pt.RootElement()));
    utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>());
  }
  ArticulatedUtils(*body).combine(bodies);
  for(int d=0; d<10; d++)
    CTR(body->joint(1+d*4)._trans)=Vec3T(0,0,d).template cast<ArticulatedBody::T>();
  //simulator
  PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(Vec3T(0,0,-9.81f));
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
