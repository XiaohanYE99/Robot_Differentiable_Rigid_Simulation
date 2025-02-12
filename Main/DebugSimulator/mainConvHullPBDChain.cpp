#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_EXP,3,0.8f,.1f,D2R(90),D2R(-75),D2R(0),0,3,0,0,0);

  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-1.2)));

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);

  utils.assemble(*(pt.RootElement()));
  //simulator
  ConvHullPBDSimulator sim(0.01f);//0.01
  //PBDSimulator sim(0.01f);
  sim.setOutput(true);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(ConvHullPBDSimulator::Vec3T(0,0,-9.81f));
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
