#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/PBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_EXP,10,0.5f,0.1f,D2R(30),0,0,0,3,0,0,0);

  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  /*std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",true,false)));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(*bodyTmp));
  ArticulatedUtils(*body).fix([&](int jid,const Joint&)->bool {
    return jid==7;
  },ArticulatedBody::Vec::Zero(body->nrDOF()));
  ArticulatedUtils(*body).simplify(ArticulatedUtils::Vec::Zero(body->nrDOF()),10);
  ArticulatedUtils(*body).tessellate(true);*/
  utils.assemble(*(pt.RootElement()));
  //utils.convexDecompose();
  //simulator
  ConvHullPBDSimulator sim(0.005f);//0.01
  //PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(ConvHullPBDSimulator::Vec3T(0,0,-9.81f));
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
