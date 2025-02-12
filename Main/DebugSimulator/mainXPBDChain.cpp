#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/XPBDSimulator.h>
#include <Simulator/SoftJoint.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),0,10,0.5f,0.1f,D2R(30),0,0,0,3,0,0,0);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  //create floor
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  //simulator
  XPBDSimulator sim(0.033f);
  sim.setArticulatedBody(body);
  sim.setGravity(XPBDSimulator::Vec3T(0,0,-9.81f));
  sim.setHeuristcGuessStiffness();
  sim.addShape(floor);
  sim.setOutput(true);
  //close loop
  SoftJoint joint;
  joint._jidA=0;
  joint._jidB=body->nrJ()-1;
  joint._linearLimit.row(1).setZero();
  joint._linearLimit.row(2).setConstant(1);
  sim.addJoint(joint);
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
