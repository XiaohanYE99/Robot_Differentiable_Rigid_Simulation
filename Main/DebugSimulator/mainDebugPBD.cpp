#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Simulator/PBDSimulator.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),0,10,0.5f,0.1f,D2R(10),0,0,0,3,0,0,0);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-1)));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  
  utils.assemble(*(pt.RootElement()));
  {
    //PBDSimulator
    ConvHullPBDSimulator sim(0.01f);
    sim.setArticulatedBody(body);
    sim.setGravity(Vec3T(0,0,-9.81f));
    sim.addShape(floor);
    sim.setJTJ(false);
    sim.setCrossTerm(true);
    sim.debugEnergy(10.7);
  }
  {
    //XPBDSimulator
    /*PBDSimulator sim(0.01f);
    sim.setArticulatedBody(body);
    sim.setGravity(Vec3T(0,0,-9.81f));
    sim.setJTJ(true);
    sim.setCrossTerm(true);
    sim.debugEnergy(1);*/
  }
  return 0;
}
