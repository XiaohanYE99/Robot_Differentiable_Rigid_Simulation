#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Simulator/PBDSimulator.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  mpfr_float::default_precision(256);
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  std::vector<ArticulatedBody> bodies(2);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_EXP,10,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[0]);
    utils2.assemble(*(pt.RootElement()));
  }
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_EXP,10,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[1]);
    utils2.assemble(*(pt.RootElement()));
  }
  ArticulatedUtils(*body).combine(bodies);
  CTR(body->joint(12)._trans)=Vec3T(0,1,0).template cast<ArticulatedBody::T>();
  //std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-1,-1,-1),BBoxExact::Vec3T(1,1,-1)));
  //std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  //ArticulatedUtils utils(*body);

  //utils.assemble(*(pt.RootElement()));
  
    //PBDSimulator
    ConvHullPBDSimulator sim(0.01f);
    sim.setArticulatedBody(body);
    sim.setGravity(Vec3T(0,0,-9.81f));
    sim.addShape(floor);
    sim.setJTJ(false);
    sim.setCrossTerm(true);
    sim.debugEnergy(18.7);
  
  
    //XPBDSimulator
    /*PBDSimulator sim(0.01f);
    sim.setArticulatedBody(body);
    sim.setGravity(Vec3T(0,0,-9.81f));
    sim.setJTJ(true);
    sim.setCrossTerm(true);
    sim.debugEnergy(1);*/
  
  return 0;
}
