#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/PBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  std::vector<ArticulatedBody> bodies(2);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_EXP,10,.5f,0.24f,D2R(30),0,0,0,3,0,0,0);
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
  CTR(body->joint(12)._trans)=Vec3T(0,0,1).template cast<ArticulatedBody::T>();
  //simulator
  ConvHullPBDSimulator sim(0.005f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(Vec3T(0,0,-9.81f));
  //setup kinematic
  /*sim.getJointPhysicsParameter(1)._isKinematic=true;
  sim.getJointPhysicsParameter(1)._kin=[&](T t,int nrDOF) {
    ASSERT_MSG(nrDOF==3,"Incorrect DOF!")
    Vec ret=Vec::Zero(3);
    ret[0]=cos(t*2.f);
    ret[1]=sin(t*2.f);
    return ret;
  };*/
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}