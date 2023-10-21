#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/PBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_EXP,10,0.5f,0.1f,D2R(30),0,0,0,3,0,0,0);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  //simulator
  PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(Vec3T(0,0,-9.81f));
  //setup kinematic
  sim.getJointPhysicsParameter(0)._isKinematic=true;
  sim.getJointPhysicsParameter(0)._kin=[&](T t,int nrDOF) {
    ASSERT_MSG(nrDOF==3,"Incorrect DOF!")
    Vec ret=Vec::Zero(3);
    ret[0]=cos(t*2.f);
    ret[1]=sin(t*2.f);
    return ret;
  };
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
