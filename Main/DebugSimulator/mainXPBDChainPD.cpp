#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/XPBDSimulator.h>
#include <Simulator/SoftJoint.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_XYZ,10,0.5f,0.1f,D2R(90),0,0,0,3,0,0,0);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  //create floor
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-3)));
  //simulator
  XPBDSimulator sim(0.033f);
  sim.setArticulatedBody(body);
  sim.setGravity(XPBDSimulator::Vec3T(0,0,-9.81f));
  sim.setHeuristcGuessStiffness();
  sim.addShape(floor);
  sim.setOutput(true);
  //PD controller
  for(int k=0; k<body->nrJ(); k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::cout<<J.nrDOF()<<std::endl;
    //param
    auto& param=sim.getJointPhysicsParameter(k);
    param._kp=1e3;
    param._kd=1e3;
    param._tarP=[&](T time,int)->Vec {
      return Vec::Ones(J.nrDOF())*(sin(time)+cos(time))*0.1f;
    };
    param._tarD=[&](T time,int)->Vec {
      return Vec::Ones(J.nrDOF())*(cos(time)+sin(time))*0.1f;
    };
  }
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
