#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedBody::T footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),Joint::TRANS_3D|Joint::HINGE_JOINT,0.2f,footLen,0.08f,D2R(10),D2R(60),D2R(60),true);
  //ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_2D|Joint::ROT_3D_EXP,3,0.8f,.1f,D2R(90),D2R(0),D2R(0),0,3,0,0,0);
  //std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-1.0)));
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-15,0.62,-5),BBoxExact::Vec3T(5,5,5)));

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);

  ArticulatedUtils utils(*body);

  utils.assemble(*(pt.RootElement()));
  std::vector<int> JointId={3,5};
  utils.convexDecompose(JointId,8);
  //simulator
  ConvHullPBDSimulator sim(0.01f);//0.01
  //PBDSimulator sim(0.01f);
  sim.setOutput(true);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(ConvHullPBDSimulator::Vec3T(0,9.81f,0));
  //PD controller
  for(int k=1; k<body->nrJ(); k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::cout<<k<<" "<<J.nrDOF()<<" "<<J._class<<std::endl;
    //param
    auto& param=sim.getJointPhysicsParameter(k);
    param._kp=1e1;
    param._kd=1e0;
    param._tarP=[&](T time,int)->Vec {
      return Vec::Ones(J.nrDOF())*((sin(time*10))*.2f);
    };
    param._tarD=[&](T time,int)->Vec {
      return Vec::Ones(J.nrDOF())*((cos(time*10))*2.f);
    };
  }
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
