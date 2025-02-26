#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/MeshBasedPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  std::vector<ArticulatedBody> bodies(2);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-15,-15,-8),BBoxExact::Vec3T(30,30,-1)));
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_EXP,5,.5f,0.1f,D2R(60),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[0]);
    utils2.assemble(*(pt.RootElement()));
  }
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createChain(*(pt.RootElement()),Joint::ROT_3D_EXP,5,.5f,0.1f,D2R(60),0,0,0,3,0,0,0);
    ArticulatedUtils utils2(bodies[1]);
    utils2.assemble(*(pt.RootElement()));
  }
  ArticulatedUtils(*body).combine(bodies);
  CTR(body->joint(7)._trans)=Vec3T(0,0,.5).template cast<ArticulatedBody::T>();
  //simulator
  MeshBasedPBDSimulator sim(0.01f);
  //ConvHullPBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(Vec3T(0,0,-9.81f));
  sim.setOutput(true);
  //visualizer
  SimulatorVisualizer render(sim);
  render.setLightSize(0);
  std::vector<Eigen::Matrix<float,3,1>>ArticulateDiffuse;
  Eigen::Matrix<float,3,1> LightDiffuse=Eigen::Matrix<float,3,1>(.5,.5,.5);
  render.setLightDiffuse(LightDiffuse);
  //setup kinematic
  /*for(int k=1; k<2; k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::cout<<k<<" "<<J.nrDOF()<<std::endl;
    //param
    auto& param=sim.getJointPhysicsParameter(k);
    param._kp=1e3;
    param._kd=1e1;
    param._tarP=[&](T time,int)->Vec {
      return Vec3T(cos(time*2),sin(time*2),0)*4.5;
    };
    param._tarD=[&](T time,int)->Vec {
      return Vec3T(-sin(time*2),cos(time*2),0)*9;
    };
  }
  for(int k=12; k<13; k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::cout<<k<<" "<<J.nrDOF()<<std::endl;
    //param
    auto& param=sim.getJointPhysicsParameter(k);
    param._kp=1e3;
    param._kd=1e1;
    param._tarP=[&](T time,int)->Vec {
      return Vec3T(cos(time*2),sin(time*2),0)*3;
    };
    param._tarD=[&](T time,int)->Vec {
      return Vec3T(-sin(time*2),cos(time*2),0)*6;
    };
  }*/
  visualizeSimulator(argc,argv,sim);
  //run app
  /*std::vector<Vec> traj;
  traj.resize(10);
  for(int i=0;i<1;i++) {
    std::cout<<i<<std::endl;
    sim.step();
    traj.push_back(sim.pos());
  }
  int frame=0;
  render.visualize(&traj,&frame);*/
  return 0;
}
