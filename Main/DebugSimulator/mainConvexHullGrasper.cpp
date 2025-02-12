#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<ArticulatedBody> bodies(2);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createGrasper(*(pt.RootElement()),0.1);
    ArticulatedUtils utils2(bodies[0]);
    utils2.assemble(*(pt.RootElement()));
  }
  {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createBall(*(pt.RootElement()),Joint::TRANS_3D|Joint::ROT_3D_EXP,0.1);
    ArticulatedUtils utils2(bodies[1]);
    utils2.assemble(*(pt.RootElement()));
  }

  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-1.0)));
  ArticulatedUtils(*body).combine(bodies);
  //simulator
  ConvHullPBDSimulator sim(0.01f);//0.01
  //PBDSimulator sim(0.01f);
  sim.setOutput(true);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(ConvHullPBDSimulator::Vec3T(0,0,-9.81f));
  Vec pos=Vec::Zero(body->nrDOF());
  /*pos[3]=D2R(300);
  pos[4]=D2R(120);
  pos[5]=D2R(240);
  pos[6]=D2R(120);
  sim.resetWithPos(pos);*/
  //PD controller
  for(int k=3; k<body->nrJ()-2; k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::cout<<J.nrDOF()<<" "<<J._class<<std::endl;
    //param
    auto& param=sim.getJointPhysicsParameter(k);
    param._kp=1e1;
    param._kd=1e0;
    if(k==3){
        param._tarP=[&](T time,int)->Vec {
        //std::cout<<pos[k+1]<<std::endl;
        return Vec::Ones(J.nrDOF())*((sin(time*10))*.4f);
        };
        param._tarD=[&](T time,int)->Vec {
        return Vec::Ones(J.nrDOF())*((cos(time*10))*4.f);
        };
    }
    if(k==4){
        param._tarP=[&](T time,int)->Vec {
        //std::cout<<pos[k+1]<<std::endl;
        return Vec::Ones(J.nrDOF())*((sin(time*10))*-.4f);
        };
        param._tarD=[&](T time,int)->Vec {
        return Vec::Ones(J.nrDOF())*((cos(time*10))*-4.f);
        };
    }
  }
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
