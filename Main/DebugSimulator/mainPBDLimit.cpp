#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/XPBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  double stiff=1e4f;
  std::vector<ArticulatedBody> bodies;
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  for(int i=0; i<4; i++) {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createBox(*(pt.RootElement()),.5f,.5f,.5f,0);
    bodies.push_back(ArticulatedBody());
    ArticulatedUtils utils(bodies.back());
    utils.assemble(*(pt.RootElement()));
    utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>());
    utils.transformTorso(ArticulatedBody::Vec3T(ArticulatedBody::Vec3T::UnitX()*2*i));
    if(i<3) {
      //linear limit
      Joint::Mat3XT& lmt=bodies.back().joint(0)._limits;
      lmt.setZero(3,3);
      lmt(0,i)=-1;
      lmt(1,i)=1;
      lmt.row(2).setConstant(stiff);
      Joint::Mat3XT& almt=bodies.back().joint(1)._limits;
      almt.setZero(3,3);
      almt.row(2).setConstant(stiff);
    } else if(i==3) {
      //ball limit
      Joint::Mat3XT& lmt=bodies.back().joint(0)._limits;
      lmt.setConstant(3,3,std::numeric_limits<double>::infinity());
      lmt.row(1).setConstant(1);
      lmt.row(2).setConstant(stiff);
      Joint::Mat3XT& almt=bodies.back().joint(1)._limits;
      almt.setZero(3,3);
      almt.row(2).setConstant(stiff);
    }
  }
  for(int i=0; i<3; i++) {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createBox(*(pt.RootElement()),.5f,.5f,.5f,0);
    bodies.push_back(ArticulatedBody());
    ArticulatedUtils utils(bodies.back());
    utils.assemble(*(pt.RootElement()));
    utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>());
    utils.transformTorso(ArticulatedBody::Vec3T(ArticulatedBody::Vec3T::UnitX()*(2*i+8)));
    if(i<3) {
      //linear angular limit
      Joint::Mat3XT& lmt=bodies.back().joint(0)._limits;
      lmt.setZero(3,3);
      lmt.row(2).setConstant(stiff);
      Joint::Mat3XT& almt=bodies.back().joint(1)._limits;
      almt.setZero(3,3);
      almt(0,i)=-1;
      almt(1,i)=1;
      almt.row(2).setConstant(stiff);
    }
  }
  for(int i=0; i<3; i++) {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    ArticulatedLoader::createBox(*(pt.RootElement()),.5f,.5f,.5f,0);
    bodies.push_back(ArticulatedBody());
    ArticulatedUtils utils(bodies.back());
    utils.assemble(*(pt.RootElement()));
    utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>(),true);
    utils.transformTorso(ArticulatedBody::Vec3T(ArticulatedBody::Vec3T::UnitX()*(2*i+14)));
    {
      Joint::Mat3XT& lmt=bodies.back().joint(0)._limits;
      lmt.setZero(3,3);
      lmt.row(2).setConstant(stiff);
    }
    if(i<2) {
      //dual cone limit
      Joint::Mat3XT& almt=bodies.back().joint(1)._limits;
      almt.setConstant(3,3,std::numeric_limits<double>::infinity());
      almt.col(0).setZero();
      almt(2,0)=stiff;
      almt(0,i+1)=-M_PI/18;
      almt(1,i+1)=M_PI/18;
      almt(2,i+1)=stiff;
    } else {
      //ellipse limit
      Joint::Mat3XT& almt=bodies.back().joint(1)._limits;
      almt.setConstant(3,3,std::numeric_limits<double>::infinity());
      almt.col(0).setZero();
      almt(2,0)=stiff;
      almt(1,1)=M_PI/18;
      almt(2,1)=stiff;
      almt(1,2)=M_PI/18;
      almt(2,2)=stiff;
    }
  }
  ArticulatedUtils(*body).combine(bodies);
  //simulator
  PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.setGravity(Vec3T(0,0,-9.81f));
  sim.setHeuristcGuessStiffness();
  sim.setJTJ(false);
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
