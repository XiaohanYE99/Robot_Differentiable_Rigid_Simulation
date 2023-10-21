#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/PBDSimulator.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //create body
  std::vector<ArticulatedBody> bodies(3);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  for(int d=0; d<3; d++) {
    tinyxml2::XMLDocument pt;
    pt.InsertEndChild(pt.NewElement("root"));
    if(d==2)
      ArticulatedLoader::createBox(*(pt.RootElement()),5.f,5.f,2.f,0);
    else ArticulatedLoader::createBox(*(pt.RootElement()),.5f,.5f,.5f,0);
    ArticulatedUtils utils(bodies[d]);
    utils.assemble(*(pt.RootElement()));
    utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>());
    if(d==2)
      utils.scaleMass(0.01f);
  }
  ArticulatedUtils(*body).combine(bodies);
  CTR(body->joint(1)._trans)=Vec3T(0,0,-1.5f).template cast<ArticulatedBody::T>();
  CTR(body->joint(5)._trans)=Vec3T(1,0,-1.5f).template cast<ArticulatedBody::T>();
  //simulator
  PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(Vec3T(0,0,-9.81f));
  //setup kinematic
  sim.getJointPhysicsParameter(9)._isKinematic=true;
  sim.getJointPhysicsParameter(9)._kin=[&](T t,int nrDOF) {
    ASSERT_MSG(nrDOF==3,"Incorrect DOF!")
    Vec ret=Vec::Zero(3);
    ret[0]=0;
    ret[1]=-2.5f;
    ret[2]=sin(t*2);
    return ret;
  };
  sim.getJointPhysicsParameter(10)._isKinematic=true;
  sim.getJointPhysicsParameter(10)._kin=[&](T t,int nrDOF) {
    ASSERT_MSG(nrDOF==3,"Incorrect DOF!")
    return Vec::Zero(3);
  };
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
