#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/PBDSimulator.h>
#include <Utils/RotationUtils.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3X4T trans;
  trans.setIdentity();
  ROT(trans)=expWGradV<T,Vec3T>(Vec3T::UnitX()*M_PI/4,NULL,NULL);
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createBox(*(pt.RootElement()),2.5f,2.5f,2.5f,0);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-2.5f,-2.5f,-2.5f),BBoxExact::Vec3T(2.5f,2.5f,2.5f)));
  std::shared_ptr<CompositeShapeExact> shape(new CompositeShapeExact({floor}, {trans.template cast<GEOMETRY_SCALAR>()}));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  utils.addBase(3,Vec3T::Zero().template cast<ArticulatedBody::T>());
  CTR(body->joint(1)._trans)=Vec3T(0,-1.25f,5).template cast<ArticulatedBody::T>();
  //simulator
  PBDSimulator sim(0.01f);
  sim.setArticulatedBody(body);
  sim.setGravity(Vec3T(0,0,-9.81f));
  sim.setHeuristcGuessStiffness();
  sim.addShape(shape);
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
