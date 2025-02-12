#include "Environment/SphericalBBoxExact.h"
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/XPBDSimulator.h>
#include <Utils/Utils.h>

using namespace PHYSICSMOTION;

void readJoints(std::vector<Joint>& joints, int parentId, const tinyxml2::XMLElement* g) {
  Joint joint;
  joint._parent=parentId;
  joint._name=get<std::string>(*g,"<xmlattr>.name");
  //compute depth
  joint._depth=parentId==-1?0:joints[parentId]._depth+1;
  //read trans
  {
    ArticulatedBody::Vec3T v=parsePtreeDef<ArticulatedBody::Vec3T>(*g,"<xmlattr>.pos","0 0 0");
    ArticulatedBody::Vec4T q=parsePtreeDef<ArticulatedBody::Vec4T>(*g,"<xmlattr>.quat","1 0 0 0");
    CTR(joint._trans)=v;
    ROT(joint._trans)=Eigen::Quaternion<ArticulatedBody::T>(q[0],q[1],q[2],q[3]).toRotationMatrix();
  }
  //read geometry
  if(g->FirstChildElement("geom")->FindAttribute("type") == NULL) {
    //capsule
    ArticulatedBody::Mat3X4T trans;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    ArticulatedBody::Vec6T ft=parsePtreeDef<ArticulatedBody::Vec6T>(*gg,"<xmlattr>.fromto","0 0 0 0 0 0");
    //ArticulatedBody::Vec4T q=parsePtreeDef<ArticulatedBody::Vec4T>(*gg,"<xmlattr>.quat","1 0 0 0");
    ArticulatedBody::T size=get<ArticulatedBody::T>(*gg,"<xmlattr>.size");
    ArticulatedBody::T len=(ft.segment<3>(3)-ft.segment<3>(0)).norm();
    CTR(trans)=(ft.segment<3>(3)+ft.segment<3>(0))/2;
    ROT(trans)=Eigen::Quaternion<ArticulatedBody::T>::FromTwoVectors(ArticulatedBody::Vec3T::UnitX(),(ft.segment<3>(3)-ft.segment<3>(0)).normalized()).toRotationMatrix();
    //trans=(Eigen::Quaternion<ArticulatedBody::T>(q[0],q[1],q[2],q[3]).toRotationMatrix()*trans).eval();
    std::shared_ptr<SphericalBBoxExact> shape(new SphericalBBoxExact(len/2,size));
    joint._mesh=std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({shape}, {trans.template cast<GEOMETRY_SCALAR>()}));
    joint.assemble(get<ArticulatedBody::T>(*gg,"<xmlattr>.density"));
  } else {
    //box
    ArticulatedBody::Mat3X4T trans;
    const tinyxml2::XMLElement* gg=g->FirstChildElement("geom");
    ArticulatedBody::Vec3T p=parsePtreeDef<ArticulatedBody::Vec3T>(*gg,"<xmlattr>.pos","0 0 0");
    ArticulatedBody::Vec4T q=parsePtreeDef<ArticulatedBody::Vec4T>(*gg,"<xmlattr>.quat","1 0 0 0");
    ArticulatedBody::Vec3T size=parsePtreeDef<ArticulatedBody::Vec3T>(*gg,"<xmlattr>.size","0 0 0");
    CTR(trans)=p;
    ROT(trans)=Eigen::Quaternion<ArticulatedBody::T>(q[0],q[1],q[2],q[3]).toRotationMatrix();
    std::shared_ptr<BBoxExact> shape(new BBoxExact(size[0],size[1],size[2]));
    joint._mesh=std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({shape}, {trans.template cast<GEOMETRY_SCALAR>()}));
    joint.assemble(get<ArticulatedBody::T>(*gg,"<xmlattr>.density"));
  }
  //read joints
  {
    std::vector<const tinyxml2::XMLElement*> joints;
    for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
      if(std::string(gc->Name()) == "joint")
        joints.push_back(gc);
    joint._typeJoint=joints.empty()?Joint::FIX_JOINT:Joint::ROT_3D_XYZ;
    joint._limits.setZero(3,3);
    for(int i=0; i<(int)joints.size(); i++) {
      int id=-1;
      ArticulatedBody::Vec3T a=parsePtreeDef<ArticulatedBody::Vec3T>(*joints[0],"<xmlattr>.axis","0 0 0");
      ArticulatedBody::Vec2T range=parsePtreeDef<ArticulatedBody::Vec2T>(*joints[i],"<xmlattr>.range","0 0");
      a.cwiseAbs().maxCoeff(&id);
      joint._limits(0,id)=range[0]*M_PI/180;
      joint._limits(1,id)=range[1]*M_PI/180;
    }
    joint._limits.row(2).setConstant(1000);
    joint._control.resize(0);
    joint._damping.resize(0);
  }
  //read children
  parentId=(int)joints.size();
  joints.push_back(joint);
  for(const tinyxml2::XMLElement* gc=g->FirstChildElement(); gc; gc=gc->NextSiblingElement())
    if(std::string(gc->Name()) == "body")
      readJoints(joints,parentId,gc);
}
ArticulatedBody readMJCF(const std::string& file) {
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  tinyxml2::XMLElement* link=pt.RootElement();
  //joints
  std::vector<Joint> joints;
  for(const tinyxml2::XMLElement* g=link->FirstChildElement(); g; g=g->NextSiblingElement())
    if(std::string(g->Name()) == "worldbody")
      readJoints(joints,-1,g->FirstChildElement("body"));
  //body
  int nrDOF=0;
  int nrDDT=0;
  ArticulatedBody body;
  body=body.resizeJoints(joints.size());
  for(int i=0; i<(int)joints.size(); i++) {
    joints[i]._offDOF=nrDOF;
    joints[i]._offDDT=nrDDT;
    nrDOF+=joints[i].nrDOF();
    nrDDT+=joints[i].nrDDT();
    body.joint(i)=joints[i];
  }
  body.fillChildren();
  ArticulatedUtils utils(body);
  utils.addBase(3,ArticulatedBody::Vec3T::UnitZ());
  return body;
}
int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(readMJCF("human.xml")));
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-5,-5,-5),BBoxExact::Vec3T(5,5,-2)));
  //simulator
  XPBDSimulator sim(0.033f);
  sim.setArticulatedBody(body);
  sim.setGravity(Vec3T(0,0,-9.81f));
  sim.setHeuristcGuessStiffness();
  sim.addShape(floor);
  sim.setOutput(true);
  //run app
  visualizeSimulator(argc,argv,sim);
  return 0;
}
