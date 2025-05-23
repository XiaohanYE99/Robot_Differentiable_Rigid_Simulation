#include "ArticulatedUtils.h"
#include "ArticulatedLoader.h"
#include "JointFunc.h"
#include "PBDArticulatedGradientInfo.h"
#include <Environment/MeshExact.h>
#include <Environment/BBoxExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/CompositeShapeExact.h>
#include <Environment/ConvexDecomposition.h>
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>
#include <stack>
#include <random>

namespace PHYSICSMOTION {
ArticulatedUtils::ArticulatedUtils(ArticulatedBody& body):_body(body) {}
void ArticulatedUtils::assembleGlobalVars(const tinyxml2::XMLElement& pt) {
  _globalRho=get<T>(pt,"globalRho",3);
  _globalLimitCoef=get<T>(pt,"globalLimitCoef",1000);
  _globalControlCoef=get<T>(pt,"globalControlCoef",1);
  _globalDampingCoef=get<T>(pt,"globalDampingCoef",0);
  _globalSampleDist=get<T>(pt,"globalSampleDist",std::numeric_limits<T>::max());
  _globalSampleDistLocal=get<T>(pt,"globalSampleDistLocal",1);
}
void ArticulatedUtils::assemble(const tinyxml2::XMLElement& pt) {
  assembleGlobalVars(pt);
  //build joints
  std::vector<TransInfo> infos;
  assembleJoints(*getChild(pt,"joint"),infos);
  //apply X0
  for(int i=0; i<(int) _body.nrJ(); i++) {
    Joint& J=_body.joint(i);
    J._class=i;
    J.transformMesh(infos[i]._X0);
    if(J._parent>=0) {
      const TransInfo& infoP=infos[J._parent];
      APPLY_TRANS(J._trans,infoP._X0,J._trans)
    }
  }
  //assemble joints
  int offDOF=0,offDDT=0;
  for(int i=0; i<(int) _body.nrJ(); i++) {
    const TransInfo& I=infos[i];
    Joint& J=_body.joint(i);
    J._offDOF=offDOF;
    J._offDDT=offDDT;
    offDOF+=J.nrDOF();
    offDDT+=J.nrDDT();
    J.assemble(I._rho);
  }
  _body.fillChildren();
}
void ArticulatedUtils::assembleJoints(const tinyxml2::XMLElement& pt,std::vector<TransInfo>& infos) {
  std::stack<std::pair<const tinyxml2::XMLElement*,int>> ss;
  ss.push(std::make_pair(&pt,-1));
  while(!ss.empty()) {
    std::pair<const tinyxml2::XMLElement*,int> p=ss.top();
    ss.pop();
    //assemble this joint
    assembleJoint(*(p.first),p.second,infos);
    //find children
    for(const tinyxml2::XMLElement* v=p.first->FirstChildElement(); v; v=v->NextSiblingElement()) {
      if(std::string(v->Name())=="joint")
        ss.push(std::make_pair(v,_body.nrJ()-1));
    }
  }
}
void ArticulatedUtils::assembleJoint(const tinyxml2::XMLElement& pt,int parent,std::vector<TransInfo>& infos) {
  Joint J;
  Mat3T R;
  Mat3X4T JT;
  Mat3XT DT,DDT;
  TransInfo info;
  std::string defaultName;
  TransInfo infoParent;
  if(parent>=0)
    infoParent=infos[parent];
  else {
    infoParent._TLast.setIdentity();
    infoParent._X0.setIdentity();
    infoParent._lenX.setZero();
    infoParent._lenY.setZero();
    infoParent._lenZ.setZero();
    infoParent._origin.setZero();
  }
  //basics
  if(parent<0)
    defaultName="joint";
  else defaultName=_body._joints[parent]._name+"_"+std::to_string(_body._joints.size());
  J._name=get<std::string>(pt,"name",defaultName);
  J.setType(get<std::string>(pt,"type"));
  J._parent=parent;
  J._depth=parent<0? 0: _body._joints[parent]._depth+1;
  J._offDDT=J._offDOF=0;
  J._mult=J._offset=0;
  DT.resize(3,J.nrDOF());
  DDT.resize(3,J.nrDDT());
  Mat3XTM DTMap=mapM(DT);
  Mat3XTM DDTMap=mapM(DDT);
  //geometry
  info._lenX.setZero();
  info._lenY.setZero();
  info._lenZ.setZero();
  std::vector<std::shared_ptr<ShapeExact>> geoms;
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(!beginsWith(std::string(v->Name()),"geom"))
      continue;
    const tinyxml2::XMLElement& c=*v;
    if(std::string(v->Name())=="geom3D") {
      info._rad=get<T>(c,"rad");
      info._lenX[0]=get<T>(c,"lenX");
      info._lenY[1]=get<T>(c,"lenY");
      info._lenZ[2]=get<T>(c,"lenZ");
      if(info._rad==0)
        J._mesh.reset(new BBoxExact(info._lenX.norm()/2,info._lenY.norm()/2,info._lenZ.norm()/2));
      else J._mesh.reset(new SphericalBBoxExact(info._lenX.norm()/2,info._lenY.norm()/2,info._lenZ.norm()/2,info._rad));
      J.transformMesh(Vec3T(0,info._lenY.norm()/2,0));
    } else if(std::string(v->Name())=="geom2D") {
      info._rad=get<T>(c,"rad");
      info._lenX[0]=get<T>(c,"lenX");
      info._lenY[1]=get<T>(c,"lenY");
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new SphericalBBoxExact(info._lenX.norm()/2,info._lenY.norm()/2,0,info._rad));
      J.transformMesh(Vec3T(0,info._lenY.norm()/2,0));
    } else if(std::string(v->Name())=="geom1D") {
      info._rad=get<T>(c,"rad");
      info._lenX[0]=info._rad;
      info._lenY[1]=get<T>(c,"lenY");
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new SphericalBBoxExact(0,info._lenY.norm()/2,0,info._rad));
      J.transformMesh(Vec3T(0,info._lenY.norm()/2,0));
    } else if(std::string(v->Name())=="geom0D") {
      info._rad=get<T>(c,"rad");
      info._lenX[0]=info._rad;
      info._lenY[1]=info._rad;
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new SphericalBBoxExact(info._rad));
    }
    geoms.push_back(J._mesh);
  }
  if(geoms.size()>1)
    J._mesh.reset(new CompositeShapeExact(geoms));
  //transformation: translation
  CTR(J._trans)=infoParent._origin;
  CTR(J._trans)+=infoParent._lenX*get<T>(pt,"transX",0);
  CTR(J._trans)+=infoParent._lenY*get<T>(pt,"transY",0);
  CTR(J._trans)+=infoParent._lenZ*get<T>(pt,"transZ",0);
  CTR(J._trans)+=parsePtreeDef<Vec3T>(pt,"trans",0);
  //transformation: rotation-joint
  if(hasAttribute(pt,"rotQuat")) {
    Vec4T r=parsePtree<Vec4T>(pt,"rotQuat");
    R=QuatT(r[0],r[1],r[2],r[3]).toRotationMatrix();
  } else if(J._typeJoint==Joint::ROT_3D_EXP || J._typeJoint==Joint::ROT_3D_XYZ) {
    R.col(1)=parsePtreeDef<Vec3T>(pt,"jointDirY",Vec3T::UnitY()).normalized();
    R.col(2)=parsePtreeDef<Vec3T>(pt,"jointDirZ",Vec3T::UnitZ()).normalized();
    R.col(0)=R.col(1).cross(R.col(2)).normalized();
    R.col(2)=R.col(0).cross(R.col(1)).normalized();
  } else if(J._typeJoint==Joint::BALL_JOINT) {
    if(!hasAttribute(pt,"jointDirY")) {
      Vec3T Zto=parsePtreeDef<Vec3T>(pt,"jointDirZ",Vec3T::UnitZ()).normalized();
      R=QuatT::FromTwoVectors(Vec3T::UnitZ(),Zto).toRotationMatrix();
    } else {
      R.col(1)=parsePtree<Vec3T>(pt,"jointDirY").normalized();
      R.col(2)=parsePtree<Vec3T>(pt,"jointDirZ").normalized();
      R.col(0)=R.col(1).cross(R.col(2)).normalized();
      R.col(2)=R.col(0).cross(R.col(1)).normalized();
    }
  } else if(J._typeJoint==Joint::HINGE_JOINT) {
    R=QuatT::FromTwoVectors(Vec3T::UnitZ(),parsePtreeDef<Vec3T>(pt,"jointDirZ",Vec3T::UnitZ()).normalized()).toRotationMatrix();
  } else if(J._typeJoint==Joint::TRANS_1D) {
    R=QuatT::FromTwoVectors(Vec3T::UnitX(),parsePtreeDef<Vec3T>(pt,"jointDirX",Vec3T::UnitX()).normalized()).toRotationMatrix();
  } else {
    R=ROT(infoParent._TLast);
  }
  ROT(J._trans)=ROT(infoParent._TLast).transpose()*R;
  if(hasAttribute(pt,"rotQuatAbs")) {
    Vec4T r=parsePtree<Vec4T>(pt,"rotQuatAbs");
    ROT(J._trans)=QuatT(r[0],r[1],r[2],r[3]).toRotationMatrix();
  }
  //transformation: rotation-mesh
  R.setIdentity();
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(!beginsWith(std::string(v->Name()),"geom"))
      continue;
    if(std::string(v->Name())=="geom2D" || std::string(v->Name())=="geom3D") {
      R.col(0)=parsePtree<Vec3T>(pt,"meshDirX").normalized();
      R.col(1)=parsePtree<Vec3T>(pt,"meshDirY").normalized();
      R.col(2)=R.col(0).cross(R.col(1)).normalized();
      R.col(1)=R.col(2).cross(R.col(0)).normalized();
    } else if(std::string(v->Name())=="geom1D") {
      R=QuatT::FromTwoVectors(Vec3T::UnitY(),parsePtree<Vec3T>(pt,"meshDirY").normalized()).toRotationMatrix();
    } else if(std::string(v->Name())=="geom0D" || std::string(v->Name())=="geom2Sphere" || std::string(v->Name())=="geom3Sphere") {
      R.setIdentity();
    }
  }
  R=((ROT(infoParent._TLast)*ROT(J._trans)).transpose()*R).eval();
  if(hasAttribute(pt,"rotMeshQuatAbs")) {
    Vec4T r=parsePtree<Vec4T>(pt,"rotMeshQuatAbs");
    R=QuatT(r[0],r[1],r[2],r[3]).toRotationMatrix();
  }
  J.transformMesh(R,Vec3T::Zero());
  info._lenX=R*info._lenX;
  info._lenY=R*info._lenY;
  info._lenZ=R*info._lenZ;
  //transformation: translation-mesh
  Vec3T MT=parsePtreeDef<Vec3T>(pt,"meshTrans",Vec3T::Zero());
  info._origin=info._lenX/std::max<T>(Epsilon<T>::defaultEps(),info._lenX.norm())*MT[0]+
               info._lenY/std::max<T>(Epsilon<T>::defaultEps(),info._lenY.norm())*MT[1]+
               info._lenZ/std::max<T>(Epsilon<T>::defaultEps(),info._lenZ.norm())*MT[2];
  J.transformMesh(info._origin);
  //transformation: origin
  if(hasAttribute(pt,"x0")) {
    Vec x0=parsePtree<Vec>(pt,"x0",J.nrDOF());
    VecCM x0Map(x0.data(),x0.size());
    info._X0=JointFunc<T>::TDTDDT(J,x0Map,DTMap,DDTMap);
  } else info._X0.setIdentity();
  //limit/control/damping
  J._limits.setConstant(3,J.nrDOF(),std::numeric_limits<T>::infinity());
  J._control.setConstant(J.nrDOF(),0);
  J._damping.setConstant(J.nrDOF(),0);
  if(!J.isRoot(_body) && J.nrDOF()>0) {
    const tinyxml2::XMLElement* c=getChild(pt,"limit");
    if(c) {
      J._limits.row(0)=parsePtree<Vec>(*c,"lower",J.nrDOF());
      J._limits.row(1)=parsePtree<Vec>(*c,"upper",J.nrDOF());
      J._limits.row(2)=parsePtreeDef<Vec>(*c,"coef",_globalLimitCoef,J.nrDOF());
    }
    J._control=parsePtreeDef<Vec>(pt,"controlCoef",_globalControlCoef,J.nrDOF());
    J._damping=parsePtreeDef<Vec>(pt,"dampingCoef",_globalDampingCoef,J.nrDOF());
    for(int i=0; i<J.nrDOF(); i++) {
      if(std::isfinite(J._limits(2,i))) {
        ASSERT_MSGV(J._limits(0,i)<=J._limits(1,i),"Joint: %s, limitRange error: (%f>%f)!",J._name.c_str(),J._limits(0,i),J._limits(1,i))
      }
      ASSERT_MSGV(J._control[i]>=0,"Joint: %s, controlCoef negative (%f<0)!",J._name.c_str(),J._control[i])
      ASSERT_MSGV(J._damping[i]>=0,"Joint: %s, dampingCoef negative (%f<0)!",J._name.c_str(),J._damping[i])
    }
  }
  //info
  Vec x0=Vec::Zero(J.nrDOF());
  VecCM x0Map(x0.data(),x0.size());
  JT=JointFunc<T>::TDTDDT(J,x0Map,DTMap,DDTMap);
  APPLY_TRANS(info._TLast,infoParent._TLast,J._trans);
  APPLY_TRANS(info._TLast,info._TLast,JT)
  //descend
  if(J._mesh)
    info._rho=get<T>(pt,"rho",_globalRho);
  _body._joints.push_back(J);
  infos.push_back(info);
}
//articulated body topology operations
void ArticulatedUtils::mergeMesh(Joint& joint,
                                 const std::vector<std::shared_ptr<ShapeExact>>& geomsMerged,
                                 const std::vector<CompositeShapeExact::Mat3X4T>& transMerged) {
  if(geomsMerged.empty())
    return;
  else if((int)transMerged.size()==1 && transMerged[0].isIdentity())
    joint._mesh=geomsMerged[0];
  else joint._mesh.reset(new CompositeShapeExact(geomsMerged,transMerged));
}
void ArticulatedUtils::initMesh(int k,int sz) {
  std::vector<Vec3T> vss;
  for(int i=0;i<sz;i++) {
    T theta = M_PI * (i + 0.5) / sz;
    T phi = 2 * M_PI * (i + 0.5) / sz;

    T x = sin(theta) * cos(phi);
    T y = sin(theta) * sin(phi);
    T z = cos(theta);
    vss.push_back(Vec3T(x,y,z));
  }
  ConvexHullExact m(vss);
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(m.copy());
  _body.joint(k)._mesh=mesh;
  _body.joint(k).assemble(1);
}
void ArticulatedUtils::addBase(int dim,const Vec3T& planeNormal,bool exponential) {
  //check if we already have root
  for(int i=1; i<_body.nrJ(); i++)
    if(_body.joint(i).isRoot(_body)) {
      std::cout << "Joint #" << i << " is already a root, skip addBase!" << std::endl;
      return;
    }
  if(_body.joint(0)._typeJoint!=Joint::FIX_JOINT) {
    std::cout << "Joint #0 is not fixed, skip addBase!" << std::endl;
    return;
  }
  //update index,depth
  _body._joints.insert(_body._joints.begin(),Joint());
  _body._joints.insert(_body._joints.begin(),Joint());
  _body._joints[0]._name="baseTrans";
  _body._joints[1]._name="baseRot";
  for(int i=2; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    joint._parent+=2;
    joint._depth+=2;
    if(joint._mimic>=0)
      joint._mimic+=2;
  }
  //type joint
  if(dim==2) {
    _body.joint(1)._typeJoint=Joint::HINGE_JOINT;
    _body.joint(0)._typeJoint=Joint::TRANS_2D;
  } else {
    if(exponential)
      _body.joint(1)._typeJoint=Joint::ROT_3D_EXP;
    else _body.joint(1)._typeJoint=Joint::ROT_3D_XYZ;
    _body.joint(0)._typeJoint=Joint::TRANS_3D;
  }
  for(int i=0; i<2; i++) {
    Joint& joint=_body.joint(i);
    joint._parent=i-1;
    joint._depth=i;
    //init limit,control,damping
    joint._limits.setConstant(3,joint.nrDOF(),std::numeric_limits<T>::infinity());
    joint._control.setConstant(joint.nrDOF(),0);
    joint._damping.setConstant(joint.nrDOF(),0);
    joint._trans.setIdentity();
    //assemble mass (=0)
    joint.assemble(1);
  }
  int offDOF=0,offDDT=0;
  for(int i=0; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    if(joint._typeJoint!=Joint::FIX_JOINT) {
      joint._offDOF=offDOF;
      joint._offDDT=offDDT;
      offDOF+=joint.nrDOF();
      offDDT+=joint.nrDDT();
    } else {
      joint._offDOF=offDOF;
      joint._offDDT=offDDT;
    }
  }
  _body.fillChildren();
  //2D plane
  if(dim==2) {
    QuatT quat=QuatT::FromTwoVectors(planeNormal.normalized(),Vec3T::UnitZ());
    std::set<int> children=_body.children(1,true);
    for(int i:children) {
      Joint& J=_body.joint(i);
      J._trans=(quat.toRotationMatrix()*J._trans).eval();
    }
    Joint& J0=_body.joint(0);
    J0._trans=(quat.toRotationMatrix().transpose()*J0._trans).eval();
  }
}
void ArticulatedUtils::combine(const std::vector<ArticulatedBody>& bodies) {
  _body._joints.clear();
  //joint
  {
    Joint joint;
    joint._parent=-1;
    joint._depth=0;
    joint._typeJoint=Joint::FIX_JOINT;
    joint._mimic=-1;
    joint._limits.resize(3,0);
    joint._control.resize(0);
    joint._damping.resize(0);
    joint._trans.setIdentity();
    //mimic
    joint._mult=joint._offset=0;
    //mass
    joint._M=0;
    joint._MC.setZero();
    joint._MCCT.setZero();
    //sphere approx.
    joint._name="";
    joint.assemble(1);
    _body._joints.push_back(joint);
  }
  //add bodies
  for(const ArticulatedBody& body:bodies) {
    int off=_body.nrJ();
    for(int i=0; i<body.nrJ(); i++) {
      Joint joint=body.joint(i);
      if(joint._parent==-1)
        joint._parent=0;
      else{
        joint._parent+=off;
        joint._class+=off;
      } 
      joint._depth++;
      if(joint._mimic>=0)
        joint._mimic++;
      _body._joints.push_back(joint);
    }
  }
  //reorder
  int offDOF=0,offDDT=0;
  for(int i=0; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
  }
  _body.fillChildren();
}
ArticulatedUtils::Vec ArticulatedUtils::mergeChildren(int jid,const Vec& DOF) {
  const Joint& J=_body.joint(jid);
  if(J._children.empty())
    return DOF;
  //we choose to merge other joints to the first child
  PBDArticulatedGradientInfo<T> info(_body,DOF);
  int jidc=J._children.front();
  for(int cid=1; cid<(int) J._children.size(); cid++) //merge name
    _body.joint(jidc)._name+="+"+_body.joint(J._children[cid])._name;
  Mat3X4T ITC,TC=TRANSI(info._TM,jidc);
  INV(ITC,TC)
  //compute mesh
  std::shared_ptr<CompositeShapeExact> shape;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  std::vector<ShapeExact::Mat3X4T> trans;
  for(int j=0; j<_body.nrJ(); j++) {
    if(j==jid)
      continue;
    bool eliminate=false;
    for(int jp=j; jp>=0; jp=_body.joint(jp)._parent)
      if(jp==jid)
        eliminate=true;
    //merge mesh
    const Joint& joint=_body.joint(j);
    if(eliminate && joint._mesh) {
      shapes.push_back(joint._mesh);
      Mat3X4T TRel,TJ=TRANSI(info._TM,j);
      APPLY_TRANS(TRel,ITC,TJ)
      trans.push_back(TRel.template cast<ShapeExact::T>());
    }
  }
  if(!shapes.empty()) {
    shape.reset(new CompositeShapeExact(shapes,trans));
    _body.joint(jidc)._mesh=shape;
  }
  //update body
  return eliminate([&](int j,const Joint&) {
    if(j==jidc || j==jid)
      return false;
    bool eliminate=false;
    for(int jp=j; jp>=0; jp=_body.joint(jp)._parent)
      if(jp==jid)
        eliminate=true;
    return eliminate;
  },DOF);
}
ArticulatedUtils::Vec ArticulatedUtils::fix(std::function<bool(int,const Joint&)> canFix,const Vec& DOF) {
  _canSimplify.resize(_body.nrJ());
  for(int i=0; i<_body.nrJ(); i++)
    _canSimplify[i]=canFix(i,_body.joint(i));
  Vec DOFRet=Vec::Zero(0);
  for(int j=0; j<_body.nrJ(); j++)
    if(_canSimplify[j])
      _body.joint(j)._typeJoint=Joint::FIX_JOINT;
    else DOFRet=concatRow(DOFRet,DOF.segment(_body.joint(j)._offDOF,_body.joint(j).nrDOF()));
  //update
  int offDOF=0,offDDT=0;
  for(int j=0; j<_body.nrJ(); j++) {
    Joint& joint=_body.joint(j);
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
  }
  return DOFRet;
}
ArticulatedUtils::Vec ArticulatedUtils::eliminate(std::function<bool(int,const Joint&)> canEliminate,const Vec& DOF) {
  _canSimplify.resize(_body.nrJ());
  for(int i=0; i<_body.nrJ(); i++)
    _canSimplify[i]=canEliminate(i,_body.joint(i));
  //remap joints
  int nrJLeft=0;
  Vec DOFRet=Vec::Zero(0);
  ArticulatedBody newBody;
  std::unordered_map<int,int> remap;
  for(int j=0; j<_body.nrJ(); j++) {
    bool eliminate=false;
    for(int jp=j; jp>=0; jp=_body.joint(jp)._parent)
      if(_canSimplify[jp])
        eliminate=true;
    if(eliminate)
      remap[j]=-1;
    else {
      DOFRet=concatRow(DOFRet,DOF.segment(_body.joint(j)._offDOF,_body.joint(j).nrDOF()));
      newBody._joints.push_back(_body.joint(j));
      remap[j]=nrJLeft++;
    }
  }
  //update
  int offDOF=0,offDDT=0;
  for(int j=0; j<newBody.nrJ(); j++) {
    Joint& joint=newBody.joint(j);
    if(joint._parent>=0) {
      joint._parent=remap[joint._parent];
      ASSERT_MSG(joint._parent>=0,"The body contains a joint whose parent was eliminated!")
    }
    if(joint._mimic>=0) {
      joint._mimic=remap[joint._mimic];
      ASSERT_MSG(joint._mimic>=0,"The body contains a joint whose mimic joint was eliminated!")
    }
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
    joint._depth=joint._parent<0? 0: _body._joints[joint._parent]._depth+1;
  }
  newBody.fillChildren();
  _body=newBody;
  return DOFRet;
}
ArticulatedUtils::Vec ArticulatedUtils::simplify(std::function<bool(int,const Joint&)> canSimplify,const Vec& DOF,int nrDebug) {
  ArticulatedBody bodyOriginal=_body;
  PBDArticulatedGradientInfo<T> info(_body,DOF);
  //remove fixed
  _canSimplify.resize(_body.nrJ());
  for(int i=0; i<_body.nrJ(); i++)
    _canSimplify[i]=canSimplify(i,_body.joint(i));
  for(int i=1; i<_body.nrJ(); i++) {
    const Joint& joint=_body.joint(i);
    if(_canSimplify[i]) {
      Mat3X4T trans=TRANSI(info._TK_1KM,i);
      int parentId=joint._parent;
      while(parentId>0 && _canSimplify[parentId]) {
        Mat3X4T currTrans=trans,pTrans=TRANSI(info._TK_1KM,parentId);
        APPLY_TRANS(trans,pTrans,currTrans)
        parentId=_body.joint(parentId)._parent;
      }
      Joint& parent=_body.joint(parentId);
      ASSERT(parentId==0 || !_canSimplify[parentId])
      //merge M,MC,MCCT
      parent._M+=joint._M;
      parent._MC+=ROT(trans)*joint._MC+CTR(trans)*joint._M;
      parent._MCCT+=ROT(trans)*joint._MCCT*ROT(trans).transpose();
      parent._MCCT+=ROT(trans)*joint._MC*CTR(trans).transpose();
      parent._MCCT+=CTR(trans)*joint._MC.transpose()*ROT(trans).transpose();
      parent._MCCT+=CTR(trans)*CTR(trans).transpose()*joint._M;
      parent._MCCT=((parent._MCCT+parent._MCCT.transpose())/2).eval();
      //merge name
      parent._name+="+"+joint._name;
      //merge mesh
      if(!joint._mesh)
        continue;
      else {
        std::vector<std::shared_ptr<ShapeExact>> geomsMerged;
        std::vector<CompositeShapeExact::Mat3X4T> transMerged;
        if(parent._mesh) {
          geomsMerged.push_back(parent._mesh);
          transMerged.push_back(Mat3X4T::Identity().template cast<CompositeShapeExact::T>());
        }
        geomsMerged.push_back(joint._mesh);
        transMerged.push_back(trans.template cast<CompositeShapeExact::T>());
        mergeMesh(parent,geomsMerged,transMerged);
      }
    }
  }
  //re-index
  Vec DOFRet=Vec::Zero(0);
  ArticulatedBody body;
  std::map<int,int> idMap;
  for(int i=0; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    if(_body.children(i).empty() && !joint._mesh)
      continue;
    else if(i==0 || !_canSimplify[i]) {
      while(joint._parent>0 && _canSimplify[joint._parent]) {
        Mat3X4T currTrans=joint._trans,pTrans=TRANSI(info._TK_1KM,joint._parent).template cast<T>();
        APPLY_TRANS(joint._trans,pTrans,currTrans)
        joint._parent=_body.joint(joint._parent)._parent;
      }
      if(i>0) {
        ASSERT(joint._parent==0 || !_canSimplify[joint._parent])
      }
      idMap[i]=(int)body._joints.size();
      body._joints.push_back(joint);
      DOFRet=concatRow(DOFRet,DOF.segment(joint._offDOF,joint.nrDOF()));
    } else std::cout << "Simplified Joint: " << joint._name << std::endl;
  }
  //assemble
  int offDOF=0,offDDT=0;
  for(int i=0; i<(int) body._joints.size(); i++) {
    Joint& joint=body._joints[i];
    joint._parent=i==0? -1: idMap.find(joint._parent)->second;
    if(joint._mimic>=0)
      joint._mimic=idMap.find(joint._mimic)->second;
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
  }
  body.fillChildren();
  //replace
  _body=body;
  //debug
  DEFINE_NUMERIC_DELTA_T(T)
  for(int i=0; i<nrDebug; i++) {
    Vec xOriginal=bodyOriginal.randomPose(1),x=Vec::Zero(0);
    for(int j=0; j<bodyOriginal.nrJ(); j++) {
      const Joint& J=bodyOriginal.joint(j);
      if(j==0 || !_canSimplify[j])
        x=concat<Vec>(x,xOriginal.segment(J._offDOF,J.nrDOF()));
      else xOriginal.segment(J._offDOF,J.nrDOF())=DOF.segment(J._offDOF,J.nrDOF());
    }
    //EOriginal
    T EOriginal=0;
    PBDArticulatedGradientInfo<T> infoOriginal(bodyOriginal,xOriginal);
    for(int j=0; j<bodyOriginal.nrJ(); j++) {
      const Joint& J=bodyOriginal.joint(j);
      const Mat3T& PPT=J._MCCT;
      const Vec3T& P=J._MC;
      const Mat3X4T A=TRANSI(infoOriginal._TM,j);
      EOriginal+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
    }
    //E
    T E=0;
    ASSERT_MSG(x.size()==body.nrDOF(),"Error estimating simplified body DOF!")
    info=PBDArticulatedGradientInfo<T>(body,x);
    for(int j=0; j<body.nrJ(); j++) {
      const Joint& J=body.joint(j);
      const Mat3T& PPT=J._MCCT;
      const Vec3T& P=J._MC;
      const Mat3X4T A=TRANSI(info._TM,j);
      E+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
    }
    DEBUG_GRADIENT("simplifyConsistency",E,E-EOriginal)
  }
  return DOFRet;
}
void ArticulatedUtils::resetMesh(const std::string& path,int k,T sz) {
  MeshExact m(path);
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(m.copy());
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  mesh->getMesh(vss,iss);
  Vec3T center=Vec3T(0,0,0);
  //std::cout<<mesh->vss().size()<<" "<<vss.size();
  //for(int i=0;i<mesh->vss().size();i++) center+=mesh->vss()[i]/(mesh->vss().size());
  for(int i=0;i<mesh->vss().size();i++) vss[i]=vss[i]*sz;
  _body.joint(k)._mesh.reset(new MeshExact(vss,iss));
  _body.joint(k).assemble(1);//50
}
ArticulatedUtils::Vec ArticulatedUtils::simplify(const Vec& DOF,int nrDebug) {
  return simplify([](int,const Joint& J) {
    return J._typeJoint==Joint::FIX_JOINT;
  },DOF,nrDebug);
}
ArticulatedUtils::Vec ArticulatedUtils::replaceJoint(const Vec& DOF,int jid,Mat3XT axes) {
  PBDArticulatedGradientInfo<T> info(_body,DOF);
  Joint& J=_body._joints[jid];
  Mat3X4T trans;
  Mat3T R;
  if(J._parent>=0) {
    Mat3X4T TP=TRANSI(info._TM,J._parent);
    APPLY_TRANS(trans,TP,J._trans);
  } else trans=J._trans;
  axes=ROT(trans).transpose()*axes;
  if(axes.cols()==0) {
    J._typeJoint=Joint::FIX_JOINT;
    R.setIdentity();
  } else if(axes.cols()==1) {
    J._typeJoint=Joint::HINGE_JOINT;
    R=QuatT::FromTwoVectors(axes.col(0),Vec3T::UnitZ()).toRotationMatrix();
  } else if(axes.cols()==2) {
    J._typeJoint=Joint::BALL_JOINT;
    axes=concatCol<Mat3XT>(axes.col(0).cross(axes.col(1)),axes);
    R=axes;
  } else if(axes.cols()==3) {
    J._typeJoint=Joint::ROT_3D_XYZ;
    R=axes;
  } else {
    ASSERT_MSGV(false,"Unsupported joint type with %d axes!",(int) axes.cols())
  }
  //transform joint & children
  ROT(J._trans)=ROT(J._trans)*R;
  J.transformMesh(R.transpose(),Vec3T::Zero());
  for(int jidc:J._children) {
    Joint& JC=_body._joints[jidc];
    JC._trans=R.transpose()*JC._trans;
  }
  //control/damping/limits
  if(J._control.size()>0)
    J._control.setConstant(J.nrDOF(),J._control[0]);
  else J._control.setZero(J.nrDOF());
  if(J._damping.size()>0)
    J._damping.setConstant(J.nrDOF(),J._damping[0]);
  else J._damping.setZero(J.nrDOF());
  if(J._limits.size()>0)
    J._limits=J._limits.col(0)*Vec::Ones(J.nrDOF()).transpose();
  else J._limits.setConstant(3,J.nrDOF(),std::numeric_limits<T>::infinity());
  //reorder
  int offDOF=0,offDDT=0;
  Vec DOFRet=Vec::Zero(0);
  for(int i=0; i<(int) _body._joints.size(); i++) {
    Joint& joint=_body._joints[i];
    if(joint._mimic>=0) {
      ASSERT_MSGV(joint._mimic!=jid,"Cannot replace joint %d because joint %d mimics it!",jid,i)
    }
    if(i==jid)
      DOFRet=concatRow(DOFRet,Vec::Zero(joint.nrDOF()));
    else DOFRet=concatRow(DOFRet,DOF.segment(joint._offDOF,joint.nrDOF()));
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
  }
  return DOFRet;
}
void ArticulatedUtils::convexDecompose(T rho) {
  int meshJointId=-1;
  for(int i=0; i<_body.nrJ(); i++) {
    std::cout << _body.joint(i)._trans << std::endl;
    if(_body.joint(i)._typeJoint!=Joint::FIX_JOINT) {
      ASSERT_MSG(false,"Cannot decompose ArticulatedBody with non-fixed joints!")
    } else if(_body.joint(i)._mesh) {
      ASSERT_MSG(meshJointId==-1,"Can only decompose ArticulatedBody with a single mesh!")
      meshJointId=i;
    }
  }
  ASSERT_MSG(meshJointId>=0,"No mesh to decompose in ArticulatedBody!")
  //transMesh
  PBDArticulatedGradientInfo<ArticulatedBody::T> info(_body,Vec::Zero(_body.nrDOF()));
  Mat3X4T transMesh=TRANSI(info._TM,meshJointId);
  //decompose
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  _body.joint(meshJointId)._mesh->getMesh(vss,iss);
  MeshExact mesh;
  //mesh.init(vss,iss);
  ConvexDecomposition decompose(mesh);
  _body=ArticulatedLoader::createDummy();
  for(int i=0; i<(int)decompose.getConvexHulls().size(); i++) {
    Joint joint;
    //joint
    joint._parent=0;
    joint._depth=_body.joint(joint._parent)._depth+1;
    joint._typeJoint=Joint::FIX_JOINT;
    joint._offDOF=0;
    joint._offDDT=0;
    joint._limits.setZero(3,0);
    joint._control.setZero(0);
    joint._damping.setZero(0);
    joint._trans.setIdentity();
    //mimic
    joint._mult=joint._offset=0;
    //mesh approx.
    joint._mesh=decompose.getConvexHulls()[i];
    joint.transformMesh(transMesh);
    //mass
    joint._name="Convex"+std::to_string(i);
    joint.assemble(rho);
    //VTKWriter<double> os(joint._name,joint._name+".vtk",true);
    //decompose.getConvexHulls()[i]->writeVTK(os,Mat3X4T::Identity());
    //std::cout << decompose.getConvexHulls()[i]->vss().size() << std::endl;
    //std::cout << transMesh << std::endl;
    _body._joints.push_back(joint);
  }
  std::cout << "Decomposed into " << decompose.getConvexHulls().size() << " sub-joints!" << std::endl;
}
void ArticulatedUtils::convexDecompose(std::vector<int> jid, int maxConvexHulls, T rho) {
  /*std::vector<MeshExact> Mesh;
  PBDArticulatedGradientInfo<ArticulatedBody::T> info(_body,Vec::Zero(_body.nrDOF()));
  for(uint i=0;i<jid.size();i++){
    int meshJointId=jid[i];//transMesh
    //decompose
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    _body.joint(meshJointId)._mesh->getMesh(vss,iss);
    MeshExact mesh;
    mesh.init(vss,iss);
    Mesh.push_back(mesh);
  }
  for(uint k=0;k<Mesh.size();k++){
    int meshJointId=jid[k];//transMesh
    ConvexDecomposition decompose(Mesh[k],1.,Vec3T(0,0,0),Mat3T::Identity(),maxConvexHulls);
    //Mat3X4T transMesh=TRANSI(info._TM,meshJointId);
    _body.joint(meshJointId)._mesh=decompose.getConvexHulls()[0];
    _body.joint(meshJointId)._name="Convex"+std::to_string(0);
    _body.joint(meshJointId).assemble(rho);
    
    for(int i=1; i<(int)decompose.getConvexHulls().size(); i++) {
      Joint joint;
      //joint
      joint._parent=meshJointId;
      joint._depth=_body.joint(joint._parent)._depth+1;
      joint._typeJoint=Joint::FIX_JOINT;
      joint._offDOF=_body.joint(meshJointId)._offDOF;
      joint._offDDT=_body.joint(meshJointId)._offDDT;
      joint._limits=_body.joint(meshJointId)._limits;
      joint._control=_body.joint(meshJointId)._control;
      joint._damping=_body.joint(meshJointId)._damping;
      joint._trans.setIdentity();//_body.joint(meshJointId)._trans;
      //mimic
      joint._mult=_body.joint(meshJointId)._mult;
      joint._offset=_body.joint(meshJointId)._offset;
      //mesh approx.
      joint._mesh=decompose.getConvexHulls()[i];
      //joint.transformMesh(transMesh);
      //mass
      joint._name="Convex"+std::to_string(i);
      joint.assemble(rho);
      joint._class=meshJointId;
      _body._joints.push_back(joint);
      //_body._joints[_body._joints.size()-1]._class=meshJointId;
    }
    std::cout << "Decomposed into " << decompose.getConvexHulls().size() << " sub-joints!" << std::endl;
  }*/
}
void ArticulatedUtils::addBody(ArticulatedBody& body) {
  int offset=(int) _body._joints.size();
  std::vector<Joint> joints=body._joints;
  if(!_body._joints.empty()) {
    for(auto& joint:joints) {
      if(joint._parent!=-1)
        joint._parent+=offset;
      for(int& child:joint._children)
        child+=offset;
      joint._offDOF+=_body.nrDOF();
      joint._offDDT+=_body.nrDDT();
    }
    joints[0]._parent=0;
    _body._joints[0]._children.push_back(offset);
    _body._joints.insert(_body._joints.end(),joints.begin(),joints.end());
  } else {
    _body._joints.insert(_body._joints.end(),joints.begin(),joints.end());
  }
}
void ArticulatedUtils::untessellate(bool rebuildBVH, int refine) {
  for (int j = 0; j < _body.nrJ(); j++) {
    Joint& joint = _body.joint(j);
    if (joint._mesh) {
      // 获取当前网格的顶点和索引
      std::vector<Eigen::Matrix<double, 3, 1>> vss; // 顶点列表
      std::vector<Eigen::Matrix<int, 3, 1>> iss;   // 面索引列表
      joint._mesh->getMesh(vss, iss);

      if (refine > 0 && j >= 2) {
        int n = refine; // 需要移除的顶点数量
        if (n >= vss.size() - 3) { // 确保至少有 3 个顶点剩余
          std::cerr << "Error: Cannot remove more vertices than allowed." << std::endl;
          return;
        }

        // **逐步减面**
        for (int removed = 0; removed < n; ++removed) {
          // **选择一个适合移除的顶点**
          int vertexToRemove = -1;
          double minImpact = std::numeric_limits<double>::max();

          // 遍历所有顶点，选择对形状和分布影响最小的顶点
          for (size_t i = 0; i < vss.size(); ++i) {
            // 检查顶点是否可以移除（必须与至少 3 个三角形相连）
            std::set<int> connectedVertices;
            for (const auto& face : iss) {
              if (face(0) == i || face(1) == i || face(2) == i) {
                connectedVertices.insert(face(0));
                connectedVertices.insert(face(1));
                connectedVertices.insert(face(2));
              }
            }
            if (connectedVertices.size() < 3) continue; // 跳过不可移除的顶点

            // 计算移除顶点的几何代价（法向量变化）
            double impact = 0.0;
            for (const auto& face : iss) {
              if (face(0) == i || face(1) == i || face(2) == i) {
                const Eigen::Matrix<double, 3, 1>& v0 = vss[face(0)];
                const Eigen::Matrix<double, 3, 1>& v1 = vss[face(1)];
                const Eigen::Matrix<double, 3, 1>& v2 = vss[face(2)];

                Eigen::Matrix<double, 3, 1> normal = (v1 - v0).cross(v2 - v0).normalized();
                Eigen::Matrix<double, 3, 1> centroid = (v0 + v1 + v2) / 3.0;
                Eigen::Matrix<double, 3, 1> newNormal = (v1 - centroid).cross(v2 - centroid).normalized();
                impact += (normal - newNormal).norm();
              }
            }

            // 如果代价更小，则更新最优选择
            if (impact < minImpact) {
              minImpact = impact;
              vertexToRemove = i;
            }
          }

          if (vertexToRemove == -1) {
            std::cerr << "No removable vertex found." << std::endl;
            break;
          }

          // **移除选定的顶点并重新三角化**
          std::vector<int> connectedVertices;
          for (auto it = iss.begin(); it != iss.end();) {
            if ((*it)(0) == vertexToRemove || (*it)(1) == vertexToRemove || (*it)(2) == vertexToRemove) {
              for (int k = 0; k < 3; ++k) {
                if ((*it)(k) != vertexToRemove) {
                  connectedVertices.push_back((*it)(k));
                }
              }
              it = iss.erase(it); // 删除包含该顶点的三角形
            } else {
              ++it;
            }
          }

          // 去重并排序相邻顶点（按逆时针顺序）
          Eigen::Matrix<double, 3, 1> center = vss[vertexToRemove];
          std::sort(connectedVertices.begin(), connectedVertices.end());
          connectedVertices.erase(std::unique(connectedVertices.begin(), connectedVertices.end()), connectedVertices.end());
          std::sort(connectedVertices.begin(), connectedVertices.end(),
                    [&vss, &center](int a, int b) {
                      Eigen::Matrix<double, 3, 1> va = vss[a] - center;
                      Eigen::Matrix<double, 3, 1> vb = vss[b] - center;
                      return std::atan2(va.y(), va.x()) < std::atan2(vb.y(), vb.x());
                    });

          // 重新生成三角形
          for (size_t i = 0; i < connectedVertices.size(); ++i) {
            int v0 = connectedVertices[i];
            int v1 = connectedVertices[(i + 1) % connectedVertices.size()];
            iss.push_back((Eigen::Matrix<int, 3, 1>() << v0, v1, vertexToRemove).finished());
          }

          // 从顶点列表中移除该顶点
          vss.erase(vss.begin() + vertexToRemove);

          // 更新索引列表（调整后续顶点索引）
          for (auto& face : iss) {
            for (int k = 0; k < 3; ++k) {
              if (face(k) > vertexToRemove) {
                --face(k);
              }
            }
          }
        }
        std::cout<<iss.size();
        // 更新 mesh
        joint._mesh.reset(new MeshExact(vss, iss));
      } else {
        // 如果 refine <= 0 或 j < 2，则保留原始网格
        joint._mesh.reset(new MeshExact(vss, iss));
      }
    }
  }
}
//mesh operation
void ArticulatedUtils::tessellate(bool rebuildBVH,int refine) {
  for(int j=0; j<_body.nrJ(); j++) {
    Joint& joint=_body.joint(j);
    if(joint._mesh) {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      std::vector<Eigen::Matrix<double,3,1>> new_vss;
      std::vector<Eigen::Matrix<int,3,1>> new_iss;
      joint._mesh->getMesh(vss,iss);
      if(refine>0 && j>=2){
        new_vss = vss;
    new_iss = iss;

    // **计算整个网格的几何中心**
    Eigen::Matrix<double, 3, 1> mesh_center = Eigen::Matrix<double, 3, 1>::Zero();
    for (const auto& vertex : vss) {
        mesh_center += vertex;
    }
    mesh_center /= vss.size();  // 网格几何中心

    // 随机数生成器，用于在距离权重中加入少量随机性
    std::random_device rd;
    std::mt19937 gen(rd());

    int addedPoints = 0;

    // 新增点直到达到 n 个
    while (addedPoints < refine) {
        // **计算每个三角形质心到网格中心的距离**
        std::vector<std::pair<int, double>> faceDistances;  // 保存三角形索引及其权重
        for (size_t i = 0; i < new_iss.size(); ++i) {
            const auto& face = new_iss[i];
            // 获取当前三角形的顶点坐标
            const Eigen::Matrix<double, 3, 1>& v0 = new_vss[face(0)];
            const Eigen::Matrix<double, 3, 1>& v1 = new_vss[face(1)];
            const Eigen::Matrix<double, 3, 1>& v2 = new_vss[face(2)];

            // 计算三角形的质心
            Eigen::Matrix<double, 3, 1> centroid = (v0 + v1 + v2) / 3.0;

            // 计算质心到网格中心的距离
            double distance = (centroid - mesh_center).norm();
            faceDistances.emplace_back(i, distance);
        }

        // **选择距离网格中心较近的一个三角形**
        // 按距离排序（从近到远）
        std::sort(faceDistances.begin(), faceDistances.end(),
                  [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                      return a.second < b.second;
                  });

        // 引入随机性，在前 10% 中随机选择一个三角形
        int maxCandidates = std::max(1, static_cast<int>(faceDistances.size() * 0.3));
        std::uniform_int_distribution<> dist(0, maxCandidates - 1);
        int selectedFaceIndex = faceDistances[dist(gen)].first;

        // 获取选中的三角形
        const auto& face = new_iss[selectedFaceIndex];
        int v0_index = face(0);
        int v1_index = face(1);
        int v2_index = face(2);

        // 获取当前三角形的顶点坐标
        const Eigen::Matrix<double, 3, 1>& v0 = new_vss[v0_index];
        const Eigen::Matrix<double, 3, 1>& v1 = new_vss[v1_index];
        const Eigen::Matrix<double, 3, 1>& v2 = new_vss[v2_index];

        // **计算三角形的中点（质心）**
        Eigen::Matrix<double, 3, 1> centroid = (v0 + v1 + v2) / 3.0;

        // 添加新点到顶点列表，并记录新顶点索引
        int centroid_index = new_vss.size();  // 新点索引为当前顶点列表的末尾
        new_vss.push_back(centroid);

        // **细分该三角形为三个新的三角形**
        new_iss.push_back((Eigen::Matrix<int, 3, 1>() << v0_index, v1_index, centroid_index).finished());
        new_iss.push_back((Eigen::Matrix<int, 3, 1>() << v1_index, v2_index, centroid_index).finished());
        new_iss.push_back((Eigen::Matrix<int, 3, 1>() << v2_index, v0_index, centroid_index).finished());

        // 删除原三角形
        new_iss[selectedFaceIndex] = new_iss.back();
        new_iss.pop_back();  // 移除最后一个面（原三角形）

        // 更新新增点计数
        ++addedPoints;
    }
        joint._mesh.reset(new MeshExact(new_vss,new_iss,rebuildBVH));
        //joint.assemble(1);
      }
      else joint._mesh.reset(new MeshExact(vss,iss,rebuildBVH));
    }
  }
}
void ArticulatedUtils::BBApproxiate(bool rebuildBVH) {
  for(int j=0; j<_body.nrJ(); j++) {
    Joint& joint=_body.joint(j);
    if(joint._mesh) {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      addBox(vss,iss,
             joint._mesh->getBB().minCorner().template cast<double>(),
             joint._mesh->getBB().maxCorner().template cast<double>());
      joint._mesh.reset(new MeshExact(vss,iss,rebuildBVH));
    }
  }
}
void ArticulatedUtils::makeConvex() {
  for(int j=0; j<_body.nrJ(); j++) {
    Joint& joint=_body.joint(j);
    if(joint._mesh) {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      joint._mesh->getMesh(vss,iss);
      joint._mesh.reset(new ConvexHullExact(vss));
    }
  }
}
//rigid transformation
ArticulatedUtils::Mat3X4T ArticulatedUtils::transformTorso(const Mat3X4T& trans) {
  Joint& J=_body.joint(0);
  J.transformMesh(trans);
  for(int c:J._children) {
    Joint& cJ=_body.joint(c);
    Mat3X4T transC=cJ._trans;
    APPLY_TRANS(cJ._trans,trans,transC)
  }
  return trans;
}
ArticulatedUtils::Mat3X4T ArticulatedUtils::transformTorso(const Mat3T& R) {
  Mat3X4T trans=Mat3X4T::Identity();
  ROT(trans)=R;
  return transformTorso(trans);
}
ArticulatedUtils::Mat3X4T ArticulatedUtils::transformTorso(const Vec3T& t) {
  Mat3X4T trans=Mat3X4T::Identity();
  CTR(trans)=t;
  return transformTorso(trans);
}
ArticulatedUtils::Mat3X4T ArticulatedUtils::transformTorsoToCOM() {
  Joint& J=_body.joint(_body.rootJointId());
  if(J._mesh) {
    Mat3X4T trans;
    ROT(trans).setIdentity();
    CTR(trans)=-J.getC();
    transformTorso(trans);
    return trans;
  } else return Mat3X4T::Identity();
}
ArticulatedUtils::Mat3X4T ArticulatedUtils::scaleBody(T coef) {
  for(Joint& J:_body._joints) {
    CTR(J._trans)*=coef;
    J._mesh->scale(coef);
    J._M*=std::pow(coef,3);
    J._MC*=std::pow(coef,4);
    J._MCCT*=std::pow(coef,5);
    if(J._mesh)
      J._mesh->scale(coef);
    if(!J.isRotational()) {
      for(int d=0; d<J._limits.cols(); d++)
        if(std::isfinite(J._limits(2,d))) {
          J._limits(0,d)*=coef;
          J._limits(1,d)*=coef;
        }
      if(J._mimic>=0)
        J._offset*=coef;
    }
  }
  return Mat3X4T::Identity()*coef;
}
//mass transformation
void ArticulatedUtils::scaleMass(T coef) {
  for(int j=0; j<_body.nrJ(); j++) {
    Joint& J=_body.joint(j);
    J._M*=coef;
    J._MC*=coef;
    J._MCCT*=coef;
  }
}
ArticulatedUtils::T ArticulatedUtils::totalMass() const {
  T M=0;
  for(int j=0; j<_body.nrJ(); j++) {
    const Joint& J=_body.joint(j);
    M+=J._M;
  }
  return M;
}
}
