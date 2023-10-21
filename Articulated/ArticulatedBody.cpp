#include "ArticulatedBody.h"
#include "ArticulatedUtils.h"
#include "PBDArticulatedGradientInfo.h"
#include <Utils/RotationUtils.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include <random>

namespace PHYSICSMOTION {
//ArticulatedBody
ArticulatedBody::ArticulatedBody() {}
ArticulatedBody::ArticulatedBody(const tinyxml2::XMLElement& pt) {
  ArticulatedUtils utils(*this);
  utils.assemble(pt);
}
bool ArticulatedBody::read(std::istream& is,IOData* dat) {
  readBinaryData(_joints,is,dat);
  return is.good();
}
bool ArticulatedBody::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_joints,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ArticulatedBody::copy() const {
  return std::shared_ptr<SerializableBase>(new ArticulatedBody);
}
std::string ArticulatedBody::type() const {
  return typeid(ArticulatedBody).name();
}
void ArticulatedBody::randomize(int nrLink,bool chain) {
  int offDOF=0;
  int offDDT=0;
  _joints.resize(nrLink);
  std::random_device rd;
  std::mt19937 gen(rd());
  for(int i=0; i<nrLink; i++) {
    if(chain)
      _joints[i]._parent=i-1;
    else
      _joints[i]._parent=i==0?-1:std::uniform_int_distribution<>(0,i-1)(gen);
    _joints[i]._M=std::uniform_real_distribution<>(0.5,1)(gen);
    _joints[i]._MC=_joints[i]._M*Vec3T::Random();
    _joints[i]._MCCT=_joints[i]._MC*_joints[i]._MC.transpose()/_joints[i]._M;
    _joints[i]._typeJoint=1<<std::uniform_int_distribution<>(0,7)(gen);
    CTR(_joints[i]._trans).setRandom();
    ROT(_joints[i]._trans)=expWGradV<T,Vec3T>(Vec3T::Random()*M_PI);
    _joints[i]._limits.resize(3,_joints[i].nrDOF());
    for(int d=0; d<_joints[i].nrDOF(); d++) {
      _joints[i]._limits(0,d)=-M_PI;
      _joints[i]._limits(1,d)=M_PI;
      _joints[i]._limits(2,d)=1;
    }
    _joints[i]._offDOF=offDOF;
    _joints[i]._offDDT=offDDT;
    offDOF+=_joints[i].nrDOF();
    offDDT+=_joints[i].nrDDT();
  }
  fillChildren();
  //I should assemble _geom in this function,
  //but I decide not to do this to save computation, as randomize is for debug only
}
std::set<int> ArticulatedBody::children(int id,bool direct) const {
  std::set<int> ret;
  for(int i=0; i<(int)nrJ(); i++)
    if(direct) {
      if(joint(i)._parent==id)
        ret.insert(i);
    } else {
      for(int j=i; j>=0; j=joint(j)._parent)
        if(j==id)
          ret.insert(i);
    }
  return ret;
}
int ArticulatedBody::commonRoot(int id,int id2) const {
  if(id==id2)
    return id;
  //idChain
  std::vector<int> idChain;
  while(id>=0) {
    idChain.push_back(id);
    id=joint(id)._parent;
  }
  //id2Chain
  std::vector<int> id2Chain;
  while(id2>=0) {
    id2Chain.push_back(id2);
    id2=joint(id2)._parent;
  }
  //compute common root
  int ret=-1;
  while(idChain.back()==id2Chain.back()) {
    ret=idChain.back();
    idChain.pop_back();
    id2Chain.pop_back();
  }
  return ret;
}
int ArticulatedBody::hasJoint(Joint::JOINT_TYPE type) const {
  for(int i=0; i<nrJ(); i++)
    if(joint(i)._typeJoint==type)
      return i;
  return -1;
}
ArticulatedBody ArticulatedBody::resizeJoints(int nr) const {
  ArticulatedBody body;
  body._joints=_joints;
  body._joints.resize(nr);
  body.fillChildren();
  return body;
  //I should assemble _geom in this function,
  //but I decide not to do this to save computation, as resizeJoints is used by NESimplifiedDynamics only
}
bool ArticulatedBody::isLeaf(int id) const {
  return children(id).size()==1;
}
ArticulatedBody::Vec ArticulatedBody::control() const {
  Vec ret=Vec::Zero(0);
  for(int j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._control);
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::damping() const {
  Vec ret=Vec::Zero(0);
  for(int j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._damping);
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::coefLimit() const {
  Vec ret=Vec::Zero(0);
  for(int j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(2).transpose());
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::lowerLimit() const {
  Vec ret=Vec::Zero(0);
  for(int j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(0).transpose());
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::upperLimit() const {
  Vec ret=Vec::Zero(0);
  for(int j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(1).transpose());
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::lowerLimit(T infty) const {
  Vec ret=lowerLimit();
  for(int i=0; i<ret.size(); i++)
    if(!std::isfinite(ret[i]))
      ret[i]=-infty;
  return ret;
}
ArticulatedBody::Vec ArticulatedBody::upperLimit(T infty) const {
  Vec ret=upperLimit();
  for(int i=0; i<ret.size(); i++)
    if(!std::isfinite(ret[i]))
      ret[i]=infty;
  return ret;
}
const ArticulatedBody::Vec& ArticulatedBody::clampLimit(Vec& x) const {
  for(int j=0; j<nrJ(); j++) {
    const Joint& J=joint(j);
    for(int i=0; i<J.nrDOF(); i++)
      if(std::isfinite(J._limits(2,i))) {
        T& xEntry=x[J._offDOF+i];
        xEntry=std::min(std::max(xEntry,J._limits(0,i)),J._limits(1,i));
      }
  }
  return x;
}
ArticulatedBody::Vec ArticulatedBody::randomPose(T coef) const {
  Vec x=Vec::Zero(nrDOF());
  std::random_device rd;
  std::mt19937 gen(rd());
  for(int j=0; j<nrJ(); j++) {
    const Joint& J=joint(j);
    for(int i=0; i<J.nrDOF(); i++) {
      T& xEntry=x[J._offDOF+i];
      if(std::isfinite(J._limits(2,i)))
        xEntry=std::uniform_real_distribution<>(J._limits(0,i),J._limits(1,i))(gen);
      else xEntry=std::uniform_real_distribution<>(-coef,coef)(gen);
    }
  }
  return x;
}
ArticulatedBody::Mat3XT ArticulatedBody::getT(const Vec& x) const {
  PBDArticulatedGradientInfo<T> info(*this,x);
  return info._TM;
}
BBoxExact ArticulatedBody::getBB(const Vec& x) const {
  return getBB(getT(x));
}
BBoxExact ArticulatedBody::getBB(const Mat3XT& t) const {
  BBoxExact bb;
  for(int j=0; j<nrJ(); j++)
    bb.setUnion(joint(j).getBB(TRANSI(t,j)));
  return bb;
}
void ArticulatedBody::writeVTK(const std::string& path,const Mat3XT& t) const {
  VTKWriter<double> os("ArticulatedBody",path,true);
  writeVTK(os,t);
}
void ArticulatedBody::writeVTK(VTKWriter<double>& os,const Mat3XT& t) const {
  for(int JID=0; JID<nrJ(); JID++) {
    Mat3X4T trans,transJ=TRANSI(t,JID);
    APPLY_TRANS(trans,transJ,joint(JID)._transMesh);
    if(joint(JID)._mesh)
      joint(JID)._mesh->writeVTK(os,trans.template cast<ShapeExact::T>());
  }
}
void ArticulatedBody::mimic(VecM xMap) const {
  for(int i=0; i<nrJ(); i++) {
    const Joint& J=joint(i);
    if(J._mimic>=0) {
      const Joint& JM=joint(J._mimic);
      xMap.segment(J._offDOF,J.nrDOF())=xMap.segment(JM._offDOF,JM.nrDOF())*J._mult+Vec::Constant(JM.nrDOF(),J._offset);
    }
  }
}
typename ArticulatedBody::Vec ArticulatedBody::mimic(Vec x) const {
  mimic(mapV(x));
  return x;
}
void ArticulatedBody::mimic(MatT& A,Vec& b,Vec& l,Vec& u) const {
  int nrDOFReduced=0;
  std::unordered_map<int,int> DOFMap;
  for(int i=0; i<nrJ(); i++)
    if(joint(i)._mimic<0) {
      DOFMap[i]=nrDOFReduced;
      nrDOFReduced+=joint(i).nrDOF();
    }
  A.setZero(nrDOF(),nrDOFReduced);
  b.setZero(nrDOF());
  l.setZero(nrDOFReduced);
  u.setZero(nrDOFReduced);
  for(int i=0; i<nrJ(); i++) {
    const Joint& J=joint(i);
    if(J._mimic>=0) {
      int mimic=J._mimic;
      T offset=J._offset;
      T mult=J._mult;
      while(joint(mimic)._mimic>=0) {
        mimic=joint(mimic)._mimic;
        offset+=mult*joint(mimic)._offset;
        mult*=joint(mimic)._mult;
      }
      A.block(J._offDOF,DOFMap[mimic],J.nrDOF(),J.nrDOF())=MatT::Identity(J.nrDOF(),J.nrDOF())*mult;
      b.segment(J._offDOF,J.nrDOF())=Vec::Constant(J.nrDOF(),offset);
    } else {
      A.block(J._offDOF,DOFMap[i],J.nrDOF(),J.nrDOF())=MatT::Identity(J.nrDOF(),J.nrDOF());
      b.segment(J._offDOF,J.nrDOF())=Vec::Zero(J.nrDOF());
      l.segment(DOFMap[i],J.nrDOF())=lowerLimit().segment(J._offDOF,J.nrDOF());
      u.segment(DOFMap[i],J.nrDOF())=upperLimit().segment(J._offDOF,J.nrDOF());
    }
  }
}
bool ArticulatedBody::movable(int jid,int croot) const {
  if(jid==croot)
    return false;
  else if(joint(jid)._parent==-1)
    return joint(jid).nrDOF()>0;
  else if(joint(jid).nrDOF()>0)
    return true;
  else return movable(joint(jid)._parent);
}
const Joint& ArticulatedBody::jointFromDOF(int id) const {
  for(const Joint& J:_joints)
    if(id>=J._offDOF && id-J._offDOF<J.nrDOF())
      return J;
  ASSERT_MSGV(false,"Cannot find joint from DOF %d",id)
  return _joints.front();
}
Joint& ArticulatedBody::jointFromDOF(int id) {
  for(Joint& J:_joints)
    if(id>=J._offDOF && id-J._offDOF<J.nrDOF())
      return J;
  ASSERT_MSGV(false,"Cannot find joint from DOF %d",id)
  return _joints.front();
}
std::shared_ptr<Joint> ArticulatedBody::jointSmartPtr(int id) const {
  return std::shared_ptr<Joint>(new Joint(joint(id)));
}
const Joint& ArticulatedBody::joint(int id) const {
  return _joints[id];
}
Joint& ArticulatedBody::joint(int id) {
  return _joints[id];
}
int ArticulatedBody::jointId(const std::string& name) const {
  for(int k=nrJ()-1; k>=0; k--)
    if(joint(k)._name==name)
      return k;
  ASSERT_MSGV(false,"Cannot file joint with name: %s!",name.c_str())
  return -1;
}
int ArticulatedBody::rootJointId() const {
  for(int k=nrJ()-1; k>=0; k--)
    if(joint(k).isRoot(*this))
      return k;
  ASSERT_MSG(false,"Cannot find root joint!")
  return -1;
}
int ArticulatedBody::depth() const {
  int ret=0;
  for(int i=0; i<nrJ(); i++)
    ret=std::max(ret,joint(i)._depth);
  return ret;
}
int ArticulatedBody::nrDOF() const {
  int ret=0;
  for(int i=0; i<(int)nrJ(); i++)
    ret+=joint(i).nrDOF();
  return ret;
}
int ArticulatedBody::nrDDT() const {
  int ret=0;
  for(int i=0; i<(int)nrJ(); i++)
    ret+=joint(i).nrDDT();
  return ret;
}
int ArticulatedBody::nrJ() const {
  return (int)_joints.size();
}
void ArticulatedBody::fillChildren() {
  for(int i=0; i<(int)nrJ(); i++)
    joint(i)._children.clear();
  for(int i=0; i<(int)nrJ(); i++)
    if(joint(i)._parent>=0)
      joint(joint(i)._parent)._children.push_back(i);
}
void ArticulatedBody::setRootTrans(const Mat3X4T& t) {
  int root=rootJointId(),nDOF=0;
  for(int i=0; i<=root; i++)
    nDOF+=joint(i).nrDOF();
  if(nDOF>0) {
    std::cout<<"Cannot call setT on free-based ArticulatedBody"<<std::endl;
  } else {
    joint(root)._trans=t;
  }
}
void ArticulatedBody::addBase(int dim,const Vec3T& planeNormal) {
  ArticulatedUtils utils(*this);
  utils.addBase(dim,planeNormal);
}
void ArticulatedBody::simplify(int nrDebug) {
  ArticulatedUtils utils(*this);
  utils.simplify(Vec::Zero(nrDOF()),nrDebug);
}
void ArticulatedBody::scaleMass(T coef) {
  ArticulatedUtils(*this).scaleMass(coef);
}
ArticulatedBody::T ArticulatedBody::totalMass() const {
  return ArticulatedUtils(const_cast<ArticulatedBody&>(*this)).totalMass();
}
void ArticulatedBody::clearTransMesh() {
  for(auto &joint:_joints) {
    joint.clearTransMesh();
  }
}
}
