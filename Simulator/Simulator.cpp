#include "Simulator.h"
#include "SoftJoint.h"
#include <Utils/CrossSpatialUtils.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
Simulator::PhysicsParameter::PhysicsParameter()
  :_isKinematic({false,false,false}),_isDesign({false,false,false}),_kp(0),_kd(0),_friction(.75f),_kc(1e1f) {
  _kin=_tarP=_tarD=Simulator::_defaultFunction;
}
void Simulator::XPBDConstraint::update(VecM x) {
  ASSERT_MSG(_isLinear,"Cannot update nonlinear XPBDConstraint!")
  T C=x.dot(_JC);
  if(isfinite(_CL) && C<_CL)
    _C=C-_CL;
  else if(isfinite(_CU) && C>_CU)
    _C=C-_CU;
  else _C=0;
}
Simulator::Simulator(T dt):_dt(dt),_t(0),_output(false) {}
Simulator::~Simulator() {}
void Simulator::clearJoint() {
  _joints.clear();
}
void Simulator::addJoint(const SoftJoint& joint) {
  _joints.push_back(joint);
}
const std::vector<SoftJoint>& Simulator::getJoints() const {
  return _joints;
}
std::vector<SoftJoint>& Simulator::getJoints() {
  return _joints;
}
void Simulator::clearShape() {
  if(_body)
    _params.assign(_body->nrJ(),PhysicsParameter());
  else _params.clear();
  _contact=NULL;
  _t=0;
}
void Simulator::addShape(std::shared_ptr<ShapeExact> shape) {
  _shapes.push_back(shape);
  _contact=NULL;
}
void Simulator::setArticulatedBody(std::shared_ptr<ArticulatedBody> body) {
  _body=body;
  if(_body)
    _params.assign(_body->nrJ(),PhysicsParameter());
  else _params.clear();
  _contact=NULL;
}
void Simulator::setHeuristcGuessStiffness(T coef) {
  for(int k=0; k<_body->nrJ(); k++) {
    T Mg=0;
    int parentId=k;
    while(parentId>=0) {
      Mg+=CTRI(_JRCF,parentId).norm();
      parentId=_body->joint(parentId)._parent;
    }
    _params[k]._kc=coef*Mg;
  }
}
Simulator::T Simulator::getHeuristcGuessStiffness() const {
  T num=0,denom=0;
  for(int k=0; k<_body->nrJ(); k++) {
    T Mg=0;
    int parentId=k;
    while(parentId>=0) {
      Mg+=CTRI(_JRCF,parentId).norm();
      parentId=_body->joint(parentId)._parent;
    }
    num+=_params[k]._kc;
    denom+=Mg;
  }
  return num/denom;
}
const std::vector<Simulator::ContactManifold>& Simulator::getManifolds() const {
  return _manifolds;
}
const std::vector<std::shared_ptr<ShapeExact>>& Simulator::getShapes() const {
  return _shapes;
}
const Simulator::PhysicsParameter& Simulator::getJointPhysicsParameter(int jid) const {
  return _params[jid];
}
Simulator::PhysicsParameter& Simulator::getJointPhysicsParameter(int jid) {
  return _params[jid];
}
std::shared_ptr<ArticulatedBody> Simulator::getBody() const {
  return _body;
}
std::shared_ptr<ContactGenerator> Simulator::getContactSmartPtr() const {
  return _contact;
}
ContactGenerator& Simulator::getContact() {
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  PBDArticulatedGradientInfo<T> info(*_body,pos());
  ContactGenerator::Mat3XT t=info._TM.template cast<GEOMETRY_SCALAR>();
  //_contact->updateBVH(t);
  return *_contact;
}
void Simulator::setGravity(const Vec3T& g) {
  _JRCF.setZero(3,_body->nrJ()*4);
  for(int i=0; i<_body->nrJ(); i++)
    if(_body->joint(i)._M>0) {
      ROTI(_JRCF,i)=-g*_body->joint(i)._MC.transpose().template cast<T>();
      CTRI(_JRCF,i)=-_body->joint(i)._M*g;
    }
}
Simulator::Vec3T Simulator::getGravity() const {
  T denom=0;
  Vec3T g=Vec3T::Zero();
  for(int i=0; i<_body->nrJ(); i++)
    if(_body->joint(i)._M>0) {
      g-=CTRI(_JRCF,i);
      denom+=_body->joint(i)._M;
    }
  return g/denom;
}
void Simulator::setDesign(const Vec& design) {
  _design=design;
}
Simulator::Vec Simulator::getDesign() const {
  return _design;
}
void Simulator::setTime(T t) {
  _t=t;
  setPos(setKinematic(pos(),_t));
}
Simulator::T Simulator::getTime() const {
  return _t;
}
void Simulator::setOutput(bool output) {
  _output=output;
}
//helper
Simulator::Vec Simulator::setKinematic(const Vec& pos,T time) const {
  Vec posKinematic=pos;
  for(int i=0; i<_body->nrJ(); i++) {
    const Joint& J=_body->joint(i);
    int nrDJ=J.nrDOF();
    Vec curr=_params[i]._kin(time,nrDJ);
    for(int c=0; c<nrDJ; c++) {
      //kinematic target
      if(_params[i]._isKinematic[c])
        posKinematic[J._offDOF+c]=curr[c];
      //design
      if(_params[i]._isDesign[c])
        posKinematic[J._offDOF+c]=_design[J._offDOF+c];
      //constraint
      T l=J._limits(0,c);
      T h=J._limits(1,c);
      bool locked=isfinite(J._limits(2,c)) && J._limits(2,c)>0 && isfinite(l) && isfinite(h) && l==h;
      if(locked)
        posKinematic[J._offDOF+c]=l;
    }
  }
  return posKinematic;
}
Simulator::Vec Simulator::setKinematicVel(const Vec& vel,T time,T dt) const {
  Vec velKinematic=vel;
  for(int i=0; i<_body->nrJ(); i++) {
    int nrDJ=_body->joint(i).nrDOF();
    Vec next=_params[i]._kin(time+dt,nrDJ);
    Vec curr=_params[i]._kin(time,nrDJ);
    for(int c=0; c<nrDJ; c++)
      if(_params[i]._isKinematic[c])
        velKinematic[_body->joint(i)._offDOF+c]=(next[c]-curr[c])/dt;
  }
  return velKinematic;
}
void Simulator::computeLocalContactPos(const Mat3XT& t) {
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      ContactGenerator::Vec3T::Index id;
      p._nA2B.cwiseAbs().minCoeff(&id);
      p._tA2B.col(0)=p._nA2B.cross(ContactGenerator::Vec3T::Unit(id)).template cast<T>().normalized().template cast<ContactGenerator::T>();
      p._tA2B.col(1)=p._nA2B.cross(p._tA2B.col(0));
      if(m._jidA>=0)
        p._ptAL=ROTI(t,m._jidA).template cast<GEOMETRY_SCALAR>().transpose()*(p._ptA-CTRI(t,m._jidA).template cast<GEOMETRY_SCALAR>());
      if(m._jidB>=0)
        p._ptBL=ROTI(t,m._jidB).template cast<GEOMETRY_SCALAR>().transpose()*(p._ptB-CTRI(t,m._jidB).template cast<GEOMETRY_SCALAR>());
    }
}
void Simulator::detectContact(const Mat3XT& t) {
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  _contact->generateManifolds(0,false,_manifolds,t.template cast<GEOMETRY_SCALAR>());
  computeLocalContactPos(t);
}
void Simulator::mask(MatT* diag,Vec* DE,MatT* DDE,MatT* DDEX) const {
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    int nrDJ=J.nrDOF();
    for(int c=0; c<nrDJ; c++) {
      T l=J._limits(0,c);
      T h=J._limits(1,c);
      bool locked=isfinite(J._limits(2,c)) && J._limits(2,c)>0 && isfinite(l) && isfinite(h) && l==h;
      if(locked || _params[k]._isKinematic[c] || _params[k]._isDesign[c]) {
        if(diag)
          diag->diagonal()[J._offDOF+c]=std::numeric_limits<double>::infinity();
        if(DE)
          DE->coeffRef(J._offDOF+c)=0;
        if(DDE) {
          DDE->row(J._offDOF+c).setZero();
          DDE->col(J._offDOF+c).setZero();
          DDE->diagonal()[J._offDOF+c]=0;
        }
        if(DDEX)
          DDEX->row(J._offDOF+c).setZero();
      }
    }
  }
}
void Simulator::maskNonDesign(MatT& DDE) const {
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    int nrDJ=J.nrDOF();
    for(int c=0; c<nrDJ; c++) {
      T l=J._limits(0,c);
      T h=J._limits(1,c);
      bool locked=isfinite(J._limits(2,c)) && J._limits(2,c)>0 && isfinite(l) && isfinite(h) && l==h;
      if(locked || _params[k]._isKinematic[c] || _params[k]._isDesign[c])
        DDE.row(J._offDOF+c).setZero();
      if(!_params[k]._isDesign[c])
        DDE.col(J._offDOF+c).setZero();
    }
  }
}
void Simulator::markDesign(MatT& DDE) const {
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    int nrDJ=J.nrDOF();
    for(int c=0; c<nrDJ; c++) {
      T l=J._limits(0,c);
      T h=J._limits(1,c);
      bool locked=isfinite(J._limits(2,c)) && J._limits(2,c)>0 && isfinite(l) && isfinite(h) && l==h;
      if(locked || _params[k]._isKinematic[c])
        continue;
      if(_params[k]._isDesign[c])
        DDE(J._offDOF+c,J._offDOF+c)=1;
    }
  }
}
Simulator::T Simulator::energyInertial
(const GradInfo& next,const GradInfo& curr,const GradInfo& last,
 int jid,T damping,T M,const Vec3T& P,const Mat3T& PPT,
 Mat3XT* G,Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt) const {
  T E=0,coef;
  Mat3XT A;
  //dynamic force
  coef=1/(_dt*_dt);
  if(G) {
    MRR.template block<3,3>(0,jid*3)-=invDoubleCrossMatTrace<T>(ROTI(next._TM,jid)*PPT*ROTI(next._TM,jid).transpose())*coef;
    MRt.template block<3,3>(0,jid*3)+=cross<T>(ROTI(next._TM,jid)*P)*coef;
    MtR.template block<3,3>(0,jid*3)-=cross<T>(ROTI(next._TM,jid)*P)*coef;
    Mtt.template block<3,3>(0,jid*3)+=Mat3T::Identity()*M*coef;
  }
  A=TRANSI(next._TM,jid)-2*TRANSI(curr._TM,jid)+TRANSI(last._TM,jid);
  E+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*M).trace()*coef/2;
  if(G) {
    ROTI((*G),jid)+=(ROT(A)*PPT+CTR(A)*P.transpose())*coef;
    CTRI((*G),jid)+=(CTR(A)*M+ROT(A)*P)*coef;
  }
  //damping force
  if(damping>0) {
    coef=damping/_dt;
    if(G) {
      MRR.template block<3,3>(0,jid*3)-=invDoubleCrossMatTrace<T>(ROTI(next._TM,jid)*PPT*ROTI(next._TM,jid).transpose())*coef;
      MRt.template block<3,3>(0,jid*3)+=cross<T>(ROTI(next._TM,jid)*P)*coef;
      MtR.template block<3,3>(0,jid*3)-=cross<T>(ROTI(next._TM,jid)*P)*coef;
      Mtt.template block<3,3>(0,jid*3)+=Mat3T::Identity()*M*coef;
    }
    A=TRANSI(next._TM,jid)-TRANSI(curr._TM,jid);
    E+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*M).trace()*coef/2;
    if(G) {
      ROTI((*G),jid)+=(ROT(A)*PPT+CTR(A)*P.transpose())*coef;
      CTRI((*G),jid)+=(CTR(A)*M+ROT(A)*P)*coef;
    }
  }
  //external force
  E+=(TRANSI(next._TM,jid)*TRANSI(_JRCF,jid).transpose()).trace();
  if(G)
    TRANSI((*G),jid)+=TRANSI(_JRCF,jid);
  return E;
}
Simulator::T Simulator::energyPDController(VecCM x,VecCM xL,int jid,const Joint& J,int nrD,Vec* DE,MatT* DDE) const {
  T E=0;
  Eigen::Matrix<T,-1,1,0,3,1> diff,coef;
  if(J._control.size()>=nrD)
    coef=J._control.segment(0,nrD).template cast<T>();
  else coef.setOnes(nrD);
  //P controller
  if(_params[jid]._kp>0) {
    diff=x.segment(J._offDOF,nrD)-_params[jid]._tarP(_t,nrD);
    E+=(diff.array().square()*coef.array()*_params[jid]._kp).sum()/2;
    if(DE)
      DE->segment(J._offDOF,nrD).array()+=diff.array()*coef.array()*_params[jid]._kp;
    if(DDE)
      DDE->diagonal().segment(J._offDOF,nrD).array()+=coef.array()*_params[jid]._kp;
  }
  //D controller
  if(_params[jid]._kd>0) {
    diff=(x-xL).segment(J._offDOF,nrD)/_dt-_params[jid]._tarD(_t,nrD);
    E+=(diff.array().square()*coef.array()*_params[jid]._kd).sum()/2;
    if(DE)
      DE->segment(J._offDOF,nrD).array()+=diff.array()*coef.array()*_params[jid]._kd/_dt;
    if(DDE)
      DDE->diagonal().segment(J._offDOF,nrD).array()+=coef.array()*_params[jid]._kd/_dt/_dt;
  }
  return E;
}
Simulator::DOFFunction Simulator::_defaultFunction=[](Simulator::T,int nrDOF) {
  return Simulator::Vec::Zero(nrDOF);
};
}
