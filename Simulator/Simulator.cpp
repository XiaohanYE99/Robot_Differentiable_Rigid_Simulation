#include "Simulator.h"
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
Simulator::PhysicsParameter::PhysicsParameter()
  :_isKinematic(false),_kp(0),_kd(0),_friction(.75f),_kc(1e2f) {
  _kin=_tarP=_tarD=Simulator::_defaultFunction;
}
Simulator::DragEnergy::DragEnergy():_jid(-1),_k(1e2f) {}
Simulator::Simulator(T dt):_dt(dt),_t(0) {}
Simulator::~Simulator() {}
void Simulator::clearDrag() {
  _drags.clear();
}
void Simulator::addDrag(const DragEnergy& drag) {
  _drags.push_back(drag);
}
const std::vector<Simulator::DragEnergy>& Simulator::getDrags() const {
  return _drags;
}
std::vector<Simulator::DragEnergy>& Simulator::getDrags() {
  return _drags;
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
  ArticulatedUtils(*_body).tessellate(true); //tessellate the mesh
  if(_body)
    _params.assign(_body->nrJ(),PhysicsParameter());
  else _params.clear();
  _contact=NULL;
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
ContactGenerator& Simulator::getContact() {
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  PBDArticulatedGradientInfo<T> info(*_body,pos());
  ContactGenerator::Mat3XT t=info._TM.template cast<GEOMETRY_SCALAR>();
  _contact->updateBVH(t);
  return *_contact;
}
void Simulator::setTime(T t) {
  _t=t;
  setPos(setKinematic(pos(),_t));
}
Simulator::T Simulator::getTime() const {
  return _t;
}
//helper
Simulator::Vec Simulator::setKinematic(const Vec& pos,T time) const {
  Vec posKinematic=pos;
  for(int i=0; i<_body->nrJ(); i++)
    if(_params[i]._isKinematic) {
      int nrJ=_body->joint(i).nrDOF();
      posKinematic.segment(_body->joint(i)._offDOF,nrJ)=_params[i]._kin(time,nrJ);
    }
  return posKinematic;
}
Simulator::Vec Simulator::setKinematicVel(const Vec& vel,T time,T dt) const {
  Vec velKinematic=vel;
  for(int i=0; i<_body->nrJ(); i++)
    if(_params[i]._isKinematic) {
      int nrJ=_body->joint(i).nrDOF();
      velKinematic.segment(_body->joint(i)._offDOF,nrJ)=(_params[i]._kin(time+dt,nrJ)-_params[i]._kin(time,nrJ))/dt;
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
      p._lambda.setZero();
    }
}
void Simulator::detectContact(const Mat3XT& t) {
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  _contact->generateManifolds(_manifolds,t.template cast<GEOMETRY_SCALAR>());
  computeLocalContactPos(t);
}
Simulator::DOFFunction Simulator::_defaultFunction=[](Simulator::T,int nrDOF) {
  return Simulator::Vec::Zero(nrDOF);
};
}
