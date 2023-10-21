#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <Articulated/ArticulatedBody.h>
#include <Environment/ContactGenerator.h>

namespace PHYSICSMOTION {
class Simulator {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef ContactGenerator::ContactManifold ContactManifold;
  typedef ContactGenerator::ContactPoint ContactPoint;
  typedef std::function<Vec(T,int)> DOFFunction;
  struct PhysicsParameter {
    PhysicsParameter();
    bool _isKinematic;  //whether this joint has specified motion
    T _kp,_kd;          //PD controller parameter
    T _friction;        //friction coefficient
    T _kc;              //contact stiffness
    DOFFunction _kin,_tarP,_tarD;
  };
  struct DragEnergy {
    DragEnergy();
    Vec3T _ptL,_pt,_lambda;
    int _jid;
    T _k;
  };
  //API
  Simulator(T dt);
  virtual ~Simulator();
  virtual void clearDrag();
  virtual void addDrag(const DragEnergy& drag);
  virtual const std::vector<DragEnergy>& getDrags() const;
  virtual std::vector<DragEnergy>& getDrags();
  virtual void clearShape();
  virtual void addShape(std::shared_ptr<ShapeExact> shape);
  virtual void setArticulatedBody(std::shared_ptr<ArticulatedBody> body);
  virtual const std::vector<ContactManifold>& getManifolds() const;
  virtual const std::vector<std::shared_ptr<ShapeExact>>& getShapes() const;
  virtual const PhysicsParameter& getJointPhysicsParameter(int jid) const;
  virtual PhysicsParameter& getJointPhysicsParameter(int jid);
  virtual std::shared_ptr<ArticulatedBody> getBody() const;
  virtual ContactGenerator& getContact();
  virtual void setTime(T t);
  virtual T getTime() const;
  virtual void step()=0;
  virtual Vec pos() const=0;
  virtual void setPos(const Vec& pos)=0;
  virtual Vec vel() const=0;
  virtual void setVel(const Vec& vel)=0;
  virtual void setGravity(const Vec3T& g)=0;
  virtual void detectCurrentContact()=0;
 protected:
  virtual Vec setKinematic(const Vec& pos,T time) const;
  virtual Vec setKinematicVel(const Vec& vel,T time,T dt) const;
  virtual void computeLocalContactPos(const Mat3XT& t);
  virtual void detectContact(const Mat3XT& t);
  //data
  T _dt,_t;
  std::vector<PhysicsParameter> _params;
  std::vector<std::shared_ptr<ShapeExact>> _shapes;
  std::shared_ptr<ArticulatedBody> _body;
  std::shared_ptr<ContactGenerator> _contact;
  std::vector<ContactManifold> _manifolds;
  std::vector<DragEnergy> _drags;
  //empty function
  static DOFFunction _defaultFunction;
};
}

#endif
