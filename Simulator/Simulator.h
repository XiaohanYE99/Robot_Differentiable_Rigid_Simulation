#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <Articulated/ArticulatedBody.h>
#include <Environment/ContactGenerator.h>

namespace PHYSICSMOTION {
class SoftJoint;
template <typename T>
struct PBDArticulatedGradientInfoMap;
template <typename T>
struct PBDArticulatedGradientInfo;
class Simulator {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef PBDArticulatedGradientInfo<T> GradInfo;
  typedef ContactGenerator::ContactManifold ContactManifold;
  typedef ContactGenerator::ContactPoint ContactPoint;
  typedef std::function<Vec(T,int)> DOFFunction;
  struct PhysicsParameter {
    PhysicsParameter();
    std::array<char,3> _isKinematic;    //whether this joint has specified motion
    std::array<char,3> _isDesign;       //whether this joint is designed parameter
    T _kp,_kd;                          //PD controller parameter
    T _friction;                        //friction coefficient
    T _kc;                              //contact stiffness
    DOFFunction _kin,_tarP,_tarD;
  };
  struct XPBDConstraint {
    void update(VecM x);
    T _C=0;
    T _alpha=0;
    T _lambda=0;
    T _dLambda=0;
    Vec _JC,_invMJC;
    //case with linear constraint
    bool _isLinear=false;
    T _CL=std::numeric_limits<double>::infinity();
    T _CU=std::numeric_limits<double>::infinity();
  };
  //API
  Simulator(T dt);
  virtual ~Simulator();
  virtual void clearJoint();
  virtual void addJoint(const SoftJoint& joint);
  virtual const std::vector<SoftJoint>& getJoints() const;
  virtual std::vector<SoftJoint>& getJoints();
  virtual void clearShape();
  virtual void addShape(std::shared_ptr<ShapeExact> shape);
  virtual void setArticulatedBody(std::shared_ptr<ArticulatedBody> body);
  virtual void setHeuristcGuessStiffness(T coef=1e1f);
  virtual T getHeuristcGuessStiffness() const;
  virtual const std::vector<ContactManifold>& getManifolds() const;
  virtual const std::vector<std::shared_ptr<ShapeExact>>& getShapes() const;
  virtual const PhysicsParameter& getJointPhysicsParameter(int jid) const;
  virtual PhysicsParameter& getJointPhysicsParameter(int jid);
  virtual std::shared_ptr<ArticulatedBody> getBody() const;
  virtual std::shared_ptr<ContactGenerator> getContactSmartPtr() const;
  virtual ContactGenerator& getContact();
  virtual void setGravity(const Vec3T& g);
  virtual Vec3T getGravity() const;
  virtual void setDesign(const Vec& design);
  virtual Vec getDesign() const;
  virtual void setTime(T t);
  virtual T getTime() const;
  virtual void setOutput(bool output);
  virtual void step()=0;
  virtual Vec pos() const=0;
  virtual void setPos(const Vec& pos)=0;
  virtual Vec vel() const=0;
  virtual void setVel(const Vec& vel)=0;
  virtual void detectCurrentContact()=0;
 protected:
  virtual Vec setKinematic(const Vec& pos,T time) const;
  virtual Vec setKinematicVel(const Vec& vel,T time,T dt) const;
  virtual void computeLocalContactPos(const Mat3XT& t);
  virtual void detectContact(const Mat3XT& t);
  void mask(MatT* diag,Vec* DE,MatT* DDE,MatT* DDEX=NULL) const;
  void maskNonDesign(MatT& DDE) const;
  void markDesign(MatT& DDE) const;
  T energyInertial(const GradInfo& next,const GradInfo& curr,const GradInfo& last,
                   int jid,T damping,T M,const Vec3T& P,const Mat3T& PPT,
                   Mat3XT* G,Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt) const;
  T energyPDController(VecCM x,VecCM xL,int jid,const Joint& J,int nrD,Vec* DE,MatT* DDE) const;
  //data
  T _dt,_t;
  Mat3XT _JRCF;
  Vec _P,_D,_design;
  bool _output;
  std::vector<PhysicsParameter> _params;
  std::vector<std::shared_ptr<ShapeExact>> _shapes;
  std::shared_ptr<ArticulatedBody> _body;
  std::shared_ptr<ContactGenerator> _contact;
  std::vector<ContactManifold> _manifolds;
  std::vector<SoftJoint> _joints;
  //empty function
  static DOFFunction _defaultFunction;
};
}

#endif
