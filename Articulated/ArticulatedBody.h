#ifndef ARTICULATED_BODY_H
#define ARTICULATED_BODY_H

#include <set>
#include "Joint.h"
#include <tinyxml2.h>

namespace PHYSICSMOTION {
struct ArticulatedBody : public SerializableBase {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  using SerializableBase::read;
  using SerializableBase::write;
  friend class ArticulatedLoader;
  friend class ArticulatedUtils;
  ArticulatedBody();
  ArticulatedBody(const tinyxml2::XMLElement& pt);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void randomize(int nrLink,bool chain=true);
  std::set<int> children(int id,bool direct=false) const;
  int commonRoot(int id,int id2) const;
  int hasJoint(Joint::JOINT_TYPE type) const;
  ArticulatedBody resizeJoints(int nr) const;
  bool isLeaf(int id) const;
  Vec control() const;
  Vec damping() const;
  Vec coefLimit() const;
  Vec lowerLimit() const;
  Vec upperLimit() const;
  Vec lowerLimit(T infty) const;
  Vec upperLimit(T infty) const;
  const Vec& clampLimit(Vec& x) const;
  Vec randomPose(T coef) const;
  Mat3XT getT(const Vec& x) const;
  BBoxExact getBB(const Vec& x) const;
  BBoxExact getBB(const Mat3XT& t) const;
  void writeVTK(const std::string& path,const Mat3XT& t) const;
  void writeVTK(VTKWriter<double>& os,const Mat3XT& t) const;
  void mimic(VecM xMap) const;
  Vec mimic(Vec x) const;
  void mimic(MatT& A,Vec& b,Vec& l,Vec& u) const;
  bool movable(int id,int croot=-1) const;
  const Joint& jointFromDOF(int id) const;
  Joint& jointFromDOF(int id);
  std::shared_ptr<Joint> jointSmartPtr(int id) const;
  const Joint& joint(int id) const;
  Joint& joint(int id);
  int jointId(const std::string& name) const;
  int rootJointId() const;
  int depth() const;
  int nrDOF() const;
  int nrDDT() const;
  int nrJ() const;
  void fillChildren();
  void setRootTrans(const Mat3X4T& t);
  void addBase(int dim,const Vec3T& planeNormal);
  void simplify(int nrDebug);
  void scaleMass(T coef);
  T totalMass() const;
 protected:
  std::vector<Joint> _joints;
};
}

#endif
