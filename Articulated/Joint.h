#ifndef JOINT_H
#define JOINT_H

#include <Environment/ShapeExact.h>
#include <Utils/SparseUtils.h>
#include <Utils/IO.h>
#include "Environment/MeshExact.h"

namespace PHYSICSMOTION {
struct ArticulatedBody;
struct Joint : public SerializableBase {
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  friend struct ArticulatedBody;
  friend class ArticulatedLoader;
  friend class ArticulatedUtils;
  using SerializableBase::read;
  using SerializableBase::write;
  enum JOINT_TYPE {
    TRANS_3D    =1<<0,
    TRANS_2D    =1<<1,
    TRANS_1D    =1<<2,
    ROT_3D_XYZ  =1<<3,
    ROT_3D_EXP  =1<<4,
    BALL_JOINT  =1<<5,
    HINGE_JOINT =1<<6,
    FIX_JOINT   =1<<7,
    NR_JOINT_TYPE,
  };
  Joint();
  Joint(const Joint& other);
  const Joint& operator=(const Joint& other);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void setType(const std::string& type);
  static std::string typeToString(JOINT_TYPE type);
  BBoxExact getBB(const Mat3X4T& t) const;
  Mat3XT getAxes(bool& markX,bool& markY,bool& markZ,const Mat3XT* t=NULL) const;
  void assemble(T rho);
  void debugTransformMesh();
  void transformMass(const Mat3X4T& trans);
  void transformMesh(const Mat3X4T& trans);
  void transformMesh(const Mat3T& R,const Vec3T& X);
  void transformMesh(const Vec3T& X);
  Eigen::Matrix<int,2,1> CBegEnd() const;
  Eigen::Matrix<int,2,1> RBegEnd() const;
  int nrDOF() const;
  int nrDDT() const;
  Mat6T getMassC(const Mat3T& R) const;
  Mat6T getMassC() const;
  Mat6T getMass(const Mat3T& R) const;
  Mat6T getMass() const;
  Vec3T getC() const;
  template <int R,int C>
  T getMEntry() const {
    if(R==3 && C==3)
      return _M;
    else if(R==3)
      return _MC[C];
    else if(C==3)
      return _MC[R];
    return _MCCT(R,C);
  }
  bool isRotational() const;
  bool isRoot(const ArticulatedBody& body) const;
  std::shared_ptr<ShapeExact> getGeomPtr() const;
  static void loopAllJointTypes(std::function<void(Joint::JOINT_TYPE)> t);
  //Lipschitz constant
  void initL1(const ArticulatedBody& body);
  VecCM getL1(int vertexId) const;
  //joint
  std::vector<int> _children;
  int _parent,_depth,_typeJoint,_mimic;
  int _offDOF,_offDDT;
  int _class;
  Mat3XT _limits;
  Vec _control;
  Vec _damping;
  Mat3X4T _trans;
  //mimic
  T _mult,_offset;
  //mass
  T _M;
  Vec3T _MC;
  Mat3T _MCCT;
  std::string _name;
  //mesh approx.
  std::shared_ptr<ShapeExact> _mesh;
  //Lipschitz constant.
  MatT _L1;
};
}

#endif
