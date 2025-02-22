#include "Joint.h"
#include "RigidBodyMass.h"
#include "ArticulatedBody.h"
#include <Environment/BBoxExact.h>
#include <Environment/MeshExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/CompositeShapeExact.h>
#include <Utils/CrossSpatialUtils.h>
#include <Utils/RotationUtils.h>

namespace PHYSICSMOTION {
//Joint
Joint::Joint():_parent(-1),_mimic(-1) {}
Joint::Joint(const Joint& other) {
  operator=(other);
}
const Joint& Joint::operator=(const Joint& other) {
  //joint
  _children=other._children;
  _parent=other._parent;
  _class=other._class;
  _depth=other._depth;
  _typeJoint=other._typeJoint;
  _mimic=other._mimic;
  _offDOF=other._offDOF;
  _offDDT=other._offDDT;
  _limits=other._limits;
  _control=other._control;
  _damping=other._damping;
  _trans=other._trans;
  //mimic
  _mult=other._mult;
  _offset=other._offset;
  //mass
  _M=other._M;
  _MC=other._MC;
  _MCCT=other._MCCT;
  _name=other._name;
  //mesh approx.
  if(other._mesh)
    _mesh=std::dynamic_pointer_cast<ShapeExact>(other._mesh->copy());
  else _mesh=NULL;
  _L1=other._L1;
  return *this;
}
bool Joint::read(std::istream& is,IOData* dat) {
  registerType<MeshExact>(dat);
  registerType<ConvexHullExact>(dat);
  registerType<BBoxExact>(dat);
  registerType<SphericalBBoxExact>(dat);
  registerType<CompositeShapeExact>(dat);
  //basic
  readBinaryData(_children,is);
  readBinaryData(_parent,is);
  readBinaryData(_depth,is);
  readBinaryData(_typeJoint,is);
  readBinaryData(_mimic,is);
  readBinaryData(_offDOF,is);
  readBinaryData(_offDDT,is);
  readBinaryData(_limits,is);
  readBinaryData(_control,is);
  readBinaryData(_damping,is);
  readBinaryData(_trans,is);
  //mimic
  readBinaryData(_mult,is);
  readBinaryData(_offset,is);
  //mass
  readBinaryData(_M,is);
  readBinaryData(_MC,is);
  readBinaryData(_MCCT,is);
  readBinaryData(_name,is);
  //mesh approx.
  readBinaryData(_mesh,is,dat);
  //L1
  readBinaryData(_L1,is);
  return is.good();
}
bool Joint::write(std::ostream& os,IOData* dat) const {
  registerType<MeshExact>(dat);
  registerType<ConvexHullExact>(dat);
  registerType<BBoxExact>(dat);
  registerType<SphericalBBoxExact>(dat);
  registerType<CompositeShapeExact>(dat);
  //basic
  writeBinaryData(_children,os);
  writeBinaryData(_parent,os);
  writeBinaryData(_depth,os);
  writeBinaryData(_typeJoint,os);
  writeBinaryData(_mimic,os);
  writeBinaryData(_offDOF,os);
  writeBinaryData(_offDDT,os);
  writeBinaryData(_limits,os);
  writeBinaryData(_control,os);
  writeBinaryData(_damping,os);
  writeBinaryData(_trans,os);
  //mimic
  writeBinaryData(_mult,os);
  writeBinaryData(_offset,os);
  //mass
  writeBinaryData(_M,os);
  writeBinaryData(_MC,os);
  writeBinaryData(_MCCT,os);
  writeBinaryData(_name,os);
  //mesh approx.
  writeBinaryData(_mesh,os,dat);
  //L1
  writeBinaryData(_L1,os);
  return os.good();
}
std::shared_ptr<SerializableBase> Joint::copy() const {
  return std::shared_ptr<SerializableBase>(new Joint);
}
std::string Joint::type() const {
  return typeid(Joint).name();
}
void Joint::setType(const std::string& type) {
#define S2T(NAME)if(type==#NAME){_typeJoint=NAME;return;}
  S2T(TRANS_3D)
  S2T(TRANS_2D)
  S2T(TRANS_1D)
  S2T(ROT_3D_XYZ)
  S2T(ROT_3D_EXP)
  S2T(BALL_JOINT)
  S2T(HINGE_JOINT)
  S2T(FIX_JOINT)
  ASSERT_MSGV(false,"Unknown Type: %s",type.c_str())
#undef S2T
}
std::string Joint::typeToString(JOINT_TYPE type) {
#define T2S(NAME)if(type==NAME){return #NAME;}
  T2S(TRANS_3D)
  T2S(TRANS_2D)
  T2S(TRANS_1D)
  T2S(ROT_3D_XYZ)
  T2S(ROT_3D_EXP)
  T2S(BALL_JOINT)
  T2S(HINGE_JOINT)
  T2S(FIX_JOINT)
  ASSERT_MSGV(false,"Unknown Type: %d",type)
  return "";
#undef T2S
}
BBoxExact Joint::getBB(const Mat3X4T& t) const {
  /*std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  if(_mesh)
    _mesh->getMesh(vss,iss);
  
  Vec3T center(0.,0.,0.);
  for(const Eigen::Matrix<double,3,1>& v:vss)
    center+=(ROT(t)*v.template cast<T>()+CTR(t)).template cast<BBoxExact::T>();
  center/=vss.size();
  T rad=-1;
  for(const Eigen::Matrix<double,3,1>& v:vss)
    rad=std::max(rad,((ROT(t)*v.template cast<T>()+CTR(t)).template cast<BBoxExact::T>()-center).norm());
  BBoxExact bb(center,rad);
  return bb;*/
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  if(_mesh)
    _mesh->getMesh(vss,iss);
  BBoxExact bb;
  for(const Eigen::Matrix<double,3,1>& v:vss)
    bb.setUnion((ROT(t)*v.template cast<T>()+CTR(t)).template cast<BBoxExact::T>());
  return bb;
}
Joint::Mat3XT Joint::getAxes(bool& markX,bool& markY,bool& markZ,const Mat3XT* t) const {
  Mat3XT ret=Mat3T::Identity();
  if(_typeJoint==TRANS_3D) {
    markX=true;
    markY=true;
    markZ=true;
  } else if(_typeJoint==TRANS_2D) {
    ret=ret.template block<3,2>(0,0);
    markX=true;
    markY=true;
    markZ=false;
  } else if(_typeJoint==TRANS_1D) {
    ret=ret.template block<3,1>(0,0);
    markX=true;
    markY=false;
    markZ=false;
  } else if(_typeJoint==ROT_3D_EXP || _typeJoint==ROT_3D_XYZ) {
    markX=true;
    markY=true;
    markZ=true;
  } else if(_typeJoint==BALL_JOINT) {
    ret=ret.template block<3,2>(0,1);
    markX=false;
    markY=true;
    markZ=true;
  } else if(_typeJoint==HINGE_JOINT) {
    ret=ret.template block<3,1>(0,2);
    markX=false;
    markY=false;
    markZ=true;
  } else if(_typeJoint==FIX_JOINT) {
    ret.setZero(3,0);
    markX=false;
    markY=false;
    markZ=false;
  } else {
    ASSERT_MSGV(false,"Unknown joint type in %s,name=%s,_typeJoint=%d!",__FUNCTION__,_name.c_str(),_typeJoint)
  }
  if(t) {
    Mat3X4T TJ,TJP=TRANSI((*t),_parent);
    ret=ROT(_trans)*ret;
    if(_parent>=0)
      ret=ROT(TJP)*ret;
    //center
    APPLY_TRANS(TJ,TJP,_trans);
    ret=concatCol<Mat3XT>(CTR(TJ),ret);
  }
  return ret;
}
void Joint::assemble(T rho) {
  if(!_mesh) {
    _M=0;
    _MC.setZero();
    _MCCT.setZero();
  } else {
    //fill mass
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    _mesh->getMesh(vss,iss);
    RigidBodyMass<T> M(vss,iss);
    _M=(T) M.getM();
    _MC=M.getMC().cast<T>();
    _MCCT=M.getMCCT().cast<T>();
    Eigen::SelfAdjointEigenSolver<Mat6T> eig(M.getMassCOM().cast<T>());
    ASSERT_MSG(eig.eigenvalues().minCoeff()>=0,"Min Eig is negative")
    ASSERT_MSG(!std::isinf(_M) && !std::isnan(_M),"Infinite Mass")
    ASSERT_MSG(_MC.array().isFinite().all(),"MC is infinite")
    ASSERT_MSG(_MCCT.array().isFinite().all(),"_MCCT is infinite")
    //multiply by rho
    _M*=rho;
    _MC*=rho;
    _MCCT*=rho;
  }
}
void Joint::debugTransformMesh() {
  assemble(1);
  Mat3X4T t;
  ROT(t)=expWGradV<T,Vec3T>(Vec3T::Random()*M_PI);
  CTR(t)=Vec3T::Random();
  transformMesh(t);
  Mat6T m0=getMass();
  assemble(1);
  DEFINE_NUMERIC_DELTA_T(T)
  DEBUG_GRADIENT("Transformed-Mass",getMass().norm(),(m0-getMass()).norm())
}
void Joint::transformMass(const Mat3X4T& trans) {
  Mat3T tmp=CTR(trans)*_MC.transpose()*ROT(trans).transpose();
  _MCCT=ROT(trans)*_MCCT*ROT(trans).transpose()+tmp+tmp.transpose()+_M*CTR(trans)*CTR(trans).transpose();
  _MC=ROT(trans)*_MC+CTR(trans)*_M;
}
void Joint::transformMesh(const Mat3X4T& trans) {
  if(!_mesh)
    return;
  else if(std::dynamic_pointer_cast<CompositeShapeExact>(_mesh))
    std::dynamic_pointer_cast<CompositeShapeExact>(_mesh)->transform(trans.template cast<GEOMETRY_SCALAR>());
  else _mesh.reset(new CompositeShapeExact({_mesh}, {trans.template cast<GEOMETRY_SCALAR>()}));
  transformMass(trans);
}
void Joint::transformMesh(const Mat3T& R,const Vec3T& X) {
  transformMesh((Mat3X4T) concatCol(R,X));
}
void Joint::transformMesh(const Vec3T& X) {
  transformMesh(Mat3T::Identity(),X);
}
Eigen::Matrix<int,2,1> Joint::CBegEnd() const {
  if(_typeJoint==TRANS_3D)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+3);
  else if(_typeJoint==TRANS_2D)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+2);
  else if(_typeJoint==TRANS_1D)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+1);
  else return Eigen::Matrix<int,2,1>(0,0);
}
Eigen::Matrix<int,2,1> Joint::RBegEnd() const {
  if(_typeJoint==ROT_3D_EXP || _typeJoint==ROT_3D_XYZ)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+3);
  else if(_typeJoint==BALL_JOINT)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+2);
  else if(_typeJoint==HINGE_JOINT)
    return Eigen::Matrix<int,2,1>(_offDOF,_offDOF+1);
  else return Eigen::Matrix<int,2,1>(0,0);
}
int Joint::nrDOF() const {
  if(_typeJoint==TRANS_3D)
    return 3;
  else if(_typeJoint==TRANS_2D)
    return 2;
  else if(_typeJoint==TRANS_1D)
    return 1;
  else if(_typeJoint==ROT_3D_EXP || _typeJoint==ROT_3D_XYZ)
    return 3;
  else if(_typeJoint==BALL_JOINT)
    return 2;
  else if(_typeJoint==HINGE_JOINT)
    return 1;
  else if(_typeJoint==FIX_JOINT)
    return 0;
  else {
    ASSERT_MSGV(false,"Unknown joint type in %s,name=%s,_typeJoint=%d!",__FUNCTION__,_name.c_str(),_typeJoint)
    return -1;
  }
}
int Joint::nrDDT() const {
  if(_typeJoint==ROT_3D_EXP || _typeJoint==ROT_3D_XYZ)
    return 9;
  else if(_typeJoint==BALL_JOINT)
    return 1;
  else return 0;
}
Joint::Mat6T Joint::getMassC(const Mat3T& R) const {
  Mat6T ret=getMassC();
  ret.block<3,3>(3,3)=R*(ret.block<3,3>(3,3)*R.transpose()).eval();
  return ret;
}
Joint::Mat6T Joint::getMassC() const {
  Mat6T ret=Mat6T::Zero();
  if(_M<=0)
    return ret;

  Mat6T M=getMass();
  Vec3T C=getC();
  Mat3T T=C*C.transpose();
  T.diagonal().array()-=C.squaredNorm();
  ret.block<3,3>(0,0)=M.block<3,3>(0,0);
  ret.block<3,3>(3,3)=M.block<3,3>(3,3)+T*_M;
  return ret;
}
Joint::Mat6T Joint::getMass(const Mat3T& R) const {
  Mat6T ret=getMass();
  ret.block<3,3>(3,3)=R*(ret.block<3,3>(3,3)*R.transpose()).eval();
  ret.block<3,3>(0,3)=R*(ret.block<3,3>(0,3)*R.transpose()).eval();
  ret.block<3,3>(3,0)=R*(ret.block<3,3>(3,0)*R.transpose()).eval();
  ret.block<3,3>(0,0)=R*(ret.block<3,3>(0,0)*R.transpose()).eval();
  return ret;
}
Joint::Mat6T Joint::getMass() const {
  Mat6T ret=Mat6T::Zero();
  ret.block<3,3>(0,0).diagonal().setConstant(_M);
  ret.block<3,3>(3,0)=cross<T>(_MC);
  ret.block<3,3>(0,3)=cross<T>(_MC).transpose();
  ret.block<3,3>(3,3)=-_MCCT;
  ret.block<3,3>(3,3).diagonal().array()+=_MCCT.diagonal().sum();
  return ret;
}
Joint::Vec3T Joint::getC() const {
  return _MC/_M;
}
bool Joint::isRotational() const {
  Eigen::Matrix<int,2,1> begEnd=RBegEnd();
  return begEnd[1]-begEnd[0];
}
bool Joint::isRoot(const ArticulatedBody& body) const {
  if(_parent==-1) {
    return true;
  } else {
    //check there are no geometry from this Joint to root Joint
    for(int j=_parent; j>=0; j=body.joint(j)._parent)
      if(body.joint(j)._mesh)
        return false;
    //check Joint type
    Eigen::Matrix<int,2,1> RP=body.joint(_parent).RBegEnd(),R=RBegEnd();
    Eigen::Matrix<int,2,1> CP=body.joint(_parent).CBegEnd(),C=CBegEnd();
    if(RP[1]>RP[0] && C[1]>C[0])
      return true;  //this Joint is translational,parent is rotational
    else if(R[1]>R[0] && CP[1]>CP[0]) {
      return true;  //this Joint is rotational,parent is translational
    }
  }
  return false;
}
std::shared_ptr<ShapeExact> Joint::getGeomPtr() const {
  return _mesh;
}
void Joint::loopAllJointTypes(std::function<void(Joint::JOINT_TYPE)> t) {
  int type=NR_JOINT_TYPE-1;
  while(type>0) {
    t((Joint::JOINT_TYPE) type);
    type>>=1;
  }
}
void Joint::initL1(const ArticulatedBody& body) {
  if(!_mesh)
    return;
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_mesh);
  ASSERT_MSG(mesh,"L1 can only be computed for joints with triangle mesh!")
  ASSERT_MSG(body.nrDOF()!=0,"Body must have at least one movable joint!")
  _L1.setZero(body.nrDOF(),mesh->vss().size());
  for(int i=0; i<(int) mesh->vss().size(); ++i) {
    //L1 of current joint
    bool isTranslational=false;
    double accumulatedL1=sqrt((double) mesh->vss()[i].squaredNorm());
    if(_typeJoint==HINGE_JOINT || _typeJoint==ROT_3D_XYZ) {
      for(int j=0; j<nrDOF(); j++)
        _L1(_offDOF+j,i)=accumulatedL1;
    } else if(_typeJoint==TRANS_1D) {
      isTranslational=true;
      _L1(_offDOF+0,i)=1;
    } else if(_typeJoint==TRANS_2D) {
      isTranslational=true;
      _L1(_offDOF+0,i)=1;
      _L1(_offDOF+1,i)=1;
    } else if(_typeJoint==TRANS_3D) {
      isTranslational=true;
      _L1(_offDOF+0,i)=1;
      _L1(_offDOF+1,i)=1;
      _L1(_offDOF+2,i)=1;
    } else {
      ASSERT_MSG(_typeJoint==FIX_JOINT,"We only support hinge or translation joint!")
    }
    //L1 of parent joints
    int parentId=_parent,childId=0;
    while(parentId!=-1) {
      ASSERT_MSG(body.joint(parentId)._typeJoint==HINGE_JOINT || body.joint(parentId)._typeJoint==FIX_JOINT ||
                 body.joint(parentId)._typeJoint==TRANS_1D || body.joint(parentId)._typeJoint==TRANS_2D ||
                 body.joint(parentId)._typeJoint==TRANS_3D || body.joint(parentId)._typeJoint==ROT_3D_XYZ,
                 "We only support hinge or translation joint!")
      if(body.joint(parentId)._typeJoint==HINGE_JOINT || body.joint(parentId)._typeJoint==ROT_3D_XYZ) {
        if(isTranslational) {
          ASSERT_MSGV(_typeJoint==TRANS_1D,"We only accept Hinge Joint between Trans1D joint and root joint, ParentId=%d ParentName=%s CurrentName=%s",parentId,body.joint(parentId)._name.c_str(),_name.c_str())
          ASSERT_MSGV(body.joint(parentId)._typeJoint!=ROT_3D_XYZ,"We only accept Hinge Joint between Trans1D joint and root joint, not ROT_3D joint, ParentId=%d ParentName=%s CurrentName=%s",parentId,body.joint(parentId)._name.c_str(),_name.c_str())
          ASSERT_MSG(parentId==_parent,"We only support Hinge Joint being the direct parent joint of a TRANS_1D joint!")
          int rotationalParentJointNum=0;
          int tmpParentId=_parent;
          // only allow 0 or 1 rotational joint between a translational joint and root joint
          while(tmpParentId!=-1) {
            // number check for all rotational joint, accept ROT_3D here
            if(body.joint(tmpParentId)._typeJoint==HINGE_JOINT || body.joint(tmpParentId)._typeJoint==ROT_3D_XYZ) {
              rotationalParentJointNum++;
              ASSERT_MSGV(rotationalParentJointNum<=1,
                          "Number of rotational joints between current translational joint and root joint cannot be more than 1. We got %d rotational joints. Current Translational Joint Name=%s",
                          rotationalParentJointNum,
                          _name.c_str())
            }
            tmpParentId=body.joint(tmpParentId)._parent;
          }
          // axis based on the reference frame of rotational joint
          Vec3T transAxis=ROT(body.joint(_parent)._trans)*Vec3T::UnitX();
          Vec3T rotAxis=Vec3T::UnitZ();
          ASSERT_MSGV(abs(abs(transAxis.transpose()*rotAxis)-1)>1e-7,
                      "We do not support translational joint being a child of rotational joint with non-orthogonal axis! ParentID=%d, Axis dot product=%lf",
                      parentId,
                      (double) (transAxis.transpose()*rotAxis))
        }
        if(parentId==_parent) {
          for(int j=0; j<body.joint(parentId).nrDOF(); j++) {
            if(isTranslational)
              _L1(body.joint(parentId)._offDOF+j,i)=accumulatedL1+CTR(_trans).norm();
            else
              _L1(body.joint(parentId)._offDOF+j,i)=accumulatedL1+CTR(_trans).norm();
          }
        } else {
          // should take the trans of childId
          // trans is the translation from the parent's the reference frame to the current one
          for(int j=0; j<body.joint(parentId).nrDOF(); j++) {
            _L1(body.joint(parentId)._offDOF+j,i)=accumulatedL1+CTR(body.joint(childId)._trans).norm();
          }
        }
        accumulatedL1=_L1(body.joint(parentId)._offDOF,i);
      } else if(body.joint(parentId)._typeJoint==TRANS_1D) {
//        isTranslational=true;
        _L1(body.joint(parentId)._offDOF+0,i)=1;
      } else if(body.joint(parentId)._typeJoint==TRANS_2D) {
//        isTranslational=true;
        _L1(body.joint(parentId)._offDOF+0,i)=1;
        _L1(body.joint(parentId)._offDOF+1,i)=1;
      } else if(body.joint(parentId)._typeJoint==TRANS_3D) {
//        isTranslational=true;
        _L1(body.joint(parentId)._offDOF+0,i)=1;
        _L1(body.joint(parentId)._offDOF+1,i)=1;
        _L1(body.joint(parentId)._offDOF+2,i)=1;
      } else {
        ASSERT_MSG(body.joint(parentId)._typeJoint==FIX_JOINT,"We only support hinge or translation joint!")
      }
      childId=parentId;
      parentId=body.joint(parentId)._parent;
    }
  }
}
Joint::VecCM Joint::getL1(int vertexId) const {
  return Eigen::Map<const Joint::Vec>(_L1.data()+_L1.rows()*vertexId,_L1.rows());
}
}
