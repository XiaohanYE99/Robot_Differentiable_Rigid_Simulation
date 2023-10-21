#include "CollisionHandler.h"
#include "DistanceFunction.h"
#include "CollisionBarrierEnergy.h"
#include "ConvexHullDistanceEnergy.h"
#include <Environment/ConvexDecomposition.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Utils/SpatialRotationUtils.h>
#include <Utils/Timing.h>
#include <Utils/Interp.h>
#include <Utils/Utils.h>
#include <unordered_set>
#include <chrono>
#include <stack>

namespace PHYSICSMOTION {
template <typename T>
typename CollisionHandler<T>::Vec CollisionHandler<T>::initializeControlPoints
(T randomInitialize,bool neutralInit,
 std::shared_ptr<ArticulatedBody> body,
 ThetaTrajectory<T>& thetaTrajectory) {
  Vec controlPoints;
  if(!neutralInit) {
    if(randomInitialize!=0) {
      srand((unsigned int)time(0));
      controlPoints=Vec::Random(thetaTrajectory.getNumDOF())*randomInitialize;
    } else controlPoints=Vec::Zero(thetaTrajectory.getNumDOF());
  } else {
    Vec theta=Vec::Zero(body->nrDOF());
    for(int i=0; i<body->nrDOF(); ++i) {
      T upper=body->upperLimit()[i];
      T lower=body->lowerLimit()[i];
      if(isfinite(upper) && isfinite(lower)) {
        theta[i]=(upper+lower)/2;
      } else if(isfinite(upper) && (!isfinite(lower))) {
        if(upper>0)theta[i]=0;
        else theta[i]=upper-0.1;
      } else if(isfinite(lower) && (!isfinite(upper))) {
        if(lower<0) theta[i]=0;
        else theta[i]=lower+0.1;
      }
    }
    MatT neutralCP=theta.rowwise().replicate(thetaTrajectory.getNumDOF()/body->nrDOF());
    controlPoints=Eigen::Map<Vec>(neutralCP.data(),thetaTrajectory.getNumDOF());
  }
  return controlPoints;
}
template <typename T>
CollisionHandler<T>::CollisionHandler(T d0,T x0,T l2,T eta,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD)
  :_d0(d0),_d0Sqr(d0*d0),_l2(l2),_eta(eta),_x0(x0),_CCDObs(true),_CCDSelf(true),_skipJIDPairs(skipSelfCCD) {}
template <typename T>
CollisionHandler<T>::CollisionHandler(std::shared_ptr<ArticulatedBody> body,int order,int totalTime,T d0,T x0,T l2,T eta,T randomInitialize,bool neutralInit,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD)
  :CollisionHandler<T>(d0,x0,l2,eta,skipSelfCCD) {
  _body=body;
  initL1();
  _thetaTrajectory=ThetaTrajectory<T>(_body->nrDOF(),order,totalTime);
  _controlPoints=initializeControlPoints(randomInitialize,neutralInit,body,_thetaTrajectory);
}
template <typename T>
bool CollisionHandler<T>::read(std::istream& is,IOData* dat) {
  registerType<Node<EntityId<T>,BBoxExact>>(dat);
  registerType<ArticulatedBody>(dat);
  registerType<MeshExact>(dat);
  readBinaryData(_body,is,dat);
  readBinaryData(_bvhBody,is,dat);
  readBinaryData(_obstacles,is,dat);
  readBinaryData(_bvhObstacle,is,dat);
  readBinaryData(_d0,is);
  readBinaryData(_d0Sqr,is);
  readBinaryData(_l2,is);
  readBinaryData(_eta,is);
  readBinaryData(_x0,is);
  readBinaryData(_CCDObs,is);
  readBinaryData(_CCDSelf,is);
  readBinaryData(_controlPoints,is);
  readBinaryData(_thetaTrajectory,is);
  readBinaryData(_obsCCPlanes,is,dat);
  readBinaryData(_selfCCPlanes,is,dat);
  readBinaryData(_skipJIDPairs,is);
  assembleObsBVH();
  return is.good();
}
template <typename T>
bool CollisionHandler<T>::write(std::ostream& os,IOData* dat) const {
  registerType<Node<EntityId<T>,BBoxExact>>(dat);
  registerType<ArticulatedBody>(dat);
  registerType<MeshExact>(dat);
  writeBinaryData(_body,os,dat);
  writeBinaryData(_bvhBody,os,dat);
  writeBinaryData(_obstacles,os,dat);
  writeBinaryData(_bvhObstacle,os,dat);
  writeBinaryData(_d0,os);
  writeBinaryData(_d0Sqr,os);
  writeBinaryData(_l2,os);
  writeBinaryData(_eta,os);
  writeBinaryData(_x0,os);
  writeBinaryData(_CCDObs,os);
  writeBinaryData(_CCDSelf,os);
  writeBinaryData(_controlPoints,os);
  writeBinaryData(_thetaTrajectory,os);
  writeBinaryData(_obsCCPlanes,os,dat);
  writeBinaryData(_selfCCPlanes,os,dat);
  writeBinaryData(_skipJIDPairs,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> CollisionHandler<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new CollisionHandler());
}
template <typename T>
std::string CollisionHandler<T>::type() const {
  return typeid(CollisionHandler).name();
}
template <typename T>
void CollisionHandler<T>::setRobot(std::shared_ptr<ArticulatedBody> body,bool randomInitialize,bool neutralInit) {
  _body=body;
  initL1();
  _thetaTrajectory=ThetaTrajectory<T>(_body->nrDOF(),_thetaTrajectory.getOrder(),_thetaTrajectory.getNumSegment());
  _controlPoints=initializeControlPoints(randomInitialize,neutralInit,body,_thetaTrajectory);
}
template <typename T>
void CollisionHandler<T>::addRobot(std::shared_ptr<ArticulatedBody> body,const Vec& DOF) {
  PBDArticulatedGradientInfo<T> info(*body,DOF);
  for(int jid=0; jid<body->nrJ(); jid++) {
    const Joint& J=body->joint(jid);
    if(J._mesh) {
      Mat3X4T t,tJ=TRANSI(info._TM,jid),tM=J._transMesh.template cast<T>();
      APPLY_TRANS(t,tJ,tM)
      addMesh(J._mesh,t);
    }
  }
}
template <typename T>
void CollisionHandler<T>::addMesh(std::shared_ptr<ShapeExact> shape,const Mat3X4T& trans) {
  std::shared_ptr<CompositeShapeExact> shapeC=std::dynamic_pointer_cast<CompositeShapeExact>(shape);
  std::shared_ptr<MeshExact> shapeM=std::dynamic_pointer_cast<MeshExact>(shape);
  if(shapeC) {
    for(int i=0; i<(int)shapeC->getGeoms().size(); i++) {
      Mat3X4T t,tM=shapeC->getTrans()[i].template cast<T>();
      APPLY_TRANS(t,trans,tM)
      addMesh(shapeC->getGeoms()[i],t);
    }
  } else if(shapeM) {
    shapeM=std::dynamic_pointer_cast<MeshExact>(shapeM->copy());
    shapeM->transform(trans.template cast<MeshExact::T>());
    _obstacles.push_back(shapeM);
  }
}
template <typename T>
void CollisionHandler<T>::addObject(const std::string& path,T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHull) {
  std::shared_ptr<ConvexDecomposition> cd;
  if(exists(path+".convex")) {
    cd.reset(new ConvexDecomposition());
    cd->SerializableBase::readStr(path+".convex");
  } else {
    cd.reset(new ConvexDecomposition(path,(double)scale,pos.template cast<double>(),rot.template cast<double>(),maxConvexHull));
    cd->SerializableBase::writeStr(path+".convex");
  }
  for(auto& ob:cd->getConvexHulls())
    _obstacles.push_back(ob);
//  assembleObsBVH();
}
template <typename T>
void CollisionHandler<T>::addSphere(T r,int res,const Vec3T& pos) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addSphere(vss,iss,pos.template cast<double>(),(double)r,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void CollisionHandler<T>::addCuboid(T l,T w,T h,const Vec3T& pos,const Mat3T& R) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addBox(vss,iss,
                        Vec3T(-l/2,-w/2,-h/2).template cast<double>(),
                        Vec3T(l/2,w/2,h/2).template cast<double>());
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  Mat3X4T trans=Mat3X4T::Identity();
  CTR(trans)=pos;
  ROT(trans)=R;
  obstacle->transform(trans.template cast<MeshExact::T>());
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void CollisionHandler<T>::addCapsule(T l,T w,T h,Vec3T pos,T radius,int res) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addCapsule(vss,iss,
                            Vec3T(-l/2,-w/2,-h/2).template cast<double>(),
                            Vec3T(l/2,w/2,h/2).template cast<double>(),(double)radius,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  obstacle->translate(pos.template cast<MeshExact::T>());
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void CollisionHandler<T>::addCapsule(const Vec3T& a,const Vec3T& b,T radius,int res) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addCapsule(vss,iss,a.template cast<double>(),b.template cast<double>(),(double)radius,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
const std::vector<Node<EntityId<T>,BBoxExact>>& CollisionHandler<T>::getBodyBVH() const {
  return _bvhBody;
}
template <typename T>
std::vector<std::shared_ptr<MeshExact>> CollisionHandler<T>::getObstacles() const {
  return _obstacles;
}
template <typename T>
std::shared_ptr<ArticulatedBody> CollisionHandler<T>::getBody() const {
  return _body;
}
template <typename T>
void CollisionHandler<T>::assembleBodyBVH(bool useConvexHull) {
  _bvhBody.clear();
  std::vector<int> ids;
  //copy each body's bvh
  CollisionGradInfo<T> info(*_body,_thetaTrajectory.getPoint(_controlPoints,(0.+(T)_thetaTrajectory.getNumSegment())/2.));
  for(int jid=0; jid<_body->nrJ(); ++jid) {
    if(!_body->joint(jid)._mesh)
      continue;
    std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jid)._mesh);
    const std::vector<Node<int,BBoxExact>>& bvh=mesh->getBVH();
    int currentSize=_bvhBody.size();
    for(int i=0; i<(int)bvh.size(); ++i) {
      const Node<int,BBoxExact>& n=bvh[i];
      if(useConvexHull && n._parent!=-1)
        continue;
      _bvhBody.emplace_back();
      _bvhBody.back()._cell._jid=jid;
      _bvhBody.back()._cell._tid=n._cell;
      _bvhBody.back()._cell._timeFrom=0;
      _bvhBody.back()._cell._timeTo=_thetaTrajectory.getNumSegment();
      _bvhBody.back()._l=(useConvexHull || n._l==-1)? -1: n._l+currentSize;
      _bvhBody.back()._r=(useConvexHull || n._r==-1)? -1: n._r+currentSize;
      _bvhBody.back()._cell._link=_bvhBody.back()._l==-1? mesh: NULL; //only assign link to leaf node, this notify that EntityId of internal nodes == -1
      _bvhBody.back()._parent=n._parent>=0? n._parent+currentSize: -1;
      _bvhBody.back()._nrCell=useConvexHull? 1: n._nrCell;
      if(useConvexHull || n._cell>=0)
        _bvhBody.back()._bb=_bvhBody.back()._cell.computeBB(info);
      else {
        _bvhBody.back()._bb.setUnion(_bvhBody[_bvhBody.back()._l]._bb);
        _bvhBody.back()._bb.setUnion(_bvhBody[_bvhBody.back()._r]._bb);
      }
      if(n._parent==-1)
        ids.push_back((int)_bvhBody.size()-1);
    }
  }
  //link roots
  int szBVH=(int)_bvhBody.size();
  Node<EntityId<T>,BBoxExact>::buildBVHBottomUpAll(_bvhBody);
  for(int i=szBVH; i<(int)_bvhBody.size(); i++) {
    _bvhBody[i]._cell._timeFrom=0;
    _bvhBody[i]._cell._timeTo=_thetaTrajectory.getNumSegment();
  }
  parityCheck();
}
template <typename T>
void CollisionHandler<T>::assembleObsBVH() {
  _bvhObstacle.clear();
  _obstaclePolytopes.resize(_obstacles.size());
  //Just for BVH,all the triangle information is lost.
  for(std::size_t i=0; i<_obstacles.size(); ++i) {
    _bvhObstacle.emplace_back();
    _bvhObstacle.back()._bb=_obstacles[i]->getBB();
    _bvhObstacle.back()._cell=i;
    _bvhObstacle.back()._nrCell=1;
    _obstaclePolytopes[i].reset(new GJKPolytope<T>(_obstacles[i]));
  }
  Node<int,BBoxExact>::buildBVHBottomUpAll(_bvhObstacle);
}
template <typename T>
void CollisionHandler<T>::reset() {
  _infoLookup.clear();
  _obstacles.clear();
  _bvhObstacle.clear();
  assembleBodyBVH();
}
//return success or fail: when success,failedOffset is undefined
//when fail,failedOffset return the leaf node into _bvhBody,which fails the safety-check
template <typename T>
PairStatus CollisionHandler<T>::singlePairCCDSelf(int bvhBodyOffset1,int bvhBodyOffset2,bool doSubdivide,bool buildPair,bool eval) {
  //Narrow phase collision detection and safety check
  int offset1,offset2;
  bool swapped=false;
  if(_bvhBody[bvhBodyOffset1]._cell<_bvhBody[bvhBodyOffset2]._cell) {
    offset1=bvhBodyOffset1;
    offset2=bvhBodyOffset2;
  } else {
    offset1=bvhBodyOffset2;
    offset2=bvhBodyOffset1;
    swapped=true;
  }
  const EntityId<T>& id1=_bvhBody[offset1]._cell;
  const EntityId<T>& id2=_bvhBody[offset2]._cell;
  T dist;
  Vec3T bary,cpa,cpb;
  bool needBuildPair=false;
  const CollisionGradInfo<T>& info=getInfo(id1.getTimeAvg());
  //convexHull-convexHull
  if(id1.isRobotConvexHull() && id2.isRobotConvexHull()) {
    CCDistanceEnergy<T>(id1.getPolytope(info),id2.getPolytope(info)).eval(&dist,(Vec12T*)NULL,NULL);
    if(dist<_d0 && doSubdivide)
      return PS_Penetrated;//Penetrated
    if(dist<(_d0+id1._convexHullPhi+id2._convexHullPhi) && doSubdivide)
      return PS_Unsafe;    //Safety Check Failed
    if(dist-_d0<_x0)
      needBuildPair=true;
  } else {
    if(_thetaTrajectory.getNumSegment()==0) {
      //if user optimizes a static pose, we need discrete collision check for an entire triangle
      Vec3T va[3]= {id1.globalV(info,0),id1.globalV(info,1),id1.globalV(info,2)};
      Vec3T vb[3]= {id2.globalV(info,0),id2.globalV(info,1),id2.globalV(info,2)};
      if(triangleTriangleIntersect(va,vb))
        needBuildPair=true;
    }
    ASSERT_MSG(id1.isRobotTriangle() && id2.isRobotTriangle(),"Narrowphase status error!")
    //vertex-triangle
    for(int d=0; d<3; ++d)
      if(id2.checkBss(d)) {
        VTDistanceEnergy<T>(id1.globalV(info,0),id1.globalV(info,1),
                            id1.globalV(info,2),id2.globalV(info,d)).eval(&dist,(Vec12T*)NULL,NULL);
        if(dist<_d0 && doSubdivide)
          return PS_Penetrated;//Penetrated
        if(dist<(_d0+id1._facePhi+id2._vertexPhi[d]) && doSubdivide)
          return PS_Unsafe;    //Safety Check Failed
        if(dist-_d0<_x0)
          needBuildPair=true;
      }
    //triangle-vertex
    for(int d=0; d<3; ++d)
      if(id1.checkBss(d)) {
        VTDistanceEnergy<T>(id2.globalV(info,0),id2.globalV(info,1),
                            id2.globalV(info,2),id1.globalV(info,d)).eval(&dist,(Vec12T*)NULL,NULL);
        if(dist<_d0 && doSubdivide)
          return PS_Penetrated;//Penetrated
        if(dist<(_d0+id1._vertexPhi[d]+id2._facePhi) && doSubdivide)
          return PS_Unsafe;    //Safety Check Failed
        if(dist-_d0<_x0)
          needBuildPair=true;
      }
    //edge-edge
    for(int d=0; d<3; ++d)
      for(int d2=0; d2<3; ++d2)
        if(id1.checkBss(3+d) && id2.checkBss(3+d2)) {
          EEDistanceEnergy<T>(id1.globalV(info,d),id1.globalV(info,(d+1)%3),
                              id2.globalV(info,d2),id2.globalV(info,(d2+1)%3)).eval(&dist,(Vec12T*)NULL,NULL);
          if(dist<_d0 && doSubdivide)
            return PS_Penetrated;//Penetrated
          if(dist<(_d0+id1._edgePhi[d]+id2._edgePhi[d]) && doSubdivide)
            return PS_Unsafe;    //Safety Check Failed
          if(dist-_d0<_x0)
            needBuildPair=true;
        }
  }
  //add pair
  if(buildPair && needBuildPair) {
    if(eval) {
      if(id1.getTimeAvg()!=id2.getTimeAvg()) {
        if(id1.getTimeDiff()>id2.getTimeDiff()) {
          if(swapped)
            return PS_RightLargerInterval;
          else return PS_LeftLargerInterval;
        } else {
          if(swapped)
            return PS_LeftLargerInterval;
          else return PS_RightLargerInterval;
        }
      }
    }
    if(id1.isRobotConvexHull() && id2.isRobotConvexHull())
      _selfCCPlanes.insertPlane(std::make_pair(id1,id2));
    else _selfTTPairs.push_back(std::make_pair(id1,id2));
  }
  return PS_Safe;
}
template <typename T>
PairStatus CollisionHandler<T>::singlePairCCDObs(int bvhBodyOffset,int obstacleId,int bvhObstacleOffset,bool subdivide,bool buildPair) {
  //Narrow phase collision detection and safety check
  const EntityId<T>& id1=_bvhBody[bvhBodyOffset]._cell;
  EntityId<T> id2(_obstaclePolytopes[obstacleId],bvhObstacleOffset==-1? -1: _obstacles[obstacleId]->getBVH()[bvhObstacleOffset]._cell);
  T dist;
  Vec12T Grad;
  Mat12T Hessian;
  Vec3T bary,cpa,cpb;
  bool needBuildPair=false;
  const CollisionGradInfo<T>& info=getInfo(id1.getTimeAvg());
  if(id1.isRobotConvexHull() && id2.isObstacleConvexHull()) {
    //convexHull-convexHull
    CCDistanceEnergy<T>(id1.getPolytope(info),*_obstaclePolytopes[obstacleId]).eval(&dist,(Vec12T*)NULL,NULL);
    if(dist<_d0 && subdivide) {
      return PS_Penetrated;//Penetrated
    }
    if(dist<(_d0+id1._convexHullPhi) && subdivide)
      return PS_Unsafe;    //Safety Check Failed
    if(dist-_d0<_x0)
      needBuildPair=true;
  } else {
    if(_thetaTrajectory.getNumSegment()==0) {
      //if user optimizes a static pose, we need discrete collision check for an entire triangle
      Vec3T va[3]= {id1.globalV(info,0),id1.globalV(info,1),id1.globalV(info,2)};
      Vec3T vb[3]= {id2.globalV(info,0),id2.globalV(info,1),id2.globalV(info,2)};
      if(triangleTriangleIntersect(va,vb))
        needBuildPair=true;
    }
    ASSERT_MSG(id1.isRobotTriangle() && id2.isObstacleTriangle(),"Narrowphase status error!")
    //vertex-triangleObs
    for(int d=0; d<3; ++d)
      if(id2.checkBss(d)) {
        VTDistanceEnergy<T>(id1.globalV(info,0),id1.globalV(info,1),
                            id1.globalV(info,2),id2.globalV(info,d)).eval(&dist,(Vec12T*)NULL,NULL);
        if(dist<_d0 && subdivide)
          return PS_Penetrated;//Penetrated
        if(dist<(_d0+id1._vertexPhi[d]) && subdivide)
          return PS_Unsafe;    //Safety Check Failed
        if(dist-_d0<_x0)
          needBuildPair=true;
      }
    //triangleObs-vertex
    for(int d=0; d<3; ++d)
      if(id1.checkBss(d)) {
        VTDistanceEnergy<T>(id2.globalV(info,0),id2.globalV(info,1),
                            id2.globalV(info,2),id1.globalV(info,d)).eval(&dist,(Vec12T*)NULL,NULL);
        if(dist<_d0 && subdivide)
          return PS_Penetrated;//Penetrated
        if(dist<(_d0+id1._facePhi) && subdivide)
          return PS_Unsafe;    //Safety Check Failed
        if(dist-_d0<_x0)
          needBuildPair=true;
      }
    //edge-edge
    for(int d=0; d<3; ++d)
      for(int d2=0; d2<3; ++d2)
        if(id1.checkBss(3+d) && id2.checkBss(3+d2)) {
          EEDistanceEnergy<T>(id1.globalV(info,d),id1.globalV(info,(d+1)%3),
                              id2.globalV(info,d2),id2.globalV(info,(d2+1)%3)).eval(&dist,(Vec12T*)NULL,NULL);
          if(dist<_d0 && subdivide)
            return PS_Penetrated;//Penetrated
          if(dist<(_d0+id1._edgePhi[d]) && subdivide)
            return PS_Unsafe;    //Safety Check Failed
          if(dist-_d0<_x0)
            needBuildPair=true;
        }
  }
  //add pair
  if(buildPair && needBuildPair) {
    if(id1.isRobotConvexHull() && id2.isObstacleConvexHull())
      _obsCCPlanes.insertPlane(std::make_pair(id1,id2));
    else _obsTTPairs.push_back(std::make_pair(id1,id2));
  }
  return PS_Safe;
}
template <typename T>
void CollisionHandler<T>::narrowphaseCCD(int bvhBodyOffset,int obstacleId,bool subdivide,bool buildPair) {
#define IS_LEAF_BODY(offset)(_bvhBody[offset]._l==-1)
#define IS_LEAF_OBS(offset) (obs->getBVH()[offset]._l==-1)
  //we are using convex hull
  if(_bvhBody[bvhBodyOffset]._cell._tid==-1) {
    if(subdivide && _subdivideOffsets.find(bvhBodyOffset)!=_subdivideOffsets.end())
      return;   //this is already labeled for subdivision
    PairStatus s=singlePairCCDObs(bvhBodyOffset,obstacleId,-1,subdivide,buildPair);
    if(s!=PS_Safe)
      _subdivideOffsets.insert(bvhBodyOffset);
    if(s==PS_Penetrated)
      _penetratedOffsets.insert(bvhBodyOffset);
    return;
  }
  //we are not using convex hull
  std::stack<std::pair<int,int>> ss;
  std::shared_ptr<MeshExact> obs=_obstacles[obstacleId];
  ss.push(std::make_pair(bvhBodyOffset,(int)obs->getBVH().size()-1));
  while(!ss.empty()) {
    bvhBodyOffset=ss.top().first;
    int bvhObstacleOffset=ss.top().second;
    ss.pop();
    if(!_bvhBody[bvhBodyOffset]._bb.intersect(obs->getBVH()[bvhObstacleOffset]._bb))
      continue;
    if(IS_LEAF_BODY(bvhBodyOffset) && IS_LEAF_OBS(bvhObstacleOffset)) {
      if(subdivide && _subdivideOffsets.find(bvhBodyOffset)!=_subdivideOffsets.end())
        continue;   //this is already labeled for subdivision
      PairStatus s=singlePairCCDObs(bvhBodyOffset,obstacleId,bvhObstacleOffset,subdivide,buildPair);
      if(s!=PS_Safe)
        _subdivideOffsets.insert(bvhBodyOffset);
      if(s==PS_Penetrated)
        _penetratedOffsets.insert(bvhBodyOffset);
    } else if(IS_LEAF_BODY(bvhBodyOffset)) {
      ss.push(std::make_pair(bvhBodyOffset,obs->getBVH()[bvhObstacleOffset]._l));
      ss.push(std::make_pair(bvhBodyOffset,obs->getBVH()[bvhObstacleOffset]._r));
    } else if(IS_LEAF_OBS(bvhObstacleOffset)) {
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,bvhObstacleOffset));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,bvhObstacleOffset));
    } else {
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,obs->getBVH()[bvhObstacleOffset]._l));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,obs->getBVH()[bvhObstacleOffset]._r));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,obs->getBVH()[bvhObstacleOffset]._l));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,obs->getBVH()[bvhObstacleOffset]._r));
    }
  }
#undef IS_LEAF_BODY
#undef IS_LEAF_OBS
}
template <typename T>
bool CollisionHandler<T>::CCD(bool doSubdivide,bool buildPair,bool eval) {
  bool safe=true,clear=true;
  if(_CCDObs) {
    safe=CCDObs(doSubdivide,buildPair,clear) && safe;
    clear=false;
  }
  if(_CCDSelf) {
    safe=CCDSelf(doSubdivide,buildPair,clear,eval) && safe;
    clear=false;
  }
  return safe;
}
template <typename T>
bool CollisionHandler<T>::CCDObs(bool subdivide,bool buildPair,bool clear) {
#define IS_LEAF_BODY(offset)(_bvhBody[offset]._l==-1)
#define IS_LEAF_OBS(offset) (_bvhObstacle[offset]._l==-1)
  ASSERT_MSG(!_bvhBody.empty(),"Called CCDObs, but _bvhBody is empty!")
  TBEG("CCDObs");
  _obsTTPairs.clear();
  if(clear) {
    _penetratedOffsets.clear();
    _subdivideOffsets.clear();
  }
  std::stack<std::pair<int,int>> ss;
  if(!_bvhObstacle.empty()) {
    ss.push(std::make_pair((int)_bvhBody.size()-1,(int)_bvhObstacle.size()-1));
    while(!ss.empty()) {
      int bvhBodyOffset=ss.top().first;
      int bvhObstacleOffset=ss.top().second;
      ss.pop();
      if(!_bvhBody[bvhBodyOffset]._bb.intersect(_bvhObstacle[bvhObstacleOffset]._bb))
        continue;
      else if(IS_LEAF_BODY(bvhBodyOffset) && IS_LEAF_OBS(bvhObstacleOffset))
        narrowphaseCCD(bvhBodyOffset,
                       _bvhObstacle[bvhObstacleOffset]._cell,
                       subdivide,buildPair);
      else if(IS_LEAF_BODY(bvhBodyOffset)) {
        ss.push(std::make_pair(bvhBodyOffset,_bvhObstacle[bvhObstacleOffset]._l));
        ss.push(std::make_pair(bvhBodyOffset,_bvhObstacle[bvhObstacleOffset]._r));
      } else if(IS_LEAF_OBS(bvhObstacleOffset)) {
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,bvhObstacleOffset));
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,bvhObstacleOffset));
      } else {
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,_bvhObstacle[bvhObstacleOffset]._l));
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._l,_bvhObstacle[bvhObstacleOffset]._r));
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,_bvhObstacle[bvhObstacleOffset]._l));
        ss.push(std::make_pair(_bvhBody[bvhBodyOffset]._r,_bvhObstacle[bvhObstacleOffset]._r));
      }
    }
  }
  TEND();
  return _subdivideOffsets.empty() && _penetratedOffsets.empty();
#undef IS_LEAF_BODY
#undef IS_LEAF_OBS
}
template <typename T>
bool CollisionHandler<T>::CCDSelf(bool doSubdivision,bool buildPair,bool clear,bool eval) {
#define IS_LEAF_BODY(offset)(_bvhBody[offset]._l==-1)
  ASSERT_MSG(!_bvhBody.empty(),"Called CCDSelf, but _bvhBody is empty!")
  TBEG("CCDSelf");
  _selfTTPairs.clear();
  if(clear) {
    _penetratedOffsets.clear();
    _subdivideOffsets.clear();
  }
  if(eval)
    _evaluateSubdivisionOffsets.clear();
  std::stack<std::pair<int,int>> ss;
  ss.push(std::make_pair((int)_bvhBody.size()-1,(int)_bvhBody.size()-1));
  while(!ss.empty()) {
    int bvhBodyOffset1=ss.top().first;
    int bvhBodyOffset2=ss.top().second;
    ss.pop();
    if(!_bvhBody[bvhBodyOffset1]._bb.intersect(_bvhBody[bvhBodyOffset2]._bb))
      continue;
    if(_thetaTrajectory.getNumSegment()>0)
      if(_bvhBody[bvhBodyOffset1]._cell._timeTo<=_bvhBody[bvhBodyOffset2]._cell._timeFrom ||
          _bvhBody[bvhBodyOffset2]._cell._timeTo<=_bvhBody[bvhBodyOffset1]._cell._timeFrom)
        continue;
    if(excludedFromCCDSelf(_bvhBody[bvhBodyOffset1]._cell._jid,_bvhBody[bvhBodyOffset2]._cell._jid))
      continue;
    if(IS_LEAF_BODY(bvhBodyOffset1) && IS_LEAF_BODY(bvhBodyOffset2)) {
      if(bvhBodyOffset1>=bvhBodyOffset2)
        continue;
      if(doSubdivision && (_subdivideOffsets.find(bvhBodyOffset1)!=_subdivideOffsets.end() || _subdivideOffsets.find(bvhBodyOffset2)!=_subdivideOffsets.end()))
        continue;   //this is already labeled for subdivision
      PairStatus s=singlePairCCDSelf(bvhBodyOffset1,bvhBodyOffset2,doSubdivision,buildPair,eval);
      if(s!=PS_Safe) {
        if(eval) {
          if(s==PS_LeftLargerInterval)
            _evaluateSubdivisionOffsets.insert(bvhBodyOffset1);
          else if(s==PS_RightLargerInterval)
            _evaluateSubdivisionOffsets.insert(bvhBodyOffset2);
        } else {
          _subdivideOffsets.insert(bvhBodyOffset1);
          _subdivideOffsets.insert(bvhBodyOffset2);
        }
      }
      if(s==PS_Penetrated) {
        std::cout<<"Penetrated in selfCCD!"<<std::endl;
        _bvhBody[bvhBodyOffset1]._cell.print();
        _bvhBody[bvhBodyOffset2]._cell.print();
        _penetratedOffsets.insert(bvhBodyOffset1);
        _penetratedOffsets.insert(bvhBodyOffset2);
      }
    } else if(IS_LEAF_BODY(bvhBodyOffset1)) {
      ss.push(std::make_pair(bvhBodyOffset1,_bvhBody[bvhBodyOffset2]._l));
      ss.push(std::make_pair(bvhBodyOffset1,_bvhBody[bvhBodyOffset2]._r));
    } else if(IS_LEAF_BODY(bvhBodyOffset2)) {
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._l,bvhBodyOffset2));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._r,bvhBodyOffset2));
    } else {
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._l,_bvhBody[bvhBodyOffset2]._l));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._l,_bvhBody[bvhBodyOffset2]._r));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._r,_bvhBody[bvhBodyOffset2]._l));
      ss.push(std::make_pair(_bvhBody[bvhBodyOffset1]._r,_bvhBody[bvhBodyOffset2]._r));
    }
  }
  TEND();
  if(eval && !_evaluateSubdivisionOffsets.empty()) {
    subdivide(ST_Eval);
    CCDSelf(doSubdivision,buildPair,true,eval);
  }
  return _subdivideOffsets.empty() && _penetratedOffsets.empty();
#undef IS_LEAF_BODY
}
template <typename T>
bool CollisionHandler<T>::CCDBF(bool doSubdivide,bool buildPair,bool eval) {
  bool safe=true,clear=true;
  if(_CCDObs) {
    safe=CCDObsBF(doSubdivide,buildPair,clear) && safe;
    clear=false;
  }
  if(_CCDSelf) {
    safe=CCDSelfBF(doSubdivide,buildPair,clear,eval) && safe;
    clear=false;
  }
  return safe;
}
template <typename T>
bool CollisionHandler<T>::CCDObsBF(bool subdivide,bool buildPair,bool clear) {
#define IS_LEAF_BODY(offset)(_bvhBody[offset]._l==-1)
#define IS_LEAF_OBS(offset) (obs->getBVH()[offset]._l==-1)
  ASSERT_MSG(!_bvhBody.empty(),"Called CCDObsBF, but _bvhBody is empty!")
  _obsTTPairs.clear();
  if(clear) {
    _penetratedOffsets.clear();
    _subdivideOffsets.clear();
  }
  if(!_bvhObstacle.empty()) {
    for(int bvhBodyOffset=0; bvhBodyOffset<(int)_bvhBody.size(); ++bvhBodyOffset) {
      if(!IS_LEAF_BODY(bvhBodyOffset))
        continue;
      for(int obstacleId=0; obstacleId<(int)_obstacles.size(); ++obstacleId) {
        //we are using convex hull
        if(_bvhBody[bvhBodyOffset]._cell.isRobotConvexHull()) {
          if(subdivide && _subdivideOffsets.find(bvhBodyOffset)!=_subdivideOffsets.end())
            continue;   //this is already labeled for subdivision
          PairStatus s=singlePairCCDObs(bvhBodyOffset,obstacleId,-1,subdivide,buildPair);
          if(s!=PS_Safe)
            _subdivideOffsets.insert(bvhBodyOffset);
          if(s==PS_Penetrated)
            _penetratedOffsets.insert(bvhBodyOffset);
          continue;
        }
        //we are not using convex hull
        std::shared_ptr<MeshExact> obs=std::dynamic_pointer_cast<MeshExact>(_obstacles[obstacleId]);
        for(int bvhObstacleOffset=0; bvhObstacleOffset<(int)obs->getBVH().size(); ++bvhObstacleOffset) {
          if(!IS_LEAF_OBS(bvhObstacleOffset))
            continue;
          if(_bvhBody[bvhBodyOffset]._bb.intersect(obs->getBVH()[bvhObstacleOffset]._bb)) {
            if(subdivide && _subdivideOffsets.find(bvhBodyOffset)!=_subdivideOffsets.end())
              continue;   //this is already labeled for subdivision
            PairStatus s=singlePairCCDObs(bvhBodyOffset,obstacleId,bvhObstacleOffset,subdivide,buildPair);
            if(s!=PS_Safe)
              _subdivideOffsets.insert(bvhBodyOffset);
            if(s==PS_Penetrated)
              _penetratedOffsets.insert(bvhBodyOffset);
          }
        }
      }
    }
  }
  return (_subdivideOffsets.empty() && _penetratedOffsets.empty());
#undef IS_LEAF_BODY
#undef IS_LEAF_OBS
}
template <typename T>
bool CollisionHandler<T>::CCDSelfBF(bool doSubdivide,bool buildPair,bool clear,bool eval) {
#define IS_LEAF_BODY(offset)(_bvhBody[offset]._l==-1)
  ASSERT_MSG(!_bvhBody.empty(),"Called CCDSelfBF, but _bvhBody is empty!")
  TBEG("selfCCDBF");
  _selfTTPairs.clear();
  if(clear) {
    _penetratedOffsets.clear();
    _subdivideOffsets.clear();
  }
  if(eval)
    _evaluateSubdivisionOffsets.clear();
  for(int bvhBodyOffset1=0; bvhBodyOffset1<(int)_bvhBody.size(); ++bvhBodyOffset1) {
    if(!IS_LEAF_BODY(bvhBodyOffset1))
      continue;
    for(int bvhBodyOffset2=bvhBodyOffset1+1; bvhBodyOffset2<(int)_bvhBody.size(); ++bvhBodyOffset2) {
      if(!IS_LEAF_BODY(bvhBodyOffset2))
        continue;
      if(excludedFromCCDSelf(_bvhBody[bvhBodyOffset1]._cell._jid,_bvhBody[bvhBodyOffset2]._cell._jid))
        continue;
      if(!_bvhBody[bvhBodyOffset1]._bb.intersect(_bvhBody[bvhBodyOffset2]._bb))
        continue;
      if(_bvhBody[bvhBodyOffset1]._cell._timeTo<=_bvhBody[bvhBodyOffset2]._cell._timeFrom ||
          _bvhBody[bvhBodyOffset2]._cell._timeTo<=_bvhBody[bvhBodyOffset1]._cell._timeFrom)
        continue;
      if(doSubdivide && (_subdivideOffsets.find(bvhBodyOffset1)!=_subdivideOffsets.end() || _subdivideOffsets.find(bvhBodyOffset2)!=_subdivideOffsets.end()))
        continue;   //this is already labeled for subdivision
      PairStatus s=singlePairCCDSelf(bvhBodyOffset1,bvhBodyOffset2,doSubdivide,buildPair,eval);
      if(s!=PS_Safe) {
        if(eval) {
          if(s==PS_LeftLargerInterval)
            _evaluateSubdivisionOffsets.insert(bvhBodyOffset1);
          else if(s==PS_RightLargerInterval)
            _evaluateSubdivisionOffsets.insert(bvhBodyOffset2);
        } else {
          _subdivideOffsets.insert(bvhBodyOffset1);
          _subdivideOffsets.insert(bvhBodyOffset2);
        }
      }
      if(s==PS_Penetrated) {
        _penetratedOffsets.insert(bvhBodyOffset1);
        _penetratedOffsets.insert(bvhBodyOffset2);
      }
    }
  }
  TEND();
  if(eval && !_evaluateSubdivisionOffsets.empty()) {
    subdivide(ST_Eval);
    CCDSelfBF(doSubdivide,buildPair,true,eval);
  }
  return _subdivideOffsets.empty() && _penetratedOffsets.empty();
#undef IS_LEAF_BODY
}
template <typename T>
void CollisionHandler<T>::debugCCDObsBF(bool doSubdivide,bool buildPair,bool exhaustive) {
  std::vector<std::pair<EntityId<T>,EntityId<T>>> TTPairs,TTPairsBF;
  CCSeparatingPlanes<T> CCPlanes,CCPlanesBF;
  std::unordered_set<int> subdivideOffsets,subdivideOffsetsBF;
  std::unordered_set<int> penetratedOffsets,penetratedOffsetsBF;
  if(exhaustive) {
    exhaustiveSubdivide(false,false);
    TTPairs=_obsTTPairs;
    exhaustiveSubdivide(true,false);
    TTPairsBF=_obsTTPairs;
  } else {
    std::vector<Node<EntityId<T>,BBoxExact>> tmpBvhBody=_bvhBody;
    update();
    _obsCCPlanes.clear();   //normally, CCPlanes are persistent but we clear it for debug
    CCDObs(doSubdivide,buildPair,true);
    subdivideOffsets=_subdivideOffsets;
    penetratedOffsets=_penetratedOffsets;
    TTPairs=_obsTTPairs;
    CCPlanes=_obsCCPlanes;

    _bvhBody=tmpBvhBody;
    update();
    _obsCCPlanes.clear();   //normally, CCPlanes are persistent but we clear it for debug
    CCDObsBF(doSubdivide,buildPair,true);
    subdivideOffsetsBF=_subdivideOffsets;
    penetratedOffsetsBF=_penetratedOffsets;
    TTPairsBF=_obsTTPairs;
    CCPlanesBF=_obsCCPlanes;
  }
  //compare TTPairs
  compareTTPair(TTPairs,TTPairsBF);
  //compute CCPairs
  CCSeparatingPlanes<T>::compareCCPair(CCPlanes,CCPlanesBF);
  //compare subdivideOffsets
  std::cout<<"Found "<<subdivideOffsets.size()<<" subdivideOffsets "<<subdivideOffsetsBF.size()<<" subdivideOffsetsBF!"<<std::endl;
  ASSERT_MSG(subdivideOffsets==subdivideOffsetsBF,"subdivideOffsets not same!")
  //compare penetratedOffsets
  std::cout<<"Found "<<penetratedOffsets.size()<<" penetratedOffsets "<<penetratedOffsetsBF.size()<<" penetratedOffsetsBF!"<<std::endl;
  ASSERT_MSG(penetratedOffsets==penetratedOffsetsBF,"penetratedOffsets not same!")
  std::cout<<"debugCCDObsBF success!"<<std::endl;
}
template <typename T>
void CollisionHandler<T>::debugCCDSelfBF(bool doSubdivide,bool buildPair,bool exhaustive,bool eval) {
  std::vector<std::pair<EntityId<T>,EntityId<T>>> TTPairs,TTPairsBF;
  CCSeparatingPlanes<T> CCPlanes,CCPlanesBF;
  std::unordered_set<int> subdivideOffsets,subdivideOffsetsBF;
  std::unordered_set<int> penetratedOffsets,penetratedOffsetsBF;
  std::unordered_set<int> evaluateSubdivisionOffsets,evaluateSubdivisionOffsetsBF;
  if(exhaustive) {
    exhaustiveSubdivide(false,eval);
    TTPairs=_selfTTPairs;
    exhaustiveSubdivide(true,eval);
    TTPairsBF=_selfTTPairs;
  } else {
    std::vector<Node<EntityId<T>,BBoxExact>> tmpBvhBody=_bvhBody;
    update();
    _obsCCPlanes.clear();   //normally, CCPlanes are persistent but we clear it for debug
    CCDSelf(doSubdivide,buildPair,true,eval);
    subdivideOffsets=_subdivideOffsets;
    penetratedOffsets=_penetratedOffsets;
    evaluateSubdivisionOffsets=_evaluateSubdivisionOffsets;
    TTPairs=_selfTTPairs;
    CCPlanes=_selfCCPlanes;

    _bvhBody=tmpBvhBody;
    update();
    _obsCCPlanes.clear();   //normally, CCPlanes are persistent but we clear it for debug
    CCDSelfBF(doSubdivide,buildPair,true,eval);
    subdivideOffsetsBF=_subdivideOffsets;
    penetratedOffsetsBF=_penetratedOffsets;
    evaluateSubdivisionOffsetsBF=_evaluateSubdivisionOffsets;
    TTPairsBF=_selfTTPairs;
    CCPlanesBF=_selfCCPlanes;
  }
  //compare TTPairs
  compareTTPair(TTPairs,TTPairsBF);
  //compare CCPairs
  CCSeparatingPlanes<T>::compareCCPair(CCPlanes,CCPlanesBF);
  //compare subdivideOffsets
  std::cout<<"Found "<<subdivideOffsets.size()<<" subdivideOffsets "<<subdivideOffsetsBF.size()<<" subdivideOffsetsBF!"<<std::endl;
  ASSERT_MSG(subdivideOffsets==subdivideOffsetsBF,"subdivideOffsets not same!")
  //compare penetratedOffsets
  std::cout<<"Found "<<penetratedOffsets.size()<<" penetratedOffsets "<<penetratedOffsetsBF.size()<<" penetratedOffsetsBF!"<<std::endl;
  ASSERT_MSG(penetratedOffsets==penetratedOffsetsBF,"penetratedOffsets not same!")
  //compare evaluateSubdivisionOffsets
  std::cout<<"Found "<<evaluateSubdivisionOffsets.size()<<" evaluateSubdivisionOffsets "<<evaluateSubdivisionOffsetsBF.size()<<" evaluateSubdivisionOffsetsBF!"<<std::endl;
  ASSERT_MSG(evaluateSubdivisionOffsets==evaluateSubdivisionOffsetsBF,"evaluateSubdivisionOffsets not same!")
  std::cout<<"debugCCDSelfBF success!"<<std::endl;
}
//we will input the offset into the vector _bvhBody,and we guarantee _bvhBody[offset] is a leaf node
//we modify the _bvhBody vector,subdividing the leaf into two
template <typename T>
bool CollisionHandler<T>::exhaustiveSubdivide(bool bruteForce,bool eval,int iter) {
  update();
  std::cout<<"Updated bounding boxes!"<<std::endl;
  Node<EntityId<T>,BBoxExact>::parityCheck(_bvhBody);
  int nrSubdivide=0;
  int maxIters=iter==-1? INT_MAX: iter;
  for(int i=0; i<maxIters; ++i) {
    bool ccdPass;
    if(bruteForce)
      ccdPass=CCDBF(true,true,eval);
    else ccdPass=CCD(true,true,eval);
    if(!_penetratedOffsets.empty()) {
      std::cout<<"Penetration Found"<<std::endl;
      return false;
    }
    if(ccdPass) {
      std::cout<<"No subdivision needed!"<<std::endl;
      break;
    } else {
      nrSubdivide++;
      std::cout<<"Subdividing the "<<nrSubdivide<<"th time!"<<std::endl;
      subdivide(ST_All);
    }
  }
  std::cout<<"Finished exhaustiveSubdivide!"<<std::endl;
  return true;
}
template <typename T>
void CollisionHandler<T>::subdivideSingle(int offset,bool assertNoReserve) {
  //remove CCPlane related to offset
  _obsCCPlanes.subdivide(_bvhBody[offset]._cell);
  _selfCCPlanes.subdivide(_bvhBody[offset]._cell);
  //prepare two children
  ASSERT_MSG(_bvhBody[offset]._cell!=-1,"Subdividing internal nodes!")
  EntityId<T> valLeft(_bvhBody[offset]._cell);
  valLeft._timeTo=valLeft.getTimeAvg();
  EntityId<T> valRight(_bvhBody[offset]._cell);
  valRight._timeFrom=valRight.getTimeAvg();
  //subdivide
  bool isLeft=true;
  Node<EntityId<T>,BBoxExact>::insertLeaf(_bvhBody,valLeft,offset,[&](Node<EntityId<T>,BBoxExact>& n) {
    if(isLeft) {
      n._cell=valLeft;
      const GradInfo& info=getInfo(valLeft.getTimeAvg());
      n._bb=valLeft.computeBB(info);
      isLeft=false;
    } else {
      n._cell=valRight;
      const GradInfo& info=getInfo(valRight.getTimeAvg());
      n._bb=valRight.computeBB(info);
    }
  },assertNoReserve? -1: 256);
  //mark _bvhBody[offset]._cell==-1
  _bvhBody[offset]._cell._link=NULL;
}
template <typename T>
void CollisionHandler<T>::subdivide(SubdivideType type,bool debug) {
  TBEG("subdivide");
  //make sure that _subdivideOffsets independently make sense
  CCSeparatingPlanes<T> tmpObsCCPlanes=_obsCCPlanes;
  CCSeparatingPlanes<T> tmpSelfCCPlanes=_selfCCPlanes;
  int nrEmpty=Node<EntityId<T>,BBoxExact>::nrEmpty(_bvhBody),nReserve=0;
  int sizeNeeded;
  switch(type) {
  case ST_All:
    sizeNeeded=(int)(_subdivideOffsets.size()+_evaluateSubdivisionOffsets.size())*2;
    break;
  case ST_Safety:
    sizeNeeded=(int)_subdivideOffsets.size()*2;
    break;
  case ST_Eval:
    sizeNeeded=(int)_evaluateSubdivisionOffsets.size()*2;
    break;
  default:
    ASSERT_MSG(false,"Wrong subdivision type!")
  }
  if(debug)
    std::cout<<"Subdivide: #node="<<_bvhBody.size()<<" depth="<<Node<EntityId<T>,BBoxExact>::depth(_bvhBody)<<" #subd="<<sizeNeeded;
  if(nrEmpty<sizeNeeded) {
    nReserve=sizeNeeded-nrEmpty;
    Node<EntityId<T>,BBoxExact>::reserveEmpty(_bvhBody,nReserve);
  }
  //subdivide all
  if(type==ST_All || type==ST_Safety)
    for(auto offset:_subdivideOffsets)
      subdivideSingle(offset+nReserve,true);
  if(type==ST_All || type==ST_Eval)
    for(auto offset:_evaluateSubdivisionOffsets)
      subdivideSingle(offset+nReserve,true);
  //check and debug output
  parityCheck();
  if(debug) {
    std::cout<<" -> #node="<<_bvhBody.size()<<" depth="<<Node<EntityId<T>,BBoxExact>::depth(_bvhBody)<<"!"<<std::endl;
    _obsCCPlanes.parityCheck();
    _selfCCPlanes.parityCheck();
  }
  update();
  TEND();
}
template <typename T>
std::string CollisionHandler<T>::findEntityId(const EntityId<T>& id) const {
  if(id.isRobotConvexHull() || id.isRobotTriangle())
    for(int i=0; i<(int)_bvhBody.size(); i++)
      if(_bvhBody[i]._cell==id)
        return "B["+std::to_string(i)+","+std::to_string(id._jid)+","+std::to_string(id._tid)+","+std::to_string((double)id.getTimeAvg())+"]";
  if(id.isObstacleConvexHull() || id.isObstacleTriangle())
    for(int i=0; i<(int)_obstacles.size(); i++)
      if(_obstacles[i]==id._obs->mesh())
        return "O["+std::to_string(i)+","+std::to_string(id._jid)+","+std::to_string(id._tid)+"]";
  return "U";
}
template <typename T>
void CollisionHandler<T>::parityCheck() const {
  for(int i=0; i<(int)_bvhBody.size(); i++)
    if(_bvhBody[i]._cell._jid>=0) {
      if(_bvhBody[i]._l>=0 && _bvhBody[i]._r>=0) {
        ASSERT_MSGV(_bvhBody[i]._cell==-1,"internal._cell!=-1 at %dth node!",i)
      } else {
        ASSERT_MSGV(_bvhBody[i]._cell!=-1,"leaf._cell==-1 at %dth node!",i)
      }
      if(_bvhBody[i]._l>=0) {
        ASSERT_MSGV(_bvhBody[_bvhBody[i]._l]._cell._jid==_bvhBody[i]._cell._jid,"cell._l._jid mismatch at %dth node!",i)
      }
      if(_bvhBody[i]._r>=0) {
        ASSERT_MSGV(_bvhBody[_bvhBody[i]._r]._cell._jid==_bvhBody[i]._cell._jid,"cell._r._jid mismatch at %dth node!",i)
      }
    }
}
template <typename T>
void CollisionHandler<T>::clearCC() {
  _obsCCPlanes=CCSeparatingPlanes<T>();
  _selfCCPlanes=CCSeparatingPlanes<T>();
}
template <typename T>
void CollisionHandler<T>::update() {
  //pass-1 gradient info computation
  _infoLookup.clear();
  std::vector<T> timeAvgs;
  for(int i=0; i<(int)_bvhBody.size(); i++) {
    EntityId<T>& id=_bvhBody[i]._cell;
    if(_bvhBody[i]._l!=-1)
      continue;
    timeAvgs.push_back(id.getTimeAvg());
    _infoLookup.insert(std::make_pair(id.getTimeAvg(),GradInfo()));
  }
  //make unique
  std::sort(timeAvgs.begin(),timeAvgs.end());
  auto last=std::unique(timeAvgs.begin(),timeAvgs.end());
  timeAvgs.erase(last,timeAvgs.end());
  //set gradient information
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)timeAvgs.size(); i++)
    _infoLookup.find(timeAvgs[i])->second.reset(*_body,_thetaTrajectory.getPoint(_controlPoints,timeAvgs[i]));
  //pass-2 leaf node variable computation
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_bvhBody.size(); i++) {
    EntityId<T>& id=_bvhBody[i]._cell;
    if(_bvhBody[i]._l!=-1)
      continue;
    //Phi
    if(id.isRobotConvexHull())
      id._convexHullPhi=getConvexHullPhi(i);
    if(id.isRobotTriangle()) {
      id._facePhi=getFacePhi(i);
      for(int d=0; d<3; d++)
        id._edgePhi[d]=getEdgePhi(i,d);
      for(int d=0; d<3; d++)
        id._vertexPhi[d]=getVertexPhi(i,d);
    }
    //BBox
    _bvhBody[i]._bb=computeBB(id,id._timeFrom,id._timeTo);
  }
  //pass-3 merge bounding boxes
  for(int i=0; i<(int)_bvhBody.size(); i++) {
    if(_bvhBody[i]._l==-1)
      continue;
    _bvhBody[i]._bb=BBoxExact();
    _bvhBody[i]._bb.setUnion(_bvhBody[_bvhBody[i]._l]._bb);
    _bvhBody[i]._bb.setUnion(_bvhBody[_bvhBody[i]._r]._bb);
  }
}
//info query
template <typename T>
const std::unordered_map<int,std::unordered_set<int>>& CollisionHandler<T>::getSkipJIDPairs() const {
  return _skipJIDPairs;
}
template <typename T>
const std::unordered_set<int>& CollisionHandler<T>::getSubdivideOffsets() const {
  return _subdivideOffsets;
}
template <typename T>
const std::vector<Node<EntityId<T>,BBoxExact>>& CollisionHandler<T>::getBVHBody() const {
  return _bvhBody;
}
template <typename T>
const std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionHandler<T>::getObsTTPairs() const {
  return _obsTTPairs;
}
template <typename T>
std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionHandler<T>::getObsTTPairs() {
  return _obsTTPairs;
}
template <typename T>
const CCSeparatingPlanes<T>& CollisionHandler<T>::getObsCCPlanes() const {
  return _obsCCPlanes;
}
template <typename T>
CCSeparatingPlanes<T>& CollisionHandler<T>::getObsCCPlanes() {
  return _obsCCPlanes;
}
template <typename T>
const std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionHandler<T>::getSelfTTPairs() const {
  return _selfTTPairs;
}
template <typename T>
std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionHandler<T>::getSelfTTPairs() {
  return _selfTTPairs;
}
template <typename T>
const CCSeparatingPlanes<T>& CollisionHandler<T>::getSelfCCPlanes() const {
  return _selfCCPlanes;
}
template <typename T>
CCSeparatingPlanes<T>& CollisionHandler<T>::getSelfCCPlanes() {
  return _selfCCPlanes;
}
template <typename T>
const std::unordered_map<T,typename CollisionHandler<T>::GradInfo>& CollisionHandler<T>::getInfoLookup() const {
  return _infoLookup;
}
template <typename T>
std::unordered_map<T,typename CollisionHandler<T>::GradInfo>& CollisionHandler<T>::getInfoLookup() {
  return _infoLookup;
}
template <typename T>
const ThetaTrajectory<T>& CollisionHandler<T>::getThetaTrajectory() const {
  return _thetaTrajectory;
}
template <typename T>
const typename CollisionHandler<T>::Vec& CollisionHandler<T>::getControlPoints() const {
  return _controlPoints;
}
template <typename T>
void CollisionHandler<T>::setControlPoints(const Vec& x) {
  _controlPoints=x;
}
template <typename T>
void CollisionHandler<T>::setCCDObs(bool CCDObs) {
  _CCDObs=CCDObs;
}
template <typename T>
void CollisionHandler<T>::setCCDSelf(bool CCDSelf) {
  _CCDSelf=CCDSelf;
}
template <typename T>
T CollisionHandler<T>::d0() const {
  return _d0;
}
template <typename T>
T CollisionHandler<T>::d0Sqr() const {
  return _d0Sqr;
}
template <typename T>
T CollisionHandler<T>::l2() const {
  return _l2;
}
template <typename T>
T CollisionHandler<T>::eta() const {
  return _eta;
}
template <typename T>
T CollisionHandler<T>::x0() const {
  return _x0;
}
//get L1
template <typename T>
void CollisionHandler<T>::initL1() {
  for(int i=0; i<_body->nrJ(); ++i)
    _body->joint(i).initL1(*_body);
}
template <typename T>
void CollisionHandler<T>::debugL1(int res) {
  //randomize solution
  _infoLookup.clear();
  _controlPoints.setRandom();
  //choose time span
  T t0=rand()/(T)RAND_MAX*_thetaTrajectory.getNumSegment();
  T t1=rand()/(T)RAND_MAX*_thetaTrajectory.getNumSegment();
  if(t1<t0)
    std::swap(t0,t1);
  std::cout<<"Testing maxGrad for time segment: ["<<t0<<","<<t1<<"]"<<std::endl;
  for(int jointId=0; jointId<_body->nrJ(); jointId++) {
    if(!_body->joint(jointId)._mesh)
      continue;
    std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jointId)._mesh);
    int vid=rand()%mesh->vss().size();
    T L1=getVertexL1(jointId,vid,t0,t1),L1Ref=0;
    Mat3X4T DTG[3];
    for(int d=0; d<3; d++) {
      DTG[d].setZero();
      for(int c=0; c<4; c++)
        DTG[d](d,c)=c<3? (T)mesh->vss()[vid][c]: (T)1;
    }
    //we need to check that L1 is the upper bound of gradient
    for(int i=0; i<res; i++) {
      T alpha=(i+0.5)/res,t=interp1D(t0,t1,alpha);
      //compute gradient
      Vec3T grad=Vec3T::Zero();
      Vec DThetaDt=_thetaTrajectory.getDerivative(_controlPoints,t,1);
      Vec theta=_thetaTrajectory.getPoint(_controlPoints,t);
      PBDArticulatedGradientInfo<T> info(*_body,theta);
      for(int d=0; d<3; d++)
        info.DTG(jointId,*_body,DTG[d],[&](int row,T val) {
        grad[d]+=DThetaDt[row]*val;
      });
      L1Ref=std::max<T>(L1Ref,grad.norm());
    }
    std::cout<<"JID="<<jointId<<" L1Ref["<<jointId<<","<<vid<<"]="<<L1Ref<<" L1["<<jointId<<","<<vid<<"]="<<L1<<std::endl;
    ASSERT_MSGV(L1Ref<L1+Epsilon<T>::defaultEps(),"L1Ref(%f)>=L1(%f)",(double)L1Ref,(double)L1)
  }
}
template <typename T>
T CollisionHandler<T>::getConvexHullL1(int jointId,T t0,T t1) const {
  return (_body->joint(jointId)._L1.template cast<T>().transpose()*_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1)).maxCoeff();
}
template <typename T>
T CollisionHandler<T>::getVertexL1(int jointId,int j,T t0,T t1) const {
  return getVertexL1(jointId,j,_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1));
}
template <typename T>
T CollisionHandler<T>::getVertexL1(int jointId,int j,const Vec& maxGrad) const {
  return _body->joint(jointId).getL1(j).template cast<T>().dot(maxGrad);
}
template <typename T>
T CollisionHandler<T>::getEdgeL1(int jointId,int tid,int i,T t0,T t1) const {
  ASSERT_MSG(i>=0 && i<3,"The starting Vertex ID of an edge in a triangle must be >=0 and <3!")
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jointId)._mesh);
  Vec maxGrad=_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1);
  T L11=getVertexL1(jointId,mesh->iss()[tid][i],maxGrad);
  T L12=getVertexL1(jointId,mesh->iss()[tid][(i+1)%3],maxGrad);
  return std::max({L11,L12});
}
template <typename T>
T CollisionHandler<T>::getFaceL1(int jointId,int tid,T t0,T t1) const {
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jointId)._mesh);
  Vec maxGrad=_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1);
  T L11=getVertexL1(jointId,mesh->iss()[tid][0],maxGrad);
  T L12=getVertexL1(jointId,mesh->iss()[tid][1],maxGrad);
  T L13=getVertexL1(jointId,mesh->iss()[tid][2],maxGrad);
  return std::max({L11,L12,L13});
}
//get Phi
template <typename T>
T CollisionHandler<T>::getConvexHullPhi(int bvhBodyOffset) const {
  T timeDiff=_bvhBody[bvhBodyOffset]._cell.getTimeDiff();
  int jointId=_bvhBody[bvhBodyOffset]._cell._jid;
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jointId)._mesh);
  return getConvexHullL1(jointId,
                         _bvhBody[bvhBodyOffset]._cell._timeFrom,
                         _bvhBody[bvhBodyOffset]._cell._timeTo)*timeDiff+_l2*pow(timeDiff,_eta);
}
template <typename T>
T CollisionHandler<T>::getVertexPhi(int bvhBodyOffset,int vertexId) const {
  T timeDiff=_bvhBody[bvhBodyOffset]._cell.getTimeDiff();
  int jointId=_bvhBody[bvhBodyOffset]._cell._jid;
  int triangleId=_bvhBody[bvhBodyOffset]._cell._tid;
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jointId)._mesh);
  return getVertexL1(jointId,mesh->iss()[triangleId][vertexId],
                     _bvhBody[bvhBodyOffset]._cell._timeFrom,
                     _bvhBody[bvhBodyOffset]._cell._timeTo)*timeDiff+_l2*pow(timeDiff,_eta);
}
template <typename T>
T CollisionHandler<T>::getEdgePhi(int bvhBodyOffset,int startVertexId) const {
  T timeDiff=_bvhBody[bvhBodyOffset]._cell.getTimeDiff();
  return getEdgeL1(_bvhBody[bvhBodyOffset]._cell._jid,
                   _bvhBody[bvhBodyOffset]._cell._tid,
                   startVertexId,
                   _bvhBody[bvhBodyOffset]._cell._timeFrom,
                   _bvhBody[bvhBodyOffset]._cell._timeTo)*timeDiff+_l2*pow(timeDiff,_eta);
}
template <typename T>
T CollisionHandler<T>::getFacePhi(int bvhBodyOffset) const {
  T timeDiff=_bvhBody[bvhBodyOffset]._cell.getTimeDiff();
  return getFaceL1(_bvhBody[bvhBodyOffset]._cell._jid,
                   _bvhBody[bvhBodyOffset]._cell._tid,
                   _bvhBody[bvhBodyOffset]._cell._timeFrom,
                   _bvhBody[bvhBodyOffset]._cell._timeTo)*timeDiff+_l2*pow(timeDiff,_eta);
}
template <typename T>
bool CollisionHandler<T>::penetrated() const {
  return !_penetratedOffsets.empty();
}
template <typename T>
BBoxExact CollisionHandler<T>::computeBB(const EntityId<T>& id,T t0,T t1) {
  T L1,phi;
  BBoxExact bb=id.computeBB(getInfo(id.getTimeAvg(),true));
  if(id._tid>=0)
    L1=getFaceL1(id._jid,id._tid,t0,t1);
  else L1=getConvexHullL1(id._jid,t0,t1);
  phi=_d0+std::max<T>(_x0,L1*id.getTimeDiff()+_l2*pow(id.getTimeDiff(),_eta));
  return bb.enlarged(BBoxExact::Vec3T::Constant((BBoxExact::T)phi));
}
//helper
template <typename T>
const typename CollisionHandler<T>::GradInfo& CollisionHandler<T>::getInfo(T timeAvg,bool mustExist) const {
  auto it=_infoLookup.find(timeAvg);
  OMP_CRITICAL_
  if(it==_infoLookup.end()) {
    ASSERT_MSGV(!mustExist,"User requires info at time=%f must exist, but it does not!",(double)timeAvg)
    const_cast<std::unordered_map<T,GradInfo>&>(_infoLookup)[timeAvg].reset(*_body,_thetaTrajectory.getPoint(_controlPoints,timeAvg));
    it=_infoLookup.find(timeAvg);
  }
  return it->second;
}
template <typename T>
bool CollisionHandler<T>::excludedFromCCDSelf(int jid1,int jid2) const {
  //check if both joints are robot joints
  if(jid1==-1 || jid2==-1)
    return false;
  //exclude same joint
  if(jid1==jid2)
    return true;
  //exclude immediate parent joints
  if(_body->joint(jid1)._parent==jid2)
    return true;
  if(_body->joint(jid2)._parent==jid1)
    return true;
  //check user skip flags
  if(jid1>jid2)
    std::swap(jid1,jid2);
  if(_skipJIDPairs.find(jid1)!=_skipJIDPairs.end())
    if(_skipJIDPairs.at(jid1).find(jid2)!=_skipJIDPairs.at(jid1).end())
      return true;
  //check if joints are connected by only fixed joints
  while(jid1!=jid2)
    if(jid1>jid2) {
      if(_body->joint(jid1)._typeJoint!=Joint::FIX_JOINT)
        return false;
      jid1=_body->joint(jid1)._parent;
    } else {
      if(_body->joint(jid2)._typeJoint!=Joint::FIX_JOINT)
        return false;
      jid2=_body->joint(jid2)._parent;
    }
  return true;
}
//instance
template class CollisionHandler<FLOAT>;
}
