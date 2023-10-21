#include "PenetrationDetector.h"
#include "CollisionHandler.h"
#include "CollisionBarrierEnergy.h"
#include "ConvexHullDistanceEnergy.h"
#include "GJK.h"
#include <Environment/ConvexDecomposition.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Utils/SpatialRotationUtils.h>
#include <Utils/VTKWriter.h>
#include <Utils/Timing.h>
#include <Utils/Interp.h>
#include <Utils/Utils.h>
#include <unordered_set>
#include <chrono>
#include <stack>

namespace PHYSICSMOTION {
//PDEntrySpan
template <typename T>
PDEntrySpan<T>::PDEntrySpan():_jidA(-1),_jidB(-1),_PD(-1),_L1A(0),_L1B(0),_timeFrom(0),_timeTo(0) {}
template <typename T>
T PDEntrySpan<T>::getTimeAvg() const {
  return (_timeTo+_timeFrom)/2;
}
template <typename T>
T PDEntrySpan<T>::getTimeDiff() const {
  return _timeTo-_timeFrom;
}
template <typename T>
bool PDEntrySpan<T>::operator<(const PDEntrySpan<T>& other) const {
  return getTimeAvg()<other.getTimeAvg();
}
template <typename T>
bool PDEntrySpan<T>::equal(const PDEntrySpan<T>& other) const {
  return _jidA==other._jidA && _jidB==other._jidB &&
         _pBL==other._pBL && _timeFrom==other._timeFrom && _timeTo==other._timeTo;
}
template <typename T>
void PDEntrySpan<T>::writeVTK(VTKWriter<double>& os,const GradInfo& info) const {
  os.setRelativeIndex();
  std::vector<Eigen::Matrix<double,3,1>> vss;
  if(_jidA==-1)
    vss.push_back(_pAL.template cast<double>());
  else vss.push_back((ROTI(info._info._TM,_jidA)*_pAL+CTRI(info._info._TM,_jidA)).template cast<double>());
  if(_jidB==-1)
    vss.push_back(_pBL.template cast<double>());
  else vss.push_back((ROTI(info._info._TM,_jidB)*_pBL+CTRI(info._info._TM,_jidB)).template cast<double>());
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(0,0,1),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(1,0,1),
                 VTKWriter<double>::LINE,true);
}
template <typename T>
void PDEntrySpan<T>::print() const {
  std::cout<<"PDEntrySpan: timeFrom="<<_timeFrom<<" timeTo="<<_timeTo<<" jidA="<<_jidA<<" jidB"<<_jidB<<" PD="<<_PD<<std::endl;
}
//PenetrationDetector
template <typename T>
PenetrationDetector<T>::PenetrationDetector(T d0,T epsTime,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD)
  :_d0(d0),_epsTime(epsTime),_DCDObs(true),_DCDSelf(true),_useSDF(false),_skipJIDPairs(skipSelfCCD) {}
template <typename T>
PenetrationDetector<T>::PenetrationDetector(std::shared_ptr<ArticulatedBody> body,int order,int totalTime,T d0,T epsTime,T randomInitialize,bool neutralInit,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD)
  :PenetrationDetector<T>(d0,epsTime,skipSelfCCD) {
  _body=body;
  initL1();
  _thetaTrajectory=ThetaTrajectory<T>(_body->nrDOF(),order,totalTime);
  _controlPoints=CollisionHandler<T>::initializeControlPoints(randomInitialize,neutralInit,body,_thetaTrajectory);
}
template <typename T>
bool PenetrationDetector<T>::read(std::istream& is,IOData* dat) {
  registerType<Node<EntityId<T>,BBoxExact>>(dat);
  registerType<ArticulatedBody>(dat);
  registerType<PointCloudExact>(dat);
  registerType<MeshExact>(dat);
  readBinaryData(_body,is,dat);
  readBinaryData(_bodyPoints,is,dat);
  readBinaryData(_obstacles,is,dat);
  readBinaryData(_obstaclePoints,is,dat);
  readBinaryData(_bvhObstacle,is,dat);
  readBinaryData(_d0,is);
  readBinaryData(_epsTime,is);
  readBinaryData(_DCDObs,is);
  readBinaryData(_DCDSelf,is);
  readBinaryData(_useSDF,is);
  readBinaryData(_controlPoints,is);
  readBinaryData(_thetaTrajectory,is);
  readBinaryData(_skipJIDPairs,is);
  assembleObsBVH();
  return is.good();
}
template <typename T>
bool PenetrationDetector<T>::write(std::ostream& os,IOData* dat) const {
  registerType<Node<EntityId<T>,BBoxExact>>(dat);
  registerType<ArticulatedBody>(dat);
  registerType<PointCloudExact>(dat);
  registerType<MeshExact>(dat);
  writeBinaryData(_body,os,dat);
  writeBinaryData(_bodyPoints,os,dat);
  writeBinaryData(_obstacles,os,dat);
  writeBinaryData(_obstaclePoints,os,dat);
  writeBinaryData(_bvhObstacle,os,dat);
  writeBinaryData(_d0,os);
  writeBinaryData(_epsTime,os);
  writeBinaryData(_DCDObs,os);
  writeBinaryData(_DCDSelf,os);
  writeBinaryData(_useSDF,os);
  writeBinaryData(_controlPoints,os);
  writeBinaryData(_thetaTrajectory,os);
  writeBinaryData(_skipJIDPairs,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> PenetrationDetector<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new PenetrationDetector());
}
template <typename T>
std::string PenetrationDetector<T>::type() const {
  return typeid(PenetrationDetector).name();
}
template <typename T>
void PenetrationDetector<T>::setRobot(std::shared_ptr<ArticulatedBody> body,bool randomInitialize,bool neutralInit) {
  _body=body;
  initL1();
  _thetaTrajectory=ThetaTrajectory<T>(_body->nrDOF(),_thetaTrajectory.getOrder(),_thetaTrajectory.getNumSegment());
  _controlPoints=CollisionHandler<T>::initializeControlPoints(randomInitialize,neutralInit,body,_thetaTrajectory);
}
template <typename T>
void PenetrationDetector<T>::addRobot(std::shared_ptr<ArticulatedBody> body,const Vec& DOF) {
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
void PenetrationDetector<T>::addMesh(std::shared_ptr<ShapeExact> shape,const Mat3X4T& trans) {
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
void PenetrationDetector<T>::addObject(const std::string& path,T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHull) {
  //ConvexDecomposition cd(path,(double)scale,pos.template cast<double>(),rot.template cast<double>(),maxConvexHull);
  //for(auto& ob:cd.getConvexHulls())
  //  _obstacles.push_back(ob);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(path,false));
  Mat3X4T trans=Mat3X4T::Identity();
  CTR(trans)=pos;
  ROT(trans)=rot;
  obstacle->scale((GEOMETRY_SCALAR) scale);
  obstacle->transform(trans.template cast<GEOMETRY_SCALAR>());
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void PenetrationDetector<T>::addSphere(T r,int res,const Vec3T& pos) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addSphere(vss,iss,pos.template cast<double>(),(double) r,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void PenetrationDetector<T>::addCuboid(T l,T w,T h,const Vec3T& pos,const Mat3T& R) {
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
void PenetrationDetector<T>::addCapsule(T l,T w,T h,Vec3T pos,T radius,int res) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addCapsule(vss,iss,
                            Vec3T(-l/2,-w/2,-h/2).template cast<double>(),
                            Vec3T(l/2,w/2,h/2).template cast<double>(),(double) radius,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  obstacle->translate(pos.template cast<MeshExact::T>());
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
void PenetrationDetector<T>::addCapsule(const Vec3T& a,const Vec3T& b,T radius,int res) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  PHYSICSMOTION::addCapsule(vss,iss,a.template cast<double>(),b.template cast<double>(),(double) radius,res);
  std::shared_ptr<MeshExact> obstacle(new MeshExact(vss,iss,true));
  _obstacles.push_back(obstacle);
//  assembleObsBVH();
}
template <typename T>
std::vector<std::shared_ptr<PointCloudExact>> PenetrationDetector<T>::getBodyPoints() const {
  return _bodyPoints;
}
template <typename T>
std::vector<std::shared_ptr<MeshExact>> PenetrationDetector<T>::getObstacles() const {
  return _obstacles;
}
template <typename T>
std::vector<std::shared_ptr<PointCloudExact>> PenetrationDetector<T>::getObstaclePoints() const {
  return _obstaclePoints;
}
template <typename T>
std::shared_ptr<ArticulatedBody> PenetrationDetector<T>::getBody() const {
  return _body;
}
template <typename T>
void PenetrationDetector<T>::assembleBodyBVH() {
  _bodyPoints.assign(_body->nrJ(),std::shared_ptr<PointCloudExact>());
  for(int jid=0; jid<_body->nrJ(); ++jid) {
    if(!_body->joint(jid)._mesh)
      continue;
    std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(_body->joint(jid)._mesh);
    _bodyPoints[jid].reset(new PointCloudExact(*mesh,(double)_d0));
  }
}
template <typename T>
void PenetrationDetector<T>::assembleObsBVH() {
  _bvhObstacle.clear();
  _obstaclePoints.resize(_obstacles.size());
  //Just for BVH,all the triangle information is lost
  for(std::size_t i=0; i<_obstacles.size(); ++i) {
    _bvhObstacle.emplace_back();
    _bvhObstacle.back()._bb=_obstacles[i]->getBB();
    _bvhObstacle.back()._cell=i;
    _bvhObstacle.back()._nrCell=1;
    //only create point cloud and SDF, if we haven't already done so.
    if(!_obstaclePoints[i]) {
      std::cout << "Sampling point cloud for " << i+1 << "/" << _obstacles.size() << " obstacle!" << std::endl;
      _obstaclePoints[i].reset(new PointCloudExact(*_obstacles[i],(double)_d0));
      VTKWriter<double> os("pointCloud","pointCloud"+std::to_string(i)+".vtk",true);
      _obstaclePoints[i]->writeVTK(os,Mat3X4T::Identity().template cast<MeshExact::T>());
      _obstaclePoints[i]->writeSDFVTK("pointCloudSDF"+std::to_string(i)+".vtk");
    }
  }
  Node<int,BBoxExact>::buildBVHBottomUpAll(_bvhObstacle);
}
template <typename T>
void PenetrationDetector<T>::reset() {
  _obstacles.clear();
  _obstaclePoints.clear();
//  assembleObsBVH();
}
//detect penetration
template <typename T>
bool PenetrationDetector<T>::DCD(bool clear,bool BF,bool output) {
  _infoLookup.clear();
  if(clear)
    _PDEntries.clear();
  //obstacle PD
  if(_DCDObs) {
    OMP_PARALLEL_FOR_
    for(int jid=0; jid<_body->nrJ(); ++jid) {
      PDEntrySpan<T> entry;
      if(output) {
        OMP_CRITICAL_
        std::cout << "Detecting contact between: jid=" << jid << " and obstacle" << " #PD=" << _PDEntries.getVector().size() << "!" << std::endl;
      }
      entry=DCD(jid,-1,BF);
      if(entry._PD>0) {
        entry._jidB=-1;
        //if(duplicate(entry))
        //  return false;
        _PDEntries.push_back(entry);
      }
    }
  }
  //self PD
  if(_DCDSelf) {
    OMP_PARALLEL_FOR_
    for(int off=0; off<_body->nrJ()*_body->nrJ(); ++off) {
      PDEntrySpan<T> entry;
      int jid=off%_body->nrJ();
      int jid2=off/_body->nrJ();
      if(!excludedFromDCDSelf(jid,jid2)) {
        if(output) {
          OMP_CRITICAL_
          std::cout << "Detecting contact between: jid=" << jid << " jid2=" << jid2 << " #PD=" << _PDEntries.getVector().size() << "!" << std::endl;
        }
        entry=DCD(jid,jid2,BF);
        if(entry._PD>0) {
          //if(duplicate(entry))
          //  return false;
          _PDEntries.push_back(entry);
        }
      }
    }
  }
  //duplicate check
  bool hasDuplicate=false;
  const auto& PDEntries=_PDEntries.getVector();
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)PDEntries.size(); i++)
    for(int j=0; !hasDuplicate && j<i; j++)
      if(PDEntries[i].equal(PDEntries[j]))
        hasDuplicate=true;
  std::cout<<"DCD finished"<<std::endl;
  return !hasDuplicate;
}
template <typename T>
PDEntrySpan<T> PenetrationDetector<T>::DCD(int jid,int jid2,bool BF) {
  if(!std::dynamic_pointer_cast<MeshExact>(_body->joint(jid)._mesh))
    return PDEntrySpan<T>();
  if(jid2>=0 && !_bodyPoints[jid2])
    return PDEntrySpan<T>();
  if(jid2<0 && _bvhObstacle.empty())
    return PDEntrySpan<T>();
  //main loop
  T rad;
  std::stack<PDEntrySpan<T>> ss;
  PDEntrySpan<T> maxPD,entry,entryL,entryR;
  ss.push(queryPD(0,_thetaTrajectory.getNumSegment(),jid,jid2));
  maxPD._PD=-std::numeric_limits<double>::max();
  while(!ss.empty()) {
    //std::cout<<ss.size()<<std::endl;
    entry=ss.top();
    ss.pop();
    if(entry._PD>maxPD._PD)
      maxPD=entry;
    if(entry._timeTo-entry._timeFrom<_epsTime)
      continue;
    rad=(entry._L1A+entry._L1B)*entry.getTimeDiff()/2;
    if(!BF && entry._PD+rad<maxPD._PD)
      continue;
    else {
      //first explore node with larger PD
      entryL=queryPD(entry._timeFrom,entry.getTimeAvg(),jid,jid2);
      entryR=queryPD(entry.getTimeAvg(),entry._timeTo,jid,jid2);
      if(entryL._PD>entryR._PD)
        std::swap(entryL,entryR);
      ss.push(entryL);
      ss.push(entryR);
    }
  }
  return maxPD;
}
template <typename T>
PDEntrySpan<T> PenetrationDetector<T>::queryPD(T timeFrom,T timeTo,int jid,int jid2) {
  if(jid2==-1)
    return queryPDObs(timeFrom,timeTo,jid);
  else return queryPDSelf(timeFrom,timeTo,jid,jid2);
}
template <typename T>
PDEntrySpan<T> PenetrationDetector<T>::queryPDSelf(T timeFrom,T timeTo,int jid,int jid2) {
  PDEntrySpan<T> ret;
  Vec3T pAL,pBL;
  ret._timeFrom=timeFrom;
  ret._timeTo=timeTo;
  ret._jidA=jid;
  ret._jidB=jid2;
  const CollisionGradInfo<T>& info=getInfo(ret.getTimeAvg());
  if(_useSDF)
    ret._PD=GJK::runPD<T>(_bodyPoints[ret._jidA],_bodyPoints[ret._jidB],TRANSI(info._info._TM,ret._jidA),TRANSI(info._info._TM,ret._jidB),pAL,pBL);
  else ret._PD=GJK::runPD<T>(std::dynamic_pointer_cast<MeshExact>(_body->joint(ret._jidA)._mesh),_bodyPoints[ret._jidB],TRANSI(info._info._TM,ret._jidA),TRANSI(info._info._TM,ret._jidB),pAL,pBL);
  ret._L1A=getConvexHullL1(ret._jidA,ret._timeFrom,ret._timeTo);
  ret._L1B=getConvexHullL1(ret._jidB,ret._timeFrom,ret._timeTo);
  ret._pAL=pAL.template cast<T>();
  ret._pBL=pBL.template cast<T>();
  return ret;
}
template <typename T>
PDEntrySpan<T> PenetrationDetector<T>::queryPDObs(T timeFrom,T timeTo,int jid,int jid2) {
  PDEntrySpan<T> ret;
  Vec3T pAL,pBL;
  ret._timeFrom=timeFrom;
  ret._timeTo=timeTo;
  ret._jidA=jid;
  ret._jidB=jid2;
  const CollisionGradInfo<T>& info=getInfo(ret.getTimeAvg());
  if(_bvhObstacle[ret._jidB]._cell==-1) {
    //internal node, just computing PD for pruning
    Vec3T n,normal,ctr;
    Mat3T hessian;
    Eigen::Matrix<int,2,1> feat;
    Mat3X4T trans=TRANSI(info._info._TM,ret._jidA),invTrans;
    INV(invTrans,trans);
    ctr=(_bvhObstacle[ret._jidB]._bb.minCorner()+_bvhObstacle[ret._jidB]._bb.maxCorner()).template cast<T>()/2;
    if(_useSDF)
      ret._PD=-_bodyPoints[ret._jidA]->template closestSDF<T>(ROT(invTrans)*ctr+CTR(invTrans),normal);
    else ret._PD=-std::dynamic_pointer_cast<MeshExact>(_body->joint(ret._jidA)._mesh)->template closest<T>(ROT(invTrans)*ctr+CTR(invTrans),n,normal,hessian,feat);
    ret._L1A=getConvexHullL1(ret._jidA,ret._timeFrom,ret._timeTo);
    ret._L1B=0;
  } else {
    //leaf node, compute PD and pAL,pBL,L1A,L1B
    if(_useSDF)
      ret._PD=GJK::runPD<T>(_bodyPoints[ret._jidA],_obstaclePoints[_bvhObstacle[ret._jidB]._cell],TRANSI(info._info._TM,ret._jidA),Mat3X4T::Identity(),pAL,pBL);
    else ret._PD=GJK::runPD<T>(std::dynamic_pointer_cast<MeshExact>(_body->joint(ret._jidA)._mesh),_obstaclePoints[_bvhObstacle[ret._jidB]._cell],TRANSI(info._info._TM,ret._jidA),Mat3X4T::Identity(),pAL,pBL);
    ret._L1A=getConvexHullL1(ret._jidA,ret._timeFrom,ret._timeTo);
    ret._L1B=0;
    ret._pAL=pAL.template cast<T>();
    ret._pBL=pBL.template cast<T>();
  }
  return ret;
}
template <typename T>
PDEntrySpan<T> PenetrationDetector<T>::queryPDObs(T timeFrom,T timeTo,int jid) {
  T rad;
  std::stack<PDEntrySpan<T>> ss;
  PDEntrySpan<T> maxPD,entry,entryL,entryR;
  ss.push(queryPDObs(timeFrom,timeTo,jid,(int) _bvhObstacle.size()-1));
  maxPD._PD=-std::numeric_limits<double>::max();
  while(!ss.empty()) {
    entry=ss.top();
    ss.pop();
    if(_bvhObstacle[entry._jidB]._cell>=0) {
      if(entry._PD>maxPD._PD)
        maxPD=entry;
    } else {
      const Node<int,BBoxExact>& node=_bvhObstacle[entry._jidB];
      rad=(node._bb.maxCorner()-node._bb.minCorner()).template cast<T>().norm()/2;
      if(entry._PD+rad<maxPD._PD)
        continue;
      else {
        //first explore node with smaller PD
        entryL=queryPDObs(timeFrom,timeTo,jid,node._l);
        entryR=queryPDObs(timeFrom,timeTo,jid,node._r);
        if(entryL._PD>entryR._PD)
          std::swap(entryL,entryR);
        ss.push(entryL);
        ss.push(entryR);
      }
    }
  }
  return maxPD;
}
template <typename T>
void PenetrationDetector<T>::writeVTK(const std::string& path) {
  //timestamp
  std::stack<PDEntrySpan<T>> ss;
  PDEntrySpan<T> entry,entryL,entryR;
  std::vector<PDEntrySpan<T>> entries;
  {
    entry._timeFrom=0;
    entry._timeTo=_thetaTrajectory.getNumSegment();
    ss.push(entry);
  }
  while(!ss.empty()) {
    entry=ss.top();
    ss.pop();
    if(entry.getTimeDiff()<_epsTime)
      entries.push_back(entry);
    else {
      entryL._timeFrom=entry._timeFrom;
      entryL._timeTo=entry.getTimeAvg();
      ss.push(entryL);
      entryR._timeFrom=entry.getTimeAvg();
      entryR._timeTo=entry._timeTo;
      ss.push(entryR);
    }
  }
  std::sort(entries.begin(),entries.end());
  //write obstacle
  recreate(path);
  {
    VTKWriter<double> os("obstacle",path+"/obstacle.vtk",true);
    for(auto mesh:_obstacles)
      mesh->writeVTK(os,Mat3X4T::Identity().template cast<MeshExact::T>());
  }
  //write robot
  _infoLookup.clear();
  std::sort(_PDEntries.begin(),_PDEntries.end());
  for(int iter=0,iter2=0; iter<(int) entries.size(); iter++) {
    const GradInfo& info=getInfo(entries[iter].getTimeAvg());
    VTKWriter<double> os("body",path+"/iter"+std::to_string(iter)+".vtk",true);
    _body->writeVTK(os,info._info._TM.template cast<ArticulatedBody::T>());
    if(iter2<(int) _PDEntries.getVector().size() && entries[iter].getTimeAvg()==_PDEntries.getVector()[iter2].getTimeAvg()) {
      std::cout << "Writing iter=" << iter << " PDIter=" << iter2 << std::endl;
      _PDEntries.getVector()[iter2].writeVTK(os,info);
      iter2++;
    }
  }
}
template <typename T>
void PenetrationDetector<T>::debugDCDBF() {
  DCD(true,false);
  typename ParallelVector<PDEntrySpan<T>>::vector_type PDEntriesNonBF=getPDEntries().getVector();
  DCD(true,true);
  typename ParallelVector<PDEntrySpan<T>>::vector_type PDEntriesBF=getPDEntries().getVector();
  //compare
  std::cout << "#PDNonBF=" << PDEntriesNonBF.size() << " #PDBF=" << PDEntriesBF.size() << std::endl;
  ASSERT_MSGV(PDEntriesNonBF.size()==PDEntriesBF.size(),
              "Size of PDEntriesNonBF(%d)!=size of PDEntriesBF(%d)!",
              (int) PDEntriesNonBF.size(),(int) PDEntriesBF.size())
  for(int i=0; i<(int) PDEntriesNonBF.size(); i++) {
    std::cout << "timeAvgNonBF=" << PDEntriesNonBF[i].getTimeAvg() << " timeAvgBF=" << PDEntriesBF[i].getTimeAvg() << std::endl;
    ASSERT_MSGV(PDEntriesNonBF[i]._jidA==PDEntriesBF[i]._jidA,
                "PDEntriesNonBF[%d]._jidA(%d)!=size of PDEntriesBF[%d]._jidA(%d)!",
                i,PDEntriesNonBF[i]._jidA,i,PDEntriesBF[i]._jidA)
    ASSERT_MSGV(PDEntriesNonBF[i]._jidB==PDEntriesBF[i]._jidB,
                "PDEntriesNonBF[%d]._jidB(%d)!=size of PDEntriesBF[%d]._jidB(%d)!",
                i,PDEntriesNonBF[i]._jidB,i,PDEntriesBF[i]._jidB)
    ASSERT_MSGV(PDEntriesNonBF[i].getTimeAvg()==PDEntriesBF[i].getTimeAvg(),
                "PDEntriesNonBF[%d].timeAvg(%f)!=size of PDEntriesBF[%d].timeAvg(%f)!",
                i,(double)PDEntriesNonBF[i].getTimeAvg(),i,(double)PDEntriesBF[i].getTimeAvg())
  }
}
//info query
template <typename T>
ParallelVector<PDEntrySpan<T>>& PenetrationDetector<T>::getPDEntries() {
  return _PDEntries;
}
template <typename T>
const ParallelVector<PDEntrySpan<T>>& PenetrationDetector<T>::getPDEntries() const {
  return _PDEntries;
}
template <typename T>
const ThetaTrajectory<T>& PenetrationDetector<T>::getThetaTrajectory() const {
  return _thetaTrajectory;
}
template <typename T>
const typename PenetrationDetector<T>::Vec& PenetrationDetector<T>::getControlPoints() const {
  return _controlPoints;
}
template <typename T>
void PenetrationDetector<T>::setControlPoints(const Vec& x) {
  _infoLookup.clear();
  _controlPoints=x;
}
template <typename T>
void PenetrationDetector<T>::setDCDObs(bool DCDObs) {
  _DCDObs=DCDObs;
}
template <typename T>
void PenetrationDetector<T>::setDCDSelf(bool DCDSelf) {
  _DCDSelf=DCDSelf;
}
template <typename T>
void PenetrationDetector<T>::setUseSDF(bool useSDF) {
  _useSDF=useSDF;
}
template <typename T>
bool PenetrationDetector<T>::getDCDObs() const {
  return _DCDObs;
}
template <typename T>
bool PenetrationDetector<T>::getDCDSelf() const {
  return _DCDSelf;
}
template <typename T>
bool PenetrationDetector<T>::getUseSDF() const {
  return _useSDF;
}
template <typename T>
T PenetrationDetector<T>::d0() const {
  return _d0;
}
template <typename T>
T PenetrationDetector<T>::epsTime() const {
  return _epsTime;
}
//get L1
template <typename T>
void PenetrationDetector<T>::initL1() {
  for(int i=0; i<_body->nrJ(); ++i)
    _body->joint(i).initL1(*_body);
}
template <typename T>
void PenetrationDetector<T>::debugL1(int res) {
  //randomize solution
  _infoLookup.clear();
  _controlPoints.setRandom();
  //choose time span
  T t0=rand()/(T) RAND_MAX*_thetaTrajectory.getNumSegment();
  T t1=rand()/(T) RAND_MAX*_thetaTrajectory.getNumSegment();
  if(t1<t0)
    std::swap(t0,t1);
  std::cout << "Testing maxGrad for time segment: [" << t0 << "," << t1 << "]" << std::endl;
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
        DTG[d](d,c)=c<3? (T) mesh->vss()[vid][c]: (T) 1;
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
    std::cout << "JID=" << jointId << " L1Ref[" << jointId << "," << vid << "]=" << L1Ref << " L1[" << jointId << "," << vid << "]=" << L1 << std::endl;
    ASSERT_MSGV(L1Ref<L1+Epsilon<T>::defaultEps(),"L1Ref(%f)>=L1(%f)",(double) L1Ref,(double) L1)
  }
}
template <typename T>
T PenetrationDetector<T>::getConvexHullL1(int jointId,T t0,T t1) const {
  if(t0==t1)
    return 0;
  return (_body->joint(jointId)._L1.template cast<T>().transpose()*_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1)).maxCoeff();
}
template <typename T>
T PenetrationDetector<T>::getVertexL1(int jointId,int j,T t0,T t1) const {
  return getVertexL1(jointId,j,_thetaTrajectory.getMaxGrad(_controlPoints,t0,t1));
}
template <typename T>
T PenetrationDetector<T>::getVertexL1(int jointId,int j,const Vec& maxGrad) const {
  return _body->joint(jointId).getL1(j).template cast<T>().dot(maxGrad);
}
//helper
template <typename T>
const typename PenetrationDetector<T>::GradInfo& PenetrationDetector<T>::getInfo(T timeAvg,bool mustExist) const {
  auto it=_infoLookup.find(timeAvg);
  OMP_CRITICAL_
  if(it==_infoLookup.end()) {
    ASSERT_MSGV(!mustExist,"User requires info at time=%f must exist, but it does not!",(double) timeAvg)
    const_cast<std::unordered_map<T,GradInfo>&>(_infoLookup)[timeAvg].reset(*_body,_thetaTrajectory.getPoint(_controlPoints,timeAvg));
    it=_infoLookup.find(timeAvg);
  }
  return it->second;
}
template <typename T>
bool PenetrationDetector<T>::excludedFromDCDSelf(int jid1,int jid2) const {
  if(jid1==-1 || jid2==-1)
    return false;
  if(jid1==jid2)
    return true;
  if(_body->joint(jid1)._parent==jid2)
    return true;
  if(_body->joint(jid2)._parent==jid1)
    return true;
  if(jid1>jid2) {
    std::swap(jid1,jid2);
  }
  if(_skipJIDPairs.find(jid1)!=_skipJIDPairs.end())
    if(_skipJIDPairs.at(jid1).find(jid2)!=_skipJIDPairs.at(jid1).end()) {
      return true;
    }
  return false;
}
template <typename T>
void PenetrationDetector<T>::print() const {
  for(const auto& entry:_PDEntries) {
    entry.print();
  }
}
//instance
template struct PDEntrySpan<FLOAT>;
template class PenetrationDetector<FLOAT>;
}
