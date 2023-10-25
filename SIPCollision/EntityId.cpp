#include "EntityId.h"
#include "ConvexHullDistanceEnergy.h"
#include <Environment/MeshExact.h>

namespace PHYSICSMOTION {
//CollisionGradInfo
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo() {}
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo(const ArticulatedBody& body,const Vec& theta) {
  reset(body,theta);
}
template <typename T>
void CollisionGradInfo<T>::reset(const ArticulatedBody& body,const Vec& theta) {
  _HTheta.setZero(body.nrDOF(),body.nrDOF());
  _DTG.setZero(3,4*body.nrJ());
  _info.reset(body,theta);
  //data
  _globalVss.resize(body.nrJ());
  _polytopes.resize(body.nrJ());
  for(int i=0; i<body.nrJ(); i++) {
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    body.joint(i)._mesh->getMesh(vss,iss);
    std::shared_ptr<MeshExact> m;
    m.reset(new MeshExact(vss,iss,false));
    //std::shared_ptr<MeshExact> m=std::dynamic_pointer_cast<MeshExact>(body.joint(i)._mesh);
    if(!m) {
      _globalVss[i].resize(3,0);
    } else {
      _globalVss[i].resize(3,m->vss().size());
      for(int j=0; j<(int)m->vss().size(); j++)
        _globalVss[i].col(j)=ROTI(_info._TM,i)*m->vss()[j].template cast<T>()+CTRI(_info._TM,i);
      _polytopes[i].reset(new GJKPolytope<T>(i,body,*this));
    }
  }
}
//EntityId
template <typename T>
EntityId<T>::EntityId() {
  reset();
}
template <typename T>
EntityId<T>::EntityId(const EntityId& other) {
  operator=(other);
}
template <typename T>
EntityId<T>::EntityId(std::shared_ptr<GJKPolytope<T>> obs,int tid) {
  reset();
  _tid=tid;
  _obs=obs;
}
template <typename T>
bool EntityId<T>::read(std::istream& is,IOData* dat) {
  registerType<GJKPolytope<T>>(dat);
  readBinaryData(_obs,is,dat);
  readBinaryData(_link,is,dat);
  readBinaryData(_timeFrom,is);
  readBinaryData(_timeTo,is);
  readBinaryData(_jid,is);
  readBinaryData(_tid,is);
  return is.good();
}
template <typename T>
bool EntityId<T>::write(std::ostream& os,IOData* dat) const {
  registerType<GJKPolytope<T>>(dat);
  writeBinaryData(_obs,os,dat);
  writeBinaryData(_link,os,dat);
  writeBinaryData(_timeFrom,os);
  writeBinaryData(_timeTo,os);
  writeBinaryData(_jid,os);
  writeBinaryData(_tid,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> EntityId<T>::copy() const {
  return std::shared_ptr<EntityId<T>>(new EntityId<T>);
}
template <typename T>
std::string EntityId<T>::type() const {
  return typeid(EntityId).name();
}
template <typename T>
void EntityId<T>::reset() {
  _obs=NULL;
  _link=NULL;
  _timeFrom=-1;
  _timeTo=-1;
  _jid=-1;
  _tid=-1;

  _facePhi=0;
  _edgePhi[0]=_edgePhi[1]=_edgePhi[2]=0;
  _vertexPhi[0]=_vertexPhi[1]=_vertexPhi[2]=0;
  _convexHullPhi=0;
}
template <typename T>
EntityId<T>& EntityId<T>::operator=(int id) {
  reset();
  ASSERT(id==-1)
  return *this;
}
template <typename T>
EntityId<T>& EntityId<T>::operator=(const EntityId& entityId) {
  _obs=entityId._obs;
  _link=entityId._link;
  _timeTo=entityId._timeTo;
  _timeFrom=entityId._timeFrom;
  _tid=entityId._tid;
  _jid=entityId._jid;

  _facePhi=entityId._facePhi;
  _edgePhi[0]=entityId._edgePhi[0];
  _edgePhi[1]=entityId._edgePhi[1];
  _edgePhi[2]=entityId._edgePhi[2];
  _vertexPhi[0]=entityId._vertexPhi[0];
  _vertexPhi[1]=entityId._vertexPhi[1];
  _vertexPhi[2]=entityId._vertexPhi[2];
  _convexHullPhi=entityId._convexHullPhi;
  return *this;
}
template <typename T>
bool EntityId<T>::operator==(const EntityId& entityId) const {
  return _obs==entityId._obs &&
         _link==entityId._link &&
         _timeFrom==entityId._timeFrom &&
         _timeTo==entityId._timeTo &&
         _jid==entityId._jid &&
         _tid==entityId._tid;
}
template <typename T>
bool EntityId<T>::operator==(const int id) const {
  ASSERT_MSG(id==-1,"EntityId::operator== only takes parameter -1!")
  //check if this is a leaf
  return !_obs && !_link;
}
template <typename T>
bool EntityId<T>::operator!=(const EntityId& entityId) const {
  return !EntityId::operator==(entityId);
}
template <typename T>
bool EntityId<T>::operator!=(const int id) const {
  return !EntityId::operator==(id);
}
template <typename T>
bool EntityId<T>::operator<(const EntityId& other) const {
  if(_obs>other._obs)
    return false;
  else if(_obs<other._obs)
    return true;

  if(_link>other._link)
    return false;
  else if(_link<other._link)
    return true;

  if(_timeFrom>other._timeFrom)
    return false;
  else if(_timeFrom<other._timeFrom)
    return true;

  if(_timeTo>other._timeTo)
    return false;
  else if(_timeTo<other._timeTo)
    return true;

  if(_jid>other._jid)
    return false;
  else if(_jid<other._jid)
    return true;

  if(_tid>other._tid)
    return false;
  else if(_tid<other._tid)
    return true;
  return false;
}
template <typename T>
const GJKPolytope<T>& EntityId<T>::getPolytope(const CollisionGradInfo<T>& info) const {
  if(isRobotConvexHull())
    return *(info._polytopes[_jid]);
  else {
    ASSERT(isObstacleConvexHull())
    return *_obs;
  }
}
template <typename T>
BBoxExact EntityId<T>::computeBB(const CollisionGradInfo<T>& info) const {
  BBoxExact ret;
  if(isRobotTriangle() || isObstacleTriangle()) {
    ret.setUnion(globalV(info,0).template cast<BBoxExact::T>());
    ret.setUnion(globalV(info,1).template cast<BBoxExact::T>());
    ret.setUnion(globalV(info,2).template cast<BBoxExact::T>());
  } else if(isRobotConvexHull()) {
    ret.minCorner()=info._globalVss[_jid].rowwise().minCoeff().template cast<BBoxExact::T>();
    ret.maxCorner()=info._globalVss[_jid].rowwise().maxCoeff().template cast<BBoxExact::T>();
  } else {
    ASSERT(isObstacleConvexHull())
    ret=_obs->mesh()->getBB();
  }
  return ret;
}
template <typename T>
typename EntityId<T>::Vec3T EntityId<T>::globalV(const CollisionGradInfo<T>& info,int d) const {
  if(isRobotTriangle()) {
    int vid=_link->iss()[_tid][d];
    return Vec3TCM(&(info._globalVss[_jid](0,vid)));
  } else {
    ASSERT(isObstacleTriangle())
    int vid=_obs->mesh()->iss()[_tid][d];
    return _obs->mesh()->vss()[vid].template cast<T>();
  }
}
template <typename T>
typename EntityId<T>::Vec3T EntityId<T>::localV(int d) const {
  if(isRobotTriangle()) {
    int vid=_link->iss()[_tid][d];
    return _link->vss()[vid].template cast<T>();
  } else {
    ASSERT(isObstacleTriangle())
    int vid=_obs->mesh()->iss()[_tid][d];
    return _obs->mesh()->vss()[vid].template cast<T>();
  }
}
template <typename T>
typename EntityId<T>::Vec3TCM EntityId<T>::point() const {
  return Vec3TCM(_vertexPhi);
}
template <typename T>
typename EntityId<T>::Vec3TM EntityId<T>::point() {
  return Vec3TM(_vertexPhi);
}
template <typename T>
bool EntityId<T>::isObstacleConvexHull() const {
  return _jid==-1 && _tid==-1 && _obs!=NULL && _link==NULL;
}
template <typename T>
bool EntityId<T>::isRobotConvexHull() const {
  return _jid>=0 && _tid==-1 && _obs==NULL && _link!=NULL;
}
template <typename T>
bool EntityId<T>::isObstacleTriangle() const {
  return _jid==-1 && _tid>=0 && _obs!=NULL && _link==NULL;
}
template <typename T>
bool EntityId<T>::isRobotTriangle() const {
  return _jid>=0 && _tid>=0 && _obs==NULL && _link!=NULL;
}
template <typename T>
bool EntityId<T>::isObstacle() const {
  return _obs!=NULL && _link==NULL;
}
template <typename T>
void EntityId<T>::print() const {
  std::cout << "EntityId: " << " timeFrom=" << _timeFrom << " timeTo=" << _timeTo
            << " jid=" << _jid << " tid=" << _tid << " obs=" << _obs << " bss=";
  for(int i=0; i<6; ++i)
    std::cout << checkBss(i);
  std::cout << std::endl;
}
template <typename T>
rational EntityId<T>::getTimeAvg() const {
  return (_timeFrom+_timeTo)/2.;
}
template <typename T>
rational EntityId<T>::getTimeDiff() const {
  return (_timeTo-_timeFrom)/2.;
}
template <typename T>
bool EntityId<T>::checkBss(int cnt) const {
  ASSERT((_obs!=NULL)^(_link!=NULL))
  if(_tid<0)
    return 0;
  else if(_obs!=NULL)
    return (_obs->mesh()->bss()[_tid]>>cnt)&1;
  else if(_link!=NULL)
    return (_link->bss()[_tid]>>cnt)&1;
  else return 0;
}
//instance
template struct CollisionGradInfo<FLOAT>;
template struct EntityId<FLOAT>;
}
