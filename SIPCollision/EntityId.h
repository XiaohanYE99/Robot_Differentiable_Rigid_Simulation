#ifndef ENTITY_ID_H
#define ENTITY_ID_H

#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/Utils.h>

namespace PHYSICSMOTION {
template <typename T>
class GJKPolytope;
template <typename T>
struct CollisionGradInfo {
  DECL_MAT_VEC_MAP_TYPES_T
  CollisionGradInfo();
  CollisionGradInfo(const ArticulatedBody& body,const Vec& theta);
  void reset(const ArticulatedBody& body,const Vec& theta);
  //data
  MatT _HTheta;
  Mat3XT _DTG;
  std::vector<Mat3XT> _globalVss;
  std::vector<std::shared_ptr<GJKPolytope<T>>> _polytopes;
  PBDArticulatedGradientInfo<T> _info;
};
template <typename T>
struct EntityId : public SerializableBase {
  DECL_MAT_VEC_MAP_TYPES_T
  using SerializableBase::read;
  using SerializableBase::write;
  EntityId();
  EntityId(const EntityId& entityId);
  EntityId(std::shared_ptr<GJKPolytope<T>> obs,int tid);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void reset();
  EntityId& operator=(int id);
  EntityId& operator=(const EntityId& entityId);
  bool operator==(const EntityId& entityId) const;
  bool operator==(const int id) const;
  bool operator!=(const EntityId& entityId) const;
  bool operator!=(const int id) const;
  bool operator<(const EntityId& other) const;
  friend std::ostream& operator<<(std::ostream& out, const EntityId<T>& id) {
    out << "EntityId: Jid=" << id._jid << " tid=" << id._tid << " timeFrom=" << id._timeFrom << " timeTo=" << id._timeTo << " obs=" << id._obs;
    return out;
  }
  const GJKPolytope<T>& getPolytope(const CollisionGradInfo<T>& info) const;
  BBoxExact computeBB(const CollisionGradInfo<T>& info) const;
  Vec3T globalV(const CollisionGradInfo<T>& info,int d) const;
  Vec3T localV(int d) const;
  Vec3TCM point() const;
  Vec3TM point();
  bool isObstacleConvexHull() const;
  bool isRobotConvexHull() const;
  bool isObstacleTriangle() const;
  bool isRobotTriangle() const;
  bool isObstacle() const;
  void print() const;
  T getTimeAvg() const;
  T getTimeDiff() const;
  bool checkBss(int cnt) const;
  //Joint:
  //  Case 1: _jid=-1,_tid=-1,_obs==NULL means internal node
  //  Case 2: _jid>=0,_tid=-1,_obs==NULL means root of joint
  //  Case 3: _jid>=0,_tid>=0,_obs==NULL means leaf of joint, representing a triangle
  //Obstacle:
  //  Case 4: _jid=-1,_tid>=0,_obs!=NULL means leaf of obstacle
  std::shared_ptr<GJKPolytope<T>> _obs;
  std::shared_ptr<MeshExact> _link;
  T _timeFrom=-1;
  T _timeTo=-1;
  int _jid=-1;
  int _tid=-1;
  //temporary variables, these are not serialized
  T _facePhi;
  T _edgePhi[3];
  T _vertexPhi[3];
  T _convexHullPhi;
};
}
#endif
