#ifndef CC_SEPARATING_PLANES_H
#define CC_SEPARATING_PLANES_H

#include "EntityId.h"
#include "CCSeparatingPlanes.h"

namespace PHYSICSMOTION {
template <typename T>
struct LSS {
  bool operator()(const std::pair<EntityId<T>,EntityId<T>>& a,const std::pair<EntityId<T>,EntityId<T>>& b) const {
    if(a.first<b.first)
      return true;
    else if(b.first<a.first)
      return false;

    if(a.second<b.second)
      return true;
    else if(b.second<a.second)
      return false;

    return false;
  }
};
template <typename T>
struct EntityIdHash {
  size_t operator()(const EntityId<T>& key) const {
    size_t seed=0;
    //std::shared_ptr<GJKPolytope<T>> _obs;
    //std::shared_ptr<MeshExact> _link;
    //T _timeFrom=-1;
    //T _timeTo=-1;
    //int _jid=-1;
    //int _tid=-1;
    hash_combine(seed,std::hash<std::shared_ptr<GJKPolytope<T>>>()(key._obs));
    hash_combine(seed,std::hash<std::shared_ptr<MeshExact>>()(key._link));
    hash_combine(seed,std::hash<T>()(key._timeFrom));
    hash_combine(seed,std::hash<T>()(key._timeTo));
    hash_combine(seed,std::hash<int>()(key._jid));
    hash_combine(seed,std::hash<int>()(key._tid));
    return seed;
  }
};
template <typename T>
struct EntityIdPairHash {
  size_t operator()(const std::pair<EntityId<T>,EntityId<T>>& key) const {
    size_t seed=0;
    EntityIdHash<T> hash;
    hash_combine(seed,hash(key.first));
    hash_combine(seed,hash(key.second));
    return seed;
  }
};
template <typename T>
struct CCSeparatingPlane : public SerializableBase {
  DECL_MAT_VEC_MAP_TYPES_T
  CCSeparatingPlane():_initialized(false),_x(Vec4T::Zero()) {}
  virtual bool read(std::istream& is,IOData*) {
    readBinaryData(_initialized,is);
    readBinaryData(_x,is);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData*) const {
    writeBinaryData(_initialized,os);
    writeBinaryData(_x,os);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new CCSeparatingPlane<T>());
  }
  virtual std::string type() const {
    return typeid(CCSeparatingPlane<T>).name();
  }
  bool _initialized;
  Vec4T _x;
};
template <typename T>
struct CCSeparatingPlanes : public SerializableBase {
  CCSeparatingPlanes():_CCUpdated(false) {}
  virtual bool read(std::istream& is,IOData* dat) {
    readBinaryData(_idToCCPairs,is,dat);
    readBinaryData(_CCPlanes,is,dat);
    readBinaryData(_CCUpdated,is);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_idToCCPairs,os,dat);
    writeBinaryData(_CCPlanes,os,dat);
    writeBinaryData(_CCUpdated,os);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new CCSeparatingPlanes<T>());
  }
  virtual std::string type() const {
    return typeid(CCSeparatingPlanes<T>).name();
  }
  void insertPlane(const std::pair<EntityId<T>,EntityId<T>>& pair) {
    ASSERT_MSG(pair.first.isRobotConvexHull(),
               "CCSeparatingPlanes only accept convex hull as first parameter!")
    ASSERT_MSG(pair.second.isRobotConvexHull() || pair.second.isObstacleConvexHull(),
               "CCSeparatingPlanes only accept convex hull as second parameter!")
    if(_CCPlanes.find(pair)==_CCPlanes.end()) {
      if(pair.first.isRobotConvexHull())
        _idToCCPairs[pair.first].insert(pair);
      if(pair.second.isRobotConvexHull())
        _idToCCPairs[pair.second].insert(pair);
      _CCPlanes[pair]=CCSeparatingPlane<T>();
      _CCUpdated=true;
    }
  }
  const std::unordered_map<std::pair<EntityId<T>,EntityId<T>>,CCSeparatingPlane<T>,EntityIdPairHash<T>>& getCCPlanes() const {
    return _CCPlanes;
  }
  //when subdividing, remove all related pairs
  void subdivide(const EntityId<T>& node) {
    std::vector<std::pair<EntityId<T>,std::pair<EntityId<T>,EntityId<T>>>> deleteCache;
    //delete from _CCPlanes
    for(auto pair:_idToCCPairs[node]) {
      //mark for delete from _idToCCPairs
      if(pair.first.isRobotConvexHull())
        deleteCache.emplace_back(pair.first,pair);
      if(pair.second.isRobotConvexHull())
        deleteCache.emplace_back(pair.second,pair);
      //delete from _CCPlanes
      _CCPlanes.erase(pair);
    }
    //delete from _idToCCPairs
    for(auto cache:deleteCache)
      _idToCCPairs[cache.first].erase(cache.second);
  }
  void clear() {
    _idToCCPairs.clear();
    _CCPlanes.clear();
    _CCUpdated=false;
  }
  int size() const {
    return (int)_CCPlanes.size();
  }
  void clearUpdateFlag() {
    _CCUpdated=false;
  }
  bool hasUpdateFlag() const {
    return _CCUpdated;
  }
  //fetch parallel
  CCSeparatingPlane<T>& getSeparatingPlane(const std::pair<EntityId<T>,EntityId<T>>& key) {
    ASSERT_MSG(_CCPlanes.find(key)!=_CCPlanes.end(),"Cannot find key in _CCPlanes!")
    return _CCPlanes.find(key)->second;
  }
  void getCCPairs(std::vector<std::pair<EntityId<T>,EntityId<T>>>& CCPairs,
                  std::vector<CCSeparatingPlane<T>*>& planes) {
    CCPairs.clear();
    planes.clear();
    for(auto& pair:_CCPlanes) {
      CCPairs.push_back(pair.first);
      planes.push_back(&(pair.second));
    }
  }
  //void parity check
  void parityCheck() const {
    for(auto it:_idToCCPairs)
      for(auto pair:it.second)
        if(_CCPlanes.find(pair)==_CCPlanes.end()) {
          pair.first.print();
          pair.second.print();
          ASSERT_MSG(false,"_BVHOffsetToCCPairs -> _CCPlanes inconsistent!")
        }
    for(auto it:_CCPlanes) {
      if(_idToCCPairs.find(it.first.first)==_idToCCPairs.end() ||
          _idToCCPairs.find(it.first.first)->second.find(it.first)==_idToCCPairs.find(it.first.first)->second.end()) {
        it.first.first.print();
        it.first.second.print();
        ASSERT_MSG(false,"_CCPlanes -> _BVHOffsetToCCPairs inconsistent!")
      }
      if(it.first.second.isRobotConvexHull())
        if(_idToCCPairs.find(it.first.second)==_idToCCPairs.end() ||
            _idToCCPairs.find(it.first.second)->second.find(it.first)==_idToCCPairs.find(it.first.second)->second.end()) {
          it.first.first.print();
          it.first.second.print();
          ASSERT_MSG(false,"_CCPlanes -> _BVHOffsetToCCPairs inconsistent!")
        }
    }
  }
  //debugger
  static inline void compareCCPair(const CCSeparatingPlanes<T>& CCPlanes,
                                   const CCSeparatingPlanes<T>& CCPlanesBF) {
    std::cout << "Found " << CCPlanes._CCPlanes.size() << " CCPlanes " << CCPlanesBF._CCPlanes.size() << " CCPlanesBF!" << std::endl;
    for(const auto& CC:CCPlanes._CCPlanes)
      if(CCPlanesBF._CCPlanes.find(CC.first)==CCPlanesBF._CCPlanes.end()) {
        CC.first.first.print();
        CC.first.second.print();
        ASSERT_MSG(false,"CCPlane not found in CCPlanesBF!")
      }
    for(const auto& CC:CCPlanesBF._CCPlanes)
      if(CCPlanes._CCPlanes.find(CC.first)==CCPlanes._CCPlanes.end()) {
        CC.first.first.print();
        CC.first.second.print();
        ASSERT_MSG(false,"CCPlaneBF not found in CCPlanes!")
      }
  }
 private:
  std::unordered_map<EntityId<T>,std::unordered_set<std::pair<EntityId<T>,EntityId<T>>,EntityIdPairHash<T>>,EntityIdHash<T>> _idToCCPairs;
  std::unordered_map<std::pair<EntityId<T>,EntityId<T>>,CCSeparatingPlane<T>,EntityIdPairHash<T>> _CCPlanes;
  bool _CCUpdated;
};
//debug
template <typename T>
static inline void compareTTPair(const std::vector<std::pair<EntityId<T>,EntityId<T>>>& TTPairs,
                                 const std::vector<std::pair<EntityId<T>,EntityId<T>>>& TTPairsBF) {
  std::cout << "Found " << TTPairs.size() << " TTPairs " << TTPairsBF.size() << " TTPairsBF!" << std::endl;
  std::set<std::pair<EntityId<T>,EntityId<T>>,LSS<T>> TTPairsSet(TTPairs.begin(),TTPairs.end());
  std::set<std::pair<EntityId<T>,EntityId<T>>,LSS<T>> TTPairsSetBF(TTPairsBF.begin(),TTPairsBF.end());
  ASSERT_MSG(TTPairs.size()==TTPairsSet.size(),"TTPairs not Unique!")
  ASSERT_MSG(TTPairsBF.size()==TTPairsSetBF.size(),"TTPairsBF not Unique!")
  for(const auto& TTPair:TTPairsSet)
    if(TTPairsSetBF.find(TTPair)==TTPairsSetBF.end()) {
      TTPair.first.print();
      TTPair.second.print();
      ASSERT_MSG(false,"TTPair not found in TTPairsBF!")
    }
  for(const auto& TTPair:TTPairsSetBF)
    if(TTPairsSet.find(TTPair)==TTPairsSet.end()) {
      TTPair.first.print();
      TTPair.second.print();
      ASSERT_MSG(false,"TTPairBF not found in TTPairs!")
    }
}
}
#endif
