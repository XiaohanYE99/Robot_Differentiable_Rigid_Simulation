#ifndef PENETRATION_DETECTOR_H
#define PENETRATION_DETECTOR_H

#include "EntityId.h"
#include "ThetaTrajectory.h"
#include <Environment/MeshExact.h>
#include <Environment/PointCloudExact.h>
#include <Articulated/ArticulatedBody.h>

namespace PHYSICSMOTION {
template <typename T>
struct PDEntrySpan {
  DECL_MAT_VEC_MAP_TYPES_T
  typedef CollisionGradInfo<T>
  GradInfo;
  //initialize
  PDEntrySpan();
  T getTimeAvg() const;
  T getTimeDiff() const;
  bool operator<(const PDEntrySpan<T>& other) const;
  bool equal(const PDEntrySpan<T>& other) const;
  void writeVTK(VTKWriter<double>& os,const GradInfo& info) const;
  void print() const;
  //data
  int _jidA,_jidB;
  Vec3T _pAL,_pBL;
  T _PD;
  //time span
  T _L1A,_L1B;
  T _timeFrom;
  T _timeTo;
};
template <typename T>
class PenetrationDetector : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef CollisionGradInfo<T>
  GradInfo;
  //initialize
  explicit PenetrationDetector(T d0=1e-2,T epsTime=1e-3,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD= {});
  explicit PenetrationDetector(std::shared_ptr<ArticulatedBody> body,int order,int totalTime,
                               T d0=1e-2,T epsTime=1e-3,T randomInitialize=1.,bool neutralInit=false,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD= {});
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  void setRobot(std::shared_ptr<ArticulatedBody> body,bool randomInitialize,bool neutralInit);
  void addRobot(std::shared_ptr<ArticulatedBody> body,const Vec& DOF);
  void addMesh(std::shared_ptr<ShapeExact> shape,const Mat3X4T& trans);
  void addObject(const std::string& path,T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHull=10);
  void addSphere(T r,int res,const Vec3T& pos);
  void addCuboid(T l,T w,T h,const Vec3T& pos,const Mat3T& R);
  void addCapsule(const Vec3T& a,const Vec3T& b,T radius,int res);
  void addCapsule(T l,T w,T h,Vec3T pos,T radius,int res);
  std::vector<std::shared_ptr<PointCloudExact>> getBodyPoints() const;
  std::vector<std::shared_ptr<MeshExact>> getObstacles() const;
  std::vector<std::shared_ptr<PointCloudExact>> getObstaclePoints() const;
  std::shared_ptr<ArticulatedBody> getBody() const;
  void assembleBodyBVH();
  void assembleObsBVH();
  void reset();
  //detect penetration
  bool DCD(bool clear,bool BF,bool output=false);
  PDEntrySpan<T> DCD(int jid,int jid2,bool BF);
  PDEntrySpan<T> queryPD(T timeFrom,T timeTo,int jid,int jid2);
  PDEntrySpan<T> queryPDSelf(T timeFrom,T timeTo,int jid,int jid2);
  PDEntrySpan<T> queryPDObs(T timeFrom,T timeTo,int jid,int jid2);
  PDEntrySpan<T> queryPDObs(T timeFrom,T timeTo,int jid);
  void print() const;
  void writeVTK(const std::string& path);
  void debugDCDBF();
  //info query
  ParallelVector<PDEntrySpan<T>>& getPDEntries();
  const ParallelVector<PDEntrySpan<T>>& getPDEntries() const;
  const ThetaTrajectory<T>& getThetaTrajectory() const;
  const Vec& getControlPoints() const;
  void setControlPoints(const Vec& x);
  void setDCDObs(bool DCDObs);
  void setDCDSelf(bool DCDSelf);
  void setUseSDF(bool useSDF);
  bool getDCDObs() const;
  bool getDCDSelf() const;
  bool getUseSDF() const;
  T d0() const;
  T epsTime() const;
  //get L1
  void initL1();
  void debugL1(int res);
  T getConvexHullL1(int jointId,T t0,T t1) const;
  T getVertexL1(int jointId,int j,const Vec& maxGrad) const;
  T getVertexL1(int jointId,int j,T t0,T t1) const;
  const GradInfo& getInfo(T timeAvg,bool mustExist=false) const;
  bool excludedFromDCDSelf(int jid1,int jid2) const;
 protected:
  //data
  std::shared_ptr<ArticulatedBody> _body;
  std::vector<std::shared_ptr<PointCloudExact>> _bodyPoints;
  std::vector<std::shared_ptr<MeshExact>> _obstacles;
  std::vector<std::shared_ptr<PointCloudExact>> _obstaclePoints;
  std::vector<Node<int,BBoxExact>> _bvhObstacle;
  //parameters
  T _d0,_epsTime;
  bool _DCDObs,_DCDSelf,_useSDF;
  Vec _controlPoints;
  ThetaTrajectory<T> _thetaTrajectory;
  std::unordered_map<int,std::unordered_set<int>> _skipJIDPairs;
  //temporary data, not serialized
  ParallelVector<PDEntrySpan<T>> _PDEntries;
  std::unordered_map<T,GradInfo> _infoLookup;
};
}
#endif
