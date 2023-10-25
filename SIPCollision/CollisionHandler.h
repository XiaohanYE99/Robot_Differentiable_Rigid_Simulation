#ifndef COLLISION_HANDLER_H
#define COLLISION_HANDLER_H

#include "EntityId.h"
#include "ThetaTrajectory.h"
#include "CCSeparatingPlanes.h"
#include <Environment/MeshExact.h>
#include <Articulated/ArticulatedBody.h>

namespace PHYSICSMOTION {
enum PairStatus {
  PS_Safe,
  PS_Unsafe,
  PS_Penetrated,
  PS_LeftLargerInterval,
  PS_RightLargerInterval,
};
enum SubdivideType {
  ST_All,
  ST_Eval,
  ST_Safety,
};
template <typename T>
class GJKPolytope;
template <typename T,typename PFunc>
class CollisionBarrierEnergy;
template <typename T>
class CollisionHandler : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef CollisionGradInfo<T> GradInfo;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  //initialize
  static Vec initializeControlPoints
  (T randomInitialize,bool neutralInit,
   std::shared_ptr<ArticulatedBody> body,
   ThetaTrajectory<T>& thetaTrajectory);
  explicit CollisionHandler(T d0=1e-3,T x0=1.,T l2=0.01,T eta=7./24.,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD= {});
  explicit CollisionHandler(std::shared_ptr<ArticulatedBody> body,int order,int totalTime,
                            T d0=1e-3,T x0=1.,T l2=0.01,T eta=7./24.,T randomInitialize=1.,bool neutralInit=false,const std::unordered_map<int,std::unordered_set<int>>& skipSelfCCD= {});
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
  const std::vector<Node<EntityId<T>,BBoxExact>>& getBodyBVH() const;
  std::vector<std::shared_ptr<MeshExact>> getObstacles() const;
  std::shared_ptr<ArticulatedBody> getBody() const;
  void assembleBodyBVH(bool useConvexHull=false);
  void assembleObsBVH();
  void reset();
  //return success or fail: failedOffset is undefined
  //when fail, failedOffset return the leaf node into _bvhBody, which fails the safety-check
  PairStatus singlePairCCDObs(int bvhBodyOffset,int obstacleId,int bvhObstacleOffset,bool doSubdivide,bool buildPair);
  PairStatus singlePairCCDSelf(int bvhBodyOffset1,int bvhBodyOffset2,bool doSubdivide,bool buildPair,bool eval);
  void narrowphaseCCD(int bvhBodyOffset,int obstacleId,bool doSubdivide,bool buildPair);
  bool CCD(bool doSubdivide,bool buildPair,bool eval);
  bool CCDObs(bool doSubdivide,bool buildPair,bool clear);
  bool CCDSelf(bool doSubdivision, bool buildPair,bool clear,bool eval);
  bool CCDBF(bool doSubdivide,bool buildPair,bool eval);
  bool CCDObsBF(bool doSubdivide,bool buildPair,bool clear);
  bool CCDSelfBF(bool doSubdivision, bool buildPair,bool clear,bool eval);
  void debugCCDObsBF(bool doSubdivide=false,bool buildPair=false,bool exhaustive=false);
  void debugCCDSelfBF(bool doSubdivide=false,bool buildPair=false,bool exhaustive=false,bool eval=false);
  //we will input the offset into the vector _bvhBody, and we guarantee _bvhBody[offset] is a leaf node
  //we modify the _bvhBody vector, subdividing the leaf into two
  bool exhaustiveSubdivide(bool bruteForce,bool eval,int iters=-1);
  void subdivideSingle(int offset,bool assertNoReserve);
  void subdivide(SubdivideType type=ST_All,bool debug=false);
  std::string findEntityId(const EntityId<T>& id) const;
  template <typename PFunc>
  void traceError(const CollisionBarrierEnergy<T,PFunc>& energy) const {
    for(auto TT:energy.getFailedTT())
      std::cout << "Failed TTPair: " << findEntityId(TT.first) << "," << findEntityId(TT.second) << std::endl;
    for(auto CC:energy.getFailedCC())
      std::cout << "Failed CCPair: " << findEntityId(CC.first) << "," << findEntityId(CC.second) << std::endl;
  }
  void parityCheck() const;
  void clearCC();
  void update();
  //info query
  const std::unordered_map<int,std::unordered_set<int>>& getSkipJIDPairs() const;
  const std::unordered_set<int>& getSubdivideOffsets() const;
  const std::vector<Node<EntityId<T>,BBoxExact>>& getBVHBody() const;
  const std::vector<std::pair<EntityId<T>,EntityId<T>>>& getObsTTPairs() const;
  std::vector<std::pair<EntityId<T>,EntityId<T>>>& getObsTTPairs();
  const CCSeparatingPlanes<T>& getObsCCPlanes() const;
  CCSeparatingPlanes<T>& getObsCCPlanes();
  const std::vector<std::pair<EntityId<T>,EntityId<T>>>& getSelfTTPairs() const;
  std::vector<std::pair<EntityId<T>,EntityId<T>>>& getSelfTTPairs();
  const CCSeparatingPlanes<T>& getSelfCCPlanes() const;
  CCSeparatingPlanes<T>& getSelfCCPlanes();
  const std::unordered_map<rational,GradInfo>& getInfoLookup() const;
  std::unordered_map<rational,GradInfo>& getInfoLookup();
  const ThetaTrajectory<T>& getThetaTrajectory() const;
  const Vec& getControlPoints() const;
  void setControlPoints(const Vec& x);
  void setCCDObs(bool CCDObs);
  void setCCDSelf(bool CCDSelf);
  T d0() const;
  T d0Sqr() const;
  T l2() const;
  T eta() const;
  T x0() const;
  //get L1
  void initL1();
  void debugL1(int res);
  T getConvexHullL1(int jointId,T t0,T t1) const;
  T getVertexL1(int jointId,int j,const Vec& maxGrad) const;
  T getVertexL1(int jointId,int j,T t0,T t1) const;
  T getFaceL1(int jointId,int tid,T t0,T t1) const;
  T getEdgeL1(int jointId,int tid,int i,T t0,T t1) const;
  //get Phi
  T getConvexHullPhi(int bvhBodyOffset) const;
  T getVertexPhi(int bvhBodyOffset,int vertexId) const;
  T getEdgePhi(int bvhBodyOffset,int startVertexId) const;
  T getFacePhi(int bvhBodyOffset) const;
  bool penetrated() const;
 protected:
  BBoxExact computeBB(const EntityId<T>& id,T t0,T t1);
  const GradInfo& getInfo(rational timeAvg,bool mustExist=false) const;
  bool excludedFromCCDSelf(int jid1,int jid2) const;
  //data
  std::shared_ptr<ArticulatedBody> _body;
  std::vector<Node<EntityId<T>,BBoxExact>> _bvhBody;
  std::vector<std::shared_ptr<MeshExact>> _obstacles;
  std::vector<std::shared_ptr<GJKPolytope<T>>> _obstaclePolytopes;
  std::vector<Node<int,BBoxExact>> _bvhObstacle;
  //parameters
  T _d0,_d0Sqr,_l2,_eta,_x0;
  bool _CCDObs,_CCDSelf;
  Vec _controlPoints;
  ThetaTrajectory<T> _thetaTrajectory;
  std::unordered_map<int,std::unordered_set<int>> _skipJIDPairs;
  //temporary data, not serialized
  CCSeparatingPlanes<T> _obsCCPlanes;
  CCSeparatingPlanes<T> _selfCCPlanes;
  std::unordered_map<rational,GradInfo> _infoLookup;
  std::unordered_set<int> _penetratedOffsets;
  std::unordered_set<int> _subdivideOffsets;
  std::unordered_set<int> _evaluateSubdivisionOffsets;
  std::vector<std::pair<EntityId<T>,EntityId<T>>> _obsTTPairs;
  std::vector<std::pair<EntityId<T>,EntityId<T>>> _selfTTPairs;
};
}
#endif
