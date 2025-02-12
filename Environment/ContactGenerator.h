#ifndef CONTACT_GENERATOR_H
#define CONTACT_GENERATOR_H

#include "BVHNode.h"
#include "ShapeExact.h"
#include "Articulated/ArticulatedBody.h"

namespace PHYSICSMOTION {
class ContactGenerator {
 public:
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  enum STATUS {
    STATIC_STATIC=1,
    STATIC_DYNAMIC=2,
    DYNAMIC_DYNAMIC=4,
  };
  struct ContactPoint {
    ContactPoint();
    T depth() const;
    void swap();
    void transform(const Mat3X4T& t);
    Vec3T _ptA,_ptB,_nA2B;
    //simulator parameter
    T _tangentBound;
    Vec3T _ptALast,_ptBLast;
    Vec3T _ptAL,_ptBL;
    Vec3T _fA,_fB;
    Mat3X2T _tA2B;
  };
  struct ContactManifold {
    ContactManifold();
    void swap();
    bool operator<(const ContactManifold& other) const;
    std::vector<ContactPoint> _points;
    std::shared_ptr<ShapeExact> _sA,_sB;
    int _sidA,_sidB,_jidA,_jidB;
    Mat3X4T _tA,_tB;
    std::vector<Mat3X4T> _DNDX;
    Vec4T _x;
  };
  typedef ShapeExact::Facet Facet;
  typedef ShapeExact::Edge Edge;
  ContactGenerator(std::shared_ptr<ArticulatedBody> body,std::vector<std::shared_ptr<ShapeExact>> shapes,bool buildBVH=true);
  void generateManifolds(T x0,bool useCCD,std::vector<ContactManifold>& manifolds,Mat3XT t,int status=STATIC_DYNAMIC|DYNAMIC_DYNAMIC);
  const std::unordered_set<Eigen::Matrix<int,2,1>,EdgeHash>& getExclude() const;
  void updateBVH(Mat3XT& t,T x0=0.0);
  BBoxExact getBB() const;
  static T epsDist();
  static T epsDir();
 private:
  void generateManifold(std::vector<ContactManifold>& manifolds,ContactManifold m);
  static void generateManifoldSphereSphereInternal(std::vector<ContactManifold>& manifolds,ContactManifold& m,const Vec3T& cA1,const Vec3T& cA2,const Vec3T& cB1,const Vec3T& cB2);
  static bool generateManifoldSphereSphere(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  static void generateManifoldSphereCapsuleInternal(std::vector<ContactManifold>& manifolds,ContactManifold& m,const Vec3T& cA,const Vec3T& cB1,const Vec3T& cB2);
  static bool generateManifoldSphereCapsule(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  static bool generateManifoldCapsuleCapsule(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  static bool generateManifoldSphereBox(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  bool generateManifoldCapsuleBox(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  bool generateManifoldBoxBox(std::vector<ContactManifold>& manifolds,ContactManifold& m);
  //data
  std::shared_ptr<ArticulatedBody> _body;
  std::vector<std::shared_ptr<ShapeExact>> _staticShapes;
  std::vector<Node<int,BBoxExact>> _staticBVH,_dynamicBVH;
  std::unordered_map<std::shared_ptr<ShapeExact>,std::vector<Facet>> _facetCache;
  std::unordered_map<std::shared_ptr<ShapeExact>,std::vector<Edge>> _edgeCache;
  std::unordered_set<Eigen::Matrix<int,2,1>,EdgeHash> _exclude;
  static T _epsDir,_epsDist;
  bool _buildBVH;
};
}

#endif
