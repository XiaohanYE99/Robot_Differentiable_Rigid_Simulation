#ifndef SAT_H
#define SAT_H

#include "ShapeExact.h"
#include "ContactGenerator.h"

namespace PHYSICSMOTION {
struct SAT {
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef ContactGenerator::ContactPoint ContactPoint;
  typedef ContactGenerator::ContactManifold ContactManifold;
  typedef std::vector<Vec2T> Facet2D;
  typedef ShapeExact::Facet Facet;
  typedef ShapeExact::Edge Edge;
  struct ProjRange {
    int _fidA,_eidA;
    int _fidB,_eidB;
    Vec2T _rngA,_rngB;
    Vec3T _n;
    T _depth;
    //additional data for edge-edge check
    Vec3T _ptA,_ptB;
  };
  static T depth(const ProjRange& rng);
  static bool matchWitness(const ProjRange& rng,T witnessA,T* verboseB);
  static bool matchWitness(const ProjRange& rng,T* verboseA,T witnessB);
  static bool matchWitness(const ProjRange& rng,T witnessA,T witnessB);
  static Vec2T project(std::shared_ptr<ShapeExact> S,
                       const Mat3X4T& trans,
                       const Vec3T& n);
  static ProjRange runSAT(std::shared_ptr<ShapeExact> A,
                          std::shared_ptr<ShapeExact> B,
                          const std::vector<Facet>& FA,
                          const std::vector<Facet>& FB,
                          const std::vector<Edge>& EA,
                          const std::vector<Edge>& EB,
                          const Mat3X4T& transA,
                          const Mat3X4T& transB,bool* intersect=NULL);
  static ProjRange runSAT(std::shared_ptr<ShapeExact> A,
                          std::shared_ptr<ShapeExact> B,
                          const Mat3X4T& transA,
                          const Mat3X4T& transB,bool* intersect=NULL);
  //not only run contact, but also generate manifold
  static ProjRange generateManifold(std::shared_ptr<ShapeExact> A,
                                    std::shared_ptr<ShapeExact> B,
                                    const std::vector<Facet>& FA,
                                    const std::vector<Facet>& FB,
                                    const std::vector<Edge>& EA,
                                    const std::vector<Edge>& EB,
                                    const Mat3X4T& transA,
                                    const Mat3X4T& transB,
                                    ContactManifold& m,bool* intersect=NULL);
  static ProjRange generateManifold(std::shared_ptr<ShapeExact> A,
                                    std::shared_ptr<ShapeExact> B,
                                    const Mat3X4T& transA,
                                    const Mat3X4T& transB,
                                    ContactManifold& m,bool* intersect=NULL);
  //clip a facet f against a reference
  static void generateManifoldFace(std::shared_ptr<ShapeExact> A,
                                   std::shared_ptr<ShapeExact> B,
                                   const std::vector<Facet>& FA,
                                   const std::vector<Facet>& FB,
                                   const std::vector<Edge>& EA,
                                   const std::vector<Edge>& EB,
                                   const Mat3X4T& transA,
                                   const Mat3X4T& transB,
                                   int fidA,std::vector<ContactPoint>& points);
  static void clip(Facet& f,const Vec3T& pos,const Vec3T& inward);
  static void clip(Facet& f,const Facet& ref);
};
}
#endif
