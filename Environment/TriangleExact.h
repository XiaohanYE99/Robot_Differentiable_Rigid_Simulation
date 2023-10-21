#ifndef TRIANGLE_EXACT_H
#define TRIANGLE_EXACT_H

#include "BBoxExact.h"

namespace PHYSICSMOTION {
struct TriangleExact : public SerializableBase {
 public:
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  TriangleExact();
  TriangleExact(const Vec3T& a,const Vec3T& b,const Vec3T& c);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  const Vec3T& v(int i) const;
  const Vec3T& normal() const;
  BBoxExact getBB() const;
  Vec3T bary(const Vec3T& p) const;
  Vec3T interp(const Vec3T& bary) const;
  Vec2T project(const Vec3T& d) const;
  bool intersect(const BBoxExact& bb) const;
  T dirGrad(const Eigen::Matrix<int,2,1>& vid,const Vec3T& D) const;
  std::pair<int,Vec3T> moveFromEdge(const Eigen::Matrix<int,2,1>& evid,const Vec3T& p0,Vec3T D) const;
  void calcPointDist(const Vec3T& pt,T& sqrDistance,Vec3T& cp,Vec3T& b,Eigen::Matrix<int,2,1>& feat) const;
  Eigen::Matrix<int,3,1> _vid,_eNId;
 private:
  void reset();
  Vec3T _v[3],_n[3],_nO[3],_nOSqr,_nt;
  Mat3T _nTnN;
  Mat2T _invM;
};
}

#endif
