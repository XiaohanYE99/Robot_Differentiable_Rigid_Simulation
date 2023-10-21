#ifndef SPHERICAL_BBOX_EXACT_H
#define SPHERICAL_BBOX_EXACT_H

#include "BBoxExact.h"

namespace PHYSICSMOTION {
struct SphericalBBoxExact : public BBoxExact {
  SphericalBBoxExact();
  SphericalBBoxExact(T rad);
  SphericalBBoxExact(T a,T rad);
  SphericalBBoxExact(T a,T b,T rad);
  SphericalBBoxExact(T a,T b,T c,T rad);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  void getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
               std::vector<Eigen::Matrix<int,3,1>>& iss) const override;
  bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                    T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                    std::vector<Vec3T>* history=NULL) const override;
  void scale(T coef) override;
  virtual bool contain(const Vec3T& pt) const;
  bool isRoundCube() const;
  bool isRoundBoard() const;
  bool isCapsule() const;
  bool isSphere() const;
  T radius() const;
  //for SAT
  std::vector<Facet> facets() const override;
  std::vector<Edge> edges() const override;
 private:
  BBoxExact _bb;
  T _rad,_radSqr;
};
}

#endif
