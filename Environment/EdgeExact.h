#ifndef EDGE_EXACT_H
#define EDGE_EXACT_H

#include "BBoxExact.h"

namespace PHYSICSMOTION {
class EdgeExact : public SerializableBase {
 public:
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  EdgeExact();
  EdgeExact(const Vec3T& a,const Vec3T& b);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  std::pair<int,Vec3T> moveFromVertex(const Vec3T& d0,const Vec3T& D) const;
  T dirGrad(int vid,const Vec3T& D) const;
  Eigen::Matrix<int,2,1> _vid,_tNId;
 private:
  void reset();
  Vec3T _v[2],_d;
  Mat3T _dTdN;
  T _lenSqr;
};
}

#endif
