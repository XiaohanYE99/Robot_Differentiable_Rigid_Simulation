#ifndef CONVEX_DECOMPOSITION_H
#define CONVEX_DECOMPOSITION_H

#include "MeshExact.h"
#include "Articulated/ArticulatedBody.h"

namespace PHYSICSMOTION {
class ConvexDecomposition : public SerializableBase {
 public:
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  explicit ConvexDecomposition();
  explicit ConvexDecomposition(const MeshExact& mesh,T scale=1.,const Vec3T& pos=Vec3T(0,0,0),const Mat3T& rot=Mat3T::Identity(),int maxConvexHulls=4);
  explicit ConvexDecomposition(const std::string &inputFile,T scale=1.,const Vec3T& pos=Vec3T(0,0,0),const Mat3T& rot=Mat3T::Identity(),int maxConvexHulls=4);
  void convexDecomposition(T scale=1.,const Vec3T& pos=Vec3T(0,0,0),const Mat3T& rot=Mat3T::Identity(),int maxConvexHulls=4);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  int getNumConvexHull();
  const std::vector<std::shared_ptr<MeshExact>>& getConvexHulls() const;
 private:
  void decompose(const MeshExact& mesh,T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHulls);
  std::vector<float> _points;
  std::vector<uint32_t> _triangles;
  std::vector<std::shared_ptr<MeshExact>> _convexHulls;
};
}

#endif
