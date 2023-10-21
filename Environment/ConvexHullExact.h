#ifndef CONVEX_HULL_EXACT_H
#define CONVEX_HULL_EXACT_H

#include "EdgeExact.h"
#include "MeshExact.h"

namespace PHYSICSMOTION {
struct ConvexHullExact : public MeshExact {
  ConvexHullExact();
  ConvexHullExact(const std::string& path);
  ConvexHullExact(const aiScene* scene);
  ConvexHullExact(const MeshExact& m);
  ConvexHullExact(const std::vector<Eigen::Matrix<double,3,1>>& vss);
  void init(const aiScene* scene,const aiNode* node,
            std::vector<Eigen::Matrix<double,3,1>>& vss);
  template <typename T2>
  void init(const std::vector<Eigen::Matrix<T2,3,1>>& vss);
  template <typename T2>
  void init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
            const std::vector<Eigen::Matrix<int,3,1>>& iss);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  void parityCheck() const;
  bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                    T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                    std::vector<Vec3T>* history=NULL) const override;
  void scale(T coef) override;
  Eigen::Matrix<T,4,1> plane(int i) const;
  int nrPlane() const;
  //for GJK
  Vec3T support(const Vec3T& D,int& id) const override;
  //for SAT
  std::vector<Facet> facets() const override;
  std::vector<Edge> edges() const override;
#ifndef SWIG
  void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const override;
#endif
 private:
  //data
  std::vector<EdgeExact> _ess;
  std::vector<std::vector<int>> _eNss;
  static T _epsEdge;
};
}

#endif
