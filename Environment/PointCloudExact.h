#ifndef POINT_CLOUD_EXACT_H
#define POINT_CLOUD_EXACT_H

#include "BBoxExact.h"
#include "MeshExact.h"

namespace PHYSICSMOTION {
struct PointCloudExact : public ShapeExact {
  PointCloudExact();
  PointCloudExact(const MeshExact& mesh,T d0);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  virtual void getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                       std::vector<Eigen::Matrix<int,3,1>>& iss) const;
  virtual bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                            T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                            std::vector<Vec3T>* history=NULL) const;
  template <typename T2>
  T2 closestSDF(const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& normal) const {
    Eigen::Matrix<double,3,1> ptR=pt.template cast<double>(),normalR;
    double dR=closestSDFInner(ptR,normalR);
    normal=normalR.template cast<T2>();
    return (T2)dR;
  }
#ifndef SWIG
  const std::vector<Node<int,BBoxExact>>& getBVH() const;
#endif
  const std::vector<Vec3T>& vss() const;
  void scale(T coef) override;
  void translate(const Vec3T& delta);
  void transform(const Mat3X4T& trans);
  //for GJK
  virtual Vec3T support(const Vec3T& D,int& id) const override;
#ifndef SWIG
  virtual void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const;
  void writeSDFVTK(const std::string& path) const;
#endif
 protected:
  void init(bool buildBVH);
  void initSDF(const MeshExact& mesh,T d0,int extend=5);
  double closestSDFInner(const Eigen::Matrix<double,3,1>& pt,Eigen::Matrix<double,3,1>& normal) const;
  //data
  std::vector<Vec3T> _vss;
  std::vector<Node<int,BBoxExact>> _bvh;
  //distance field
  Eigen::Matrix<int,3,1> _nrNode,_stride;
  Eigen::Matrix<double,3,1> _bbSDFMinC;
  Eigen::Matrix<double,3,1> _bbSDFMaxC;
  std::vector<double> _SDF;
};
}

#endif
