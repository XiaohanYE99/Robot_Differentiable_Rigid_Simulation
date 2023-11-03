#ifndef BBOX_EXACT_H
#define BBOX_EXACT_H

#include "ShapeExact.h"

namespace PHYSICSMOTION {
struct BBoxExact : public ShapeExact {
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  BBoxExact();
  BBoxExact(T halfX,T halfY,T halfZ);
  BBoxExact(const Vec3T& a,const Vec3T& b);
  BBoxExact(const Vec3T& a,const Vec3T& b,const Vec3T& c);
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
  void setUnion(const BBoxExact& other);
  void setUnion(const Vec3T& other);
  void extendUnion(const T x0);
  const Vec3T& minCorner() const;
  const Vec3T& maxCorner() const;
  Vec3T& minCorner();
  Vec3T& maxCorner();
  bool intersect(const BBoxExact& other) const;
  virtual bool contain(const Vec3T& pt) const;
  BBoxExact enlargedEps(T eps) const;
  BBoxExact enlarged(const Vec3T& ext) const;
  T distToSqr(const BBoxExact& other) const;
  T distToSqr(const Vec3T& pt) const;
  bool operator==(BBoxExact& bbox);
  //for GJK
  Vec3T support(const Vec3T& D,int& id) const override;
  //for SAT
  Vec2T project(const Vec3T& d) const override;
  std::vector<Facet> facets() const override;
  std::vector<Edge> edges() const override;
#ifndef SWIG
  void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const override;
#endif
 protected:
  static void makeGrid(std::vector<Eigen::Matrix<double,3,1>>& vss,std::vector<Eigen::Matrix<int,3,1>>& iss,
                       int RESX,int RESY,const Eigen::Matrix<double,3,1>& ctr,const Eigen::Matrix<double,3,1>& d0,const Eigen::Matrix<double,3,1>& d1);
  static void makeBox(std::vector<Eigen::Matrix<double,3,1>>& vss,std::vector<Eigen::Matrix<int,3,1>>& iss,
                      int RES,const Eigen::Matrix<double,3,1>& halfExt);
  static int normalToVid(const Vec3T& normal);
  static int normalToEid0(const Vec3T& normal);
  static int normalToEid1(const Vec3T& normal);
  static int normalToFid(const Vec3T& normal);
  Vec3T vertex(int id) const;
  Vec3T _minC,_maxC;
};
}

#endif
