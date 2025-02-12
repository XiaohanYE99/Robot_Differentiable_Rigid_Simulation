#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "BBoxExact.h"
#include "MeshExact.h"

namespace PHYSICSMOTION {
struct EnvironmentCallback {
 public:
  virtual ~EnvironmentCallback() {}
  virtual double height(double x,double y) const {
    return 0;
  }
};
template <typename T>
class Environment : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  virtual ~Environment() {}
  virtual T phi(const Vec3T& x) const;
  virtual Vec3T phiGrad(const Vec3T& x) const;
  virtual T phi(const Vec3T& x,Vec3T* g) const=0;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h) const=0;
  virtual const std::vector<Node<int,BBoxExact>>& getBVH() const=0;
  virtual const BBoxExact& getBB() const=0;
  virtual bool empty() const=0;
  virtual void createHills(double xMin,double xMax,double yMin,double yMax,EnvironmentCallback& h,double res);
  virtual void createHills(double xMin,double xMax,double yMin,double yMax,std::function<double(double,double)> h,double res)=0;
  virtual void createStair(double x,double y,double x0,double z0,double slope,int n);
  virtual void createFloor(const Eigen::Matrix<double,4,1>& plane);
};
template <typename T>
class EnvironmentExact : public Environment<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  EnvironmentExact();
  EnvironmentExact(const MeshExact& obj);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual T phi(const Vec3T& x,Vec3T* g) const override;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h) const override;
  virtual const std::vector<Node<int,BBoxExact>>& getBVH() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  virtual void createHills(double xMin,double xMax,double yMin,double yMax,std::function<double(double,double)> h,double res) override;
  const MeshExact& getMesh() const;
 private:
  MeshExact _obj;
};
template <typename T>
class EnvironmentHeight : public Environment<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  EnvironmentHeight();
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual T phi(const Vec3T& x,Vec3T* g) const override;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h) const override;
  virtual const std::vector<Node<int,BBoxExact>>& getBVH() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  virtual void createHills(double xMin,double xMax,double yMin,double yMax,std::function<double(double,double)> h,double res) override;
  Eigen::Matrix<T,-1,-1> getHeightMatrix() const;
 private:
  void computeInvCellSize();
  T getHeight(const Eigen::Matrix<int,2,1>& id) const;
  Eigen::Matrix<int,2,1> getPosId(const Vec3T& pt,Vec2T& pos) const;
  int buildBVH(int xFrom,int xTo,int yFrom,int yTo);
  BBoxExact BBCell(int x,int y) const;
  void toBezier(T val[4][4]) const;
  void buildBVH();
  //data
  std::vector<T> _height;
  std::vector<Node<int,BBoxExact>> _bvh;
  T _invCellSzX,_invCellSzY;
  int _resX,_resY;
  BBoxExact _bb;
};
}

#endif
