#ifndef MESH_EXACT_H
#define MESH_EXACT_H

#include "ShapeExact.h"
#include "TriangleExact.h"

struct aiNode;
struct aiScene;
namespace PHYSICSMOTION {
struct MeshExact : public ShapeExact {
  MeshExact();
  MeshExact(const std::string& path,bool buildBVH=true);
  MeshExact(const aiScene* scene,bool buildBVH=true);
  MeshExact(std::vector<Eigen::Matrix<double,3,1>>& vss,
            std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH=true,bool buildBSH=true);
  MeshExact(std::vector<Eigen::Matrix<double,3,1>>& vss,
            std::vector<Eigen::Matrix<double,2,1>>& tcss,
            std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH=true);
  void init(const aiScene* scene,const aiNode* node,
            std::vector<Eigen::Matrix<double,3,1>>& vss,
            std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH=true);
  template <typename T2>
  void init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
            const std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH=true,bool buildBSH=true);
  template <typename T2>
  void init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
            const std::vector<Eigen::Matrix<T2,2,1>>& tcss,
            const std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH=true);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  const std::vector<Node<int,BBoxExact>>& getBVH() const;
  const std::vector<char>& bss()const;
  const std::vector<Vec3T>& vss() const;
  std::vector<Vec3T>& vssNonConst();
  const std::vector<Vec2T>& tcss() const;
  const std::vector<Eigen::Matrix<int,3,1>>& iss() const;
  void getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
               std::vector<Eigen::Matrix<int,3,1>>& iss) const override;
  bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                    T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                    std::vector<Vec3T>* history=NULL) const override;
  void scale(T coef) override;
  void translate(const Vec3T& delta);
  void transform(const Mat3X4T& trans);
  void moveMesh(const Vec& delta);
  void setMesh(const Vec& X);
  MeshExact addMesh(const MeshExact& mesh);
  //for GJK
  Vec3T support(const Vec3T& D,int& id) const override;
  //for SAT
  Vec2T project(const Vec3T& d) const override;
  void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const override;
 protected:
  void initBss();
  std::vector<TriangleExact> _tss;
  std::vector<Vec3T> _vss;
  std::vector<Vec2T> _tcss;
  std::vector<Eigen::Matrix<int,3,1>> _iss;
  // belongingness set
  std::vector<char> _bss;
  std::vector<Node<int,BBoxExact>> _bvh;
  std::vector<Node<int,BBoxExact>> _bsh;
};
}

#endif
