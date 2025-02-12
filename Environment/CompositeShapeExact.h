#ifndef COMPOSITE_SHAPE_EXACT_H
#define COMPOSITE_SHAPE_EXACT_H

#include "ShapeExact.h"
#include "BBoxExact.h"
#include "MeshExact.h"

struct aiNode;
struct aiScene;
struct aiMaterial;
struct aiMaterialProperty;
namespace PHYSICSMOTION {
struct CompositeShapeExact : public ShapeExact {
  struct Material : public SerializableBase {
    Material();
    virtual bool read(std::istream& is,IOData* dat) override;
    virtual bool write(std::ostream& os,IOData* dat) const override;
    virtual std::shared_ptr<SerializableBase> copy() const override;
    virtual std::string type() const override;
    Eigen::Matrix<double,4,1> _ambient,_diffuse,_specular;
    double _shininess;
    bool _useWireframe;
    std::string _texFile;
  };
  CompositeShapeExact();
  CompositeShapeExact(const std::string& path,bool buildBVH=true);
  CompositeShapeExact(const std::vector<std::shared_ptr<ShapeExact>>& geoms,
                      const std::vector<Mat3X4T>& trans);
  CompositeShapeExact(const std::vector<std::shared_ptr<ShapeExact>>& geoms);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual const BBoxExact& getBB() const override;
  virtual bool empty() const override;
  virtual void getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                       std::vector<Eigen::Matrix<int,3,1>>& iss) const override;
  virtual bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                            T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                            std::vector<Vec3T>* history=NULL) const override;
  void transform(const Mat3X4T& trans);
  void scale(T coef) override;
  const std::vector<std::shared_ptr<ShapeExact>>& getGeoms() const;
  const std::vector<Mat3X4T>& getTrans() const;
  const std::vector<Material>& getMaterials() const;
  std::vector<Material>& getMaterials();
  void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const override;
 protected:
  void initBB();
  void init(const aiScene* scene,const aiNode* node,bool buildBVH,Mat3X4T trans,const std::string& path);
  static Material initMaterial(const aiMaterial& mat,const std::string& path);
  static Eigen::Matrix<double,4,1> parseVec4(const aiMaterialProperty& matP);
  static double parseFloat(const aiMaterialProperty& matP,int index=0);
  static int parseInt(const aiMaterialProperty& matP,int index=0);
  std::vector<std::shared_ptr<ShapeExact>> _geoms;
  std::vector<Mat3X4T> _trans;
  std::vector<Material> _materials;
  BBoxExact _bb;
};
}

#endif
