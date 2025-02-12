#include "CompositeShapeExact.h"
#include "MeshExact.h"
#include "ConvexHullExact.h"
#include "SphericalBBoxExact.h"
#include <Utils/IO.h>
#include <Utils/Utils.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>

namespace PHYSICSMOTION {
//Material
CompositeShapeExact::Material::Material() {
  _ambient.setZero();
  _diffuse.setZero();
  _specular.setZero();
  _shininess=10;
  _useWireframe=false;
}
bool CompositeShapeExact::Material::read(std::istream& is,IOData*) {
  readBinaryData(_ambient,is);
  readBinaryData(_diffuse,is);
  readBinaryData(_specular,is);
  readBinaryData(_shininess,is);
  readBinaryData(_useWireframe,is);
  readBinaryData(_texFile,is);
  return is.good();
}
bool CompositeShapeExact::Material::write(std::ostream& os,IOData*) const {
  writeBinaryData(_ambient,os);
  writeBinaryData(_diffuse,os);
  writeBinaryData(_specular,os);
  writeBinaryData(_shininess,os);
  writeBinaryData(_useWireframe,os);
  writeBinaryData(_texFile,os);
  return os.good();
}
std::shared_ptr<SerializableBase> CompositeShapeExact::Material::copy() const {
  return std::shared_ptr<SerializableBase>(new Material);
}
std::string CompositeShapeExact::Material::type() const {
  return typeid(Material).name();
}
//CompositeShapeExact
CompositeShapeExact::CompositeShapeExact() {}
CompositeShapeExact::CompositeShapeExact(const std::string& path,bool buildBVH) {
  Assimp::Importer importer;
  const aiScene* scene=importer.ReadFile(path.c_str(),aiProcess_JoinIdenticalVertices);
  ASSERT_MSGV(scene,"Mesh %s is empty!",path.c_str())
  std::filesystem::path dir(path);
  init(scene,NULL,buildBVH,Mat3X4T::Identity(),dir.parent_path().string());
  initBB();
}
CompositeShapeExact::CompositeShapeExact
(const std::vector<std::shared_ptr<ShapeExact>>& geoms,
 const std::vector<Mat3X4T>& trans):_geoms(geoms),_trans(trans) {
  ASSERT(_trans.size()==_geoms.size())
  initBB();
}
CompositeShapeExact::CompositeShapeExact
(const std::vector<std::shared_ptr<ShapeExact>>& geoms):_geoms(geoms),_trans(geoms.size(),Mat3X4T::Identity()) {
  initBB();
  ASSERT(_trans.size()==_geoms.size())
}
bool CompositeShapeExact::read(std::istream& is,IOData* dat) {
  registerType<MeshExact>(dat);
  registerType<ConvexHullExact>(dat);
  registerType<BBoxExact>(dat);
  registerType<SphericalBBoxExact>(dat);
  readBinaryData(_geoms,is,dat);
  readBinaryData(_trans,is,dat);
  readBinaryData(_materials,is,dat);
  readBinaryData(_bb,is,dat);
  return is.good();
}
bool CompositeShapeExact::write(std::ostream& os,IOData* dat) const {
  registerType<MeshExact>(dat);
  registerType<ConvexHullExact>(dat);
  registerType<BBoxExact>(dat);
  registerType<SphericalBBoxExact>(dat);
  writeBinaryData(_geoms,os,dat);
  writeBinaryData(_trans,os,dat);
  writeBinaryData(_materials,os,dat);
  writeBinaryData(_bb,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> CompositeShapeExact::copy() const {
  std::shared_ptr<CompositeShapeExact> ret(new CompositeShapeExact);
  ret->_geoms.resize(_geoms.size());
  for(int i=0; i<(int)ret->_geoms.size(); i++)
    ret->_geoms[i]=std::dynamic_pointer_cast<ShapeExact>(_geoms[i]->copy());
  ret->_trans=_trans;
  ret->_materials=_materials;
  ret->_bb=_bb;
  return ret;
}
std::string CompositeShapeExact::type() const {
  return typeid(CompositeShapeExact).name();
}
const BBoxExact& CompositeShapeExact::getBB() const {
  return _bb;
}
bool CompositeShapeExact::empty() const {
  for(std::shared_ptr<ShapeExact> s:_geoms)
    if(!s->empty())
      return false;
  return true;
}
void CompositeShapeExact::getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                  std::vector<Eigen::Matrix<int,3,1>>& iss) const {
  std::vector<Eigen::Matrix<double,3,1>> vssTmp;
  std::vector<Eigen::Matrix<int,3,1>> issTmp;
  for(int i=0; i<(int)_geoms.size(); i++) {
    _geoms[i]->getMesh(vssTmp,issTmp);
    Eigen::Matrix<int,3,1> off((int)vss.size(),(int)vss.size(),(int)vss.size());
    for(const Eigen::Matrix<double,3,1>& v:vssTmp)
      vss.push_back(ROT(_trans[i]).template cast<double>()*v+CTR(_trans[i]).template cast<double>());
    for(const Eigen::Matrix<int,3,1>& i:issTmp)
      iss.push_back(i+off);
  }
}
bool CompositeShapeExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                                       T& rad,Eigen::Matrix<int,2,1>& feat,bool cache,
                                       std::vector<Vec3T>*) const {
  bool ret=true;
  Vec3T nTmp,normalTmp;
  Mat3T hessianTmp;
  T radTmp;
  Eigen::Matrix<int,2,1> featTmp;
  double minDist=std::numeric_limits<double>::max();  //this is not safe, reverting to double from rational
  for(int i=0; i<(int)_geoms.size(); i++) {
    bool retTmp=_geoms[i]->closestInner(ROT(_trans[i]).transpose()*(pt-CTR(_trans[i])),nTmp,normalTmp,hessianTmp,radTmp=0,featTmp,cache,NULL);
    double dist=nTmp.template cast<double>().norm()-(double)radTmp;  //this is not safe, reverting to double from rational
    if(dist<minDist) {
      n=ROT(_trans[i])*nTmp;
      normal=ROT(_trans[i])*normalTmp;
      hessian=ROT(_trans[i]).transpose()*hessianTmp*ROT(_trans[i]);
      rad=radTmp;
      feat=featTmp;
      ret=retTmp;
      minDist=dist;
    }
  }
  return ret;
}
void CompositeShapeExact::transform(const Mat3X4T& trans) {
  for(auto& t:_trans) {
    APPLY_TRANS(t,trans,t)
  }
  initBB();
}
void CompositeShapeExact::scale(T coef) {
  for(std::shared_ptr<ShapeExact> s:_geoms)
    s->scale(coef);
  for(Mat3X4T& t:_trans)
    CTR(t)*=coef;
}
const std::vector<std::shared_ptr<ShapeExact>>& CompositeShapeExact::getGeoms() const {
  return _geoms;
}
const std::vector<CompositeShapeExact::Mat3X4T>& CompositeShapeExact::getTrans() const {
  return _trans;
}
const std::vector<CompositeShapeExact::Material>& CompositeShapeExact::getMaterials() const {
  return _materials;
}
std::vector<CompositeShapeExact::Material>& CompositeShapeExact::getMaterials() {
  return _materials;
}
void CompositeShapeExact::writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  Mat3X4T t;
  for(int i=0; i<(int)_geoms.size(); i++) {
    APPLY_TRANS(t,trans,_trans[i]);
    _geoms[i]->writeVTK(os,t);
  }
}
//helper
void CompositeShapeExact::initBB() {
  _bb=BBoxExact();
  for(int i=0; i<(int)_geoms.size(); i++) {
    BBoxExact bb=_geoms[i]->getBB();
    for(const T& x: {
          bb.minCorner()[0],bb.maxCorner()[0]
        })
      for(const T& y: {
            bb.minCorner()[1],bb.maxCorner()[1]
          })
        for(const T& z: {
              bb.minCorner()[2],bb.maxCorner()[2]
            })
          _bb.setUnion(ROT(_trans[i])*Vec3T(x,y,z)+CTR(_trans[i]));
  }
}
void CompositeShapeExact::init(const aiScene* scene,const aiNode* node,bool buildBVH,Mat3X4T trans,const std::string& path) {
  if(!node)
    init(scene,scene->mRootNode,buildBVH,trans,path);
  else {
    Eigen::Matrix<double,4,4> transNode;
    //row0
    transNode(0,0)=node->mTransformation.a1;
    transNode(0,1)=node->mTransformation.a2;
    transNode(0,2)=node->mTransformation.a3;
    transNode(0,3)=node->mTransformation.a4;
    //row1
    transNode(1,0)=node->mTransformation.b1;
    transNode(1,1)=node->mTransformation.b2;
    transNode(1,2)=node->mTransformation.b3;
    transNode(1,3)=node->mTransformation.b4;
    //row2
    transNode(2,0)=node->mTransformation.c1;
    transNode(2,1)=node->mTransformation.c2;
    transNode(2,2)=node->mTransformation.c3;
    transNode(2,3)=node->mTransformation.c4;
    //row3
    transNode(3,0)=node->mTransformation.d1;
    transNode(3,1)=node->mTransformation.d2;
    transNode(3,2)=node->mTransformation.d3;
    transNode(3,3)=node->mTransformation.d4;
    //add children
    for(int i=0; i<(int)node->mNumChildren; i++)
      init(scene,node->mChildren[i],buildBVH,trans*transNode.template cast<T>(),path);
    //add mesh
    for(int idm=0; idm<(int)node->mNumMeshes; idm++) {
      const aiMesh* m=scene->mMeshes[node->mMeshes[idm]];
      const aiMaterial* mat=scene->mMaterials[m->mMaterialIndex];
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<double,2,1>> tcss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      for(int i=0; i<(int)m->mNumVertices; i++) {
        const aiVector3D& v=m->mVertices[i];
        vss.push_back(Eigen::Matrix<double,3,1>(v.x,v.y,v.z));
        if(m->mNumUVComponents[0]==2) {
          const aiVector3D& t=m->mTextureCoords[0][i];
          tcss.push_back(Eigen::Matrix<double,2,1>(t[0],1-t[1]));
        }
      }
      for(int i=0; i<(int)m->mNumFaces; i++) {
        const aiFace& f=m->mFaces[i];
        for(int j=0; j<(int)f.mNumIndices-2; j++)
          iss.push_back(Eigen::Matrix<int,3,1>(f.mIndices[0],f.mIndices[j+1],f.mIndices[j+2]));
      }
      _geoms.push_back(std::shared_ptr<ShapeExact>(new MeshExact(vss,tcss,iss,buildBVH)));
      _materials.push_back(initMaterial(*mat,path));
      _trans.push_back(trans*transNode.template cast<T>());
    }
  }
}
CompositeShapeExact::Material CompositeShapeExact::initMaterial(const aiMaterial& mat,const std::string& path) {
  Material ret;
  for(int idmat=0; idmat<(int)mat.mNumProperties; idmat++) {
    const aiMaterialProperty* matP=mat.mProperties[idmat];
    if(std::string(matP->mKey.C_Str()).find("clr.ambient")!=std::string::npos)
      ret._ambient=parseVec4(*matP);
    if(std::string(matP->mKey.C_Str()).find("clr.diffuse")!=std::string::npos)
      ret._diffuse=parseVec4(*matP);
    if(std::string(matP->mKey.C_Str()).find("clr.specular")!=std::string::npos)
      ret._specular=parseVec4(*matP);
    if(std::string(matP->mKey.C_Str()).find("mat.shininess")!=std::string::npos)
      ret._shininess=parseFloat(*matP);
    if(std::string(matP->mKey.C_Str()).find("mat.wireframe")!=std::string::npos)
      ret._useWireframe=parseInt(*matP);
  }
  if(mat.GetTextureCount(aiTextureType_DIFFUSE)>0) {
    aiString texPath;
    mat.GetTexture(aiTextureType_DIFFUSE,0,&texPath);
    ret._texFile=path+"/"+texPath.C_Str();
    ASSERT_MSGV(exists(ret._texFile),"Cannot file texture image: %s!",ret._texFile.c_str())
  }
  return ret;
}
Eigen::Matrix<double,4,1> CompositeShapeExact::parseVec4(const aiMaterialProperty& matP) {
  Eigen::Matrix<double,4,1> ret;
  if(matP.mType==aiPTI_Float) {
    int N=matP.mDataLength/4;
    for(int i=0; i<N; i++)
      ret[i]=parseFloat(matP,i);
    return ret;
  } else if(matP.mType==aiPTI_Double) {
    int N=matP.mDataLength/8;
    for(int i=0; i<N; i++)
      ret[i]=parseFloat(matP,i);
    return ret;
  } else {
    ASSERT_MSG(false,"Reading floating type aiMaterialProperty, but property type is not float!")
    return ret;
  }
}
double CompositeShapeExact::parseFloat(const aiMaterialProperty& matP,int index) {
  if(matP.mType==aiPTI_Float)
    return ((float*)matP.mData)[index];
  else if(matP.mType==aiPTI_Double)
    return ((double*)matP.mData)[index];
  else {
    ASSERT_MSG(false,"Reading floating type aiMaterialProperty, but property type is not float!")
    return 0;
  }
}
int CompositeShapeExact::parseInt(const aiMaterialProperty& matP,int index) {
  if(matP.mType==aiPTI_Integer)
    return ((int*)matP.mData)[index];
  else {
    ASSERT_MSG(false,"Reading integer type aiMaterialProperty, but property type is not integer!")
    return 0;
  }
}
}
