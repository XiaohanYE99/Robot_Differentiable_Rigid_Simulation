#ifndef ENVIRONMENT_VISUALIZER_H
#define ENVIRONMENT_VISUALIZER_H

#include "Environment.h"
#include "ContactGenerator.h"
#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/MeshShape.h>
#include <TinyVisualizer/TerrainShape.h>
#include <TinyVisualizer/Bullet3DShape.h>

namespace PHYSICSMOTION {
struct TriangleExact;
struct MeshExact;
struct PointCloudExact;
struct BBoxExact;
struct SphericalBBoxExact;
struct ShapeExact;
extern std::shared_ptr<DRAWER::Shape> visualizeTriangleExact(std::shared_ptr<TriangleExact> et,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizeMeshExact(std::shared_ptr<MeshExact> m,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizePointCloudExact(std::shared_ptr<PointCloudExact> m,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizeBBoxExact(std::shared_ptr<BBoxExact> m,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizeSphericalBBoxExact(std::shared_ptr<SphericalBBoxExact> m,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizeShapeExactGetMesh(std::shared_ptr<ShapeExact> m,bool wire=false);  //visualize using the getMesh function, not recommended. Currently I use it for debug.
extern std::shared_ptr<DRAWER::Shape> visualizeShapeExact(std::shared_ptr<ShapeExact> m,bool wire=false);
template <typename T>
std::shared_ptr<DRAWER::Shape> visualizeEnvironmentExact(std::shared_ptr<EnvironmentExact<T>> e,const Eigen::Matrix<GLfloat,2,1>& tcMult) {
  using namespace DRAWER;
  const MeshExact& m=e->getMesh();
  std::shared_ptr<MeshShape> tri(new MeshShape);
  for(int i=0; i<(int)m.vss().size(); i++) {
    Eigen::Matrix<GLfloat,-1,1> tc;
    tc.resize(2);
    tc[0]=(GLfloat)m.vss()[i][0]*tcMult[0];
    tc[1]=(GLfloat)m.vss()[i][1]*tcMult[1];
    tri->addVertex(m.vss()[i].cast<GLfloat>(),&tc);
  }
  for(int i=0; i<(int)m.iss().size(); i++)
    tri->addIndex(m.iss()[i].template cast<GLuint>());
  tri->setMode(GL_TRIANGLES);
  tri->computeNormals();
  return tri;
}
template <typename T>
std::shared_ptr<DRAWER::Shape> visualizeEnvironmentHeight(std::shared_ptr<EnvironmentHeight<T>> e,const Eigen::Matrix<GLfloat,2,1>& tcMult) {
  using namespace DRAWER;
  Eigen::Matrix<GLfloat,-1,-1> height=e->getHeightMatrix().template cast<GLfloat>();
  std::shared_ptr<TerrainShape> shape(new TerrainShape(height,2,Eigen::Matrix<GLfloat,3,1>(1,1,1),tcMult));
  //std::cout << shape->getBB().transpose() << std::endl;
  std::shared_ptr<Bullet3DShape> trans(new Bullet3DShape);
  trans->addShape(shape);

  BBoxExact bb=e->getBB();
  Eigen::Matrix<GLfloat,6,1> bbS=shape->getBB();
  BBoxExact::Vec3T ctr=(bb.maxCorner()+bb.minCorner())/2,rng=(bb.maxCorner()-bb.minCorner())/2;
  Eigen::Matrix<GLfloat,3,1> ctrS=(bbS.segment<3>(3)+bbS.segment<3>(0))/2,rngS=(bbS.segment<3>(3)-bbS.segment<3>(0))/2;
  for(int i=0; i<3; i++)
    if(rngS[i]<=0)
      rng[i]=rngS[i]=1;

  Eigen::Matrix<GLfloat,4,4> transM1;
  transM1.setIdentity();
  transM1.block<3,1>(0,3)=-ctrS.template cast<GLfloat>();

  Eigen::Matrix<GLfloat,4,4> transM2;
  transM2.setIdentity();
  transM2.block<3,3>(0,0).diagonal().array()=rng.array().template cast<GLfloat>()/rngS.array();

  Eigen::Matrix<GLfloat,4,4> transM3;
  transM3.setIdentity();
  transM3.block<3,1>(0,3)=ctr.template cast<GLfloat>();
  trans->setLocalTransform(transM3*transM2*transM1);
  return trans;
}
template <typename T>
std::shared_ptr<DRAWER::Shape> visualizeEnvironment(std::shared_ptr<Environment<T>> e,const Eigen::Matrix<GLfloat,2,1>& tcMult=Eigen::Matrix<GLfloat,2,1>(1,1)) {
  std::shared_ptr<DRAWER::CompositeShape> shape(new DRAWER::CompositeShape);
  if(std::dynamic_pointer_cast<EnvironmentExact<T>>(e))
    shape->addShape(visualizeEnvironmentExact(std::dynamic_pointer_cast<EnvironmentExact<T>>(e),tcMult));
  else if(std::dynamic_pointer_cast<EnvironmentHeight<T>>(e))
    shape->addShape(visualizeEnvironmentHeight(std::dynamic_pointer_cast<EnvironmentHeight<T>>(e),tcMult));
  return shape;
}
extern std::shared_ptr<DRAWER::Shape> visualizeEnvironment(std::unordered_set<std::shared_ptr<ShapeExact>> obs,bool wire=false);
extern std::shared_ptr<DRAWER::Shape> visualizeEnvironment(std::vector<std::shared_ptr<ShapeExact>> obs,bool wire=false);
extern void visualizeContact(std::shared_ptr<DRAWER::CompositeShape> shape,const std::vector<ContactGenerator::ContactManifold>& manifolds,int sz=5,GLfloat forceCoef=1.f);
extern std::shared_ptr<DRAWER::Shape> visualizeContact(const std::vector<ContactGenerator::ContactManifold>& manifolds,int sz=5);
}

#endif
