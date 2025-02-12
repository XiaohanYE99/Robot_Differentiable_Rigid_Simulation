#include "EnvironmentVisualizer.h"
#include "MeshExact.h"
#include "PointCloudExact.h"
#include "SphericalBBoxExact.h"
#include "CompositeShapeExact.h"
#include "TriangleExact.h"
#include "BBoxExact.h"

#include <TinyVisualizer/Bullet3DShape.h>
#include <TinyVisualizer/MakeMesh.h>
#include <TinyVisualizer/MakeTexture.h>
#include <Utils/Epsilon.h>

namespace PHYSICSMOTION {
std::shared_ptr<DRAWER::Shape> visualizeTriangleExact(std::shared_ptr<TriangleExact> et,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<MeshShape> tri(new MeshShape);
  tri->addVertex(et->v(0).cast<GLfloat>());
  tri->addVertex(et->v(1).cast<GLfloat>());
  tri->addVertex(et->v(2).cast<GLfloat>());
  if(wire) {
    tri->addIndex(Eigen::Matrix<GLuint,2,1>(0,1));
    tri->addIndex(Eigen::Matrix<GLuint,2,1>(0,2));
    tri->addIndex(Eigen::Matrix<GLuint,2,1>(1,2));
    tri->setMode(GL_LINES);
  } else {
    tri->addIndex(Eigen::Matrix<GLuint,3,1>(0,1,2));
    tri->setMode(GL_TRIANGLES);
  }
  return tri;
}
std::shared_ptr<DRAWER::Shape> visualizeMeshExact(std::shared_ptr<MeshExact> m,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<MeshShape> tri(new MeshShape);
  for(int i=0; i<(int)m->vss().size(); i++)
    if(!m->tcss().empty()) {
      Eigen::Matrix<GLfloat,-1,1> tc=m->tcss()[i].cast<GLfloat>();
      tri->addVertex(m->vss()[i].cast<GLfloat>(),&tc);
    } else tri->addVertex(m->vss()[i].cast<GLfloat>());
  if(wire) {
    for(int i=0; i<(int)m->iss().size(); i++) {
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(m->iss()[i][0],m->iss()[i][1]));
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(m->iss()[i][1],m->iss()[i][2]));
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(m->iss()[i][2],m->iss()[i][0]));
    }
    tri->setMode(GL_LINES);
  } else {
    for(int i=0; i<(int)m->iss().size(); i++)
      tri->addIndex(m->iss()[i].template cast<GLuint>());
    tri->setMode(GL_TRIANGLES);
  }
  if(!wire)
    tri->computeNormals();
  return tri;
}
std::shared_ptr<DRAWER::Shape> visualizePointCloudExact(std::shared_ptr<PointCloudExact> m,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<MeshShape> pc(new MeshShape);
  for(int i=0; i<(int)m->vss().size(); i++) {
    pc->addVertex(m->vss()[i].cast<GLfloat>());
    pc->addIndexSingle(i);
  }
  pc->setMode(GL_POINTS);
  return pc;
}
std::shared_ptr<DRAWER::Shape> visualizeBBoxExact(std::shared_ptr<BBoxExact> m,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<Texture> checker=drawGrid(25,0.001,0,Eigen::Matrix<GLfloat,3,1>(.8f,.9f,1.f),Eigen::Matrix<GLfloat,3,1>(.1f,.1f,.1f));
  std::shared_ptr<Bullet3DShape> shape(new Bullet3DShape);
  std::shared_ptr<MeshShape> box=makeBox(1,!wire,(m->maxCorner()-m->minCorner()).template cast<GLfloat>()/2);
  box->setCastShadow(false);
  box->setTextureDiffuse(checker);
  box->setColorDiffuse(GL_POINTS,.8f,.9f,1.f);
  box->setColorDiffuse(GL_LINES,.8f,.9f,1.f);
  box->setColorDiffuse(GL_TRIANGLES,.8f,.9f,1.f);
  shape->setLocalTranslate((m->maxCorner()+m->minCorner()).template cast<GLfloat>()/2);
  shape->addShape(box);
  return shape;
}
std::shared_ptr<DRAWER::Shape> visualizeSphericalBBoxExact(std::shared_ptr<SphericalBBoxExact> m,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<Bullet3DShape> shape(new Bullet3DShape);
  std::shared_ptr<MeshShape> sbox=makeSphericalBox(8,!wire,(GLfloat)m->radius(),(m->maxCorner()-m->minCorner()).template cast<GLfloat>()/2);
  shape->setLocalTranslate((m->maxCorner()+m->minCorner()).template cast<GLfloat>()/2);
  shape->addShape(sbox);
  return shape;
}
std::shared_ptr<DRAWER::Shape> visualizeShapeExactGetMesh(std::shared_ptr<ShapeExact> m,bool wire) {
  using namespace DRAWER;
  std::shared_ptr<MeshShape> tri(new MeshShape);
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  m->getMesh(vss,iss);
  for(int i=0; i<(int)vss.size(); i++)
    tri->addVertex(vss[i].cast<GLfloat>());
  if(wire) {
    for(int i=0; i<(int)iss.size(); i++) {
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(iss[i][0],iss[i][1]));
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(iss[i][1],iss[i][2]));
      tri->addIndex(Eigen::Matrix<GLuint,2,1>(iss[i][2],iss[i][0]));
    }
    tri->setMode(GL_LINES);
  } else {
    for(int i=0; i<(int)iss.size(); i++)
      tri->addIndex(iss[i].template cast<GLuint>());
    tri->setMode(GL_TRIANGLES);
  }
  return tri;
}
std::shared_ptr<DRAWER::Shape> visualizeShapeExact(std::shared_ptr<ShapeExact> m,bool wire) {
  using namespace DRAWER;
  if(std::dynamic_pointer_cast<MeshExact>(m))
    return visualizeMeshExact(std::dynamic_pointer_cast<MeshExact>(m),wire);
  else if(std::dynamic_pointer_cast<PointCloudExact>(m))
    return visualizePointCloudExact(std::dynamic_pointer_cast<PointCloudExact>(m),wire);
  else if(std::dynamic_pointer_cast<SphericalBBoxExact>(m))
    return visualizeSphericalBBoxExact(std::dynamic_pointer_cast<SphericalBBoxExact>(m),wire);
  else if(std::dynamic_pointer_cast<BBoxExact>(m))
    return visualizeBBoxExact(std::dynamic_pointer_cast<BBoxExact>(m),wire);
  else if(std::dynamic_pointer_cast<CompositeShapeExact>(m)) {
    Eigen::Matrix<GLfloat,4,4> localT;
    localT.setIdentity();
    std::shared_ptr<Bullet3DShape> shape(new Bullet3DShape);
    const std::vector<std::shared_ptr<ShapeExact>>& geoms=std::dynamic_pointer_cast<CompositeShapeExact>(m)->getGeoms();
    const std::vector<CompositeShapeExact::Mat3X4T>& trans=std::dynamic_pointer_cast<CompositeShapeExact>(m)->getTrans();
    const std::vector<CompositeShapeExact::Material>& materials=std::dynamic_pointer_cast<CompositeShapeExact>(m)->getMaterials();
    for(int i=0; i<(int)geoms.size(); i++) {
      std::shared_ptr<Bullet3DShape> shapeI(new Bullet3DShape);
      //set material
      if((int)materials.size()>i)
        wire=materials[i]._useWireframe;
      shapeI->addShape(visualizeShapeExact(geoms[i],wire));
      //set material
      if((int)materials.size()>i) {
        //ambient
        shapeI->setColorAmbient(GL_POINTS,materials[i]._ambient[0],materials[i]._ambient[1],materials[i]._ambient[2]);
        shapeI->setColorAmbient(GL_LINES,materials[i]._ambient[0],materials[i]._ambient[1],materials[i]._ambient[2]);
        shapeI->setColorAmbient(GL_TRIANGLES,materials[i]._ambient[0],materials[i]._ambient[1],materials[i]._ambient[2]);
        //diffuse
        shapeI->setColorDiffuse(GL_POINTS,materials[i]._diffuse[0],materials[i]._diffuse[1],materials[i]._diffuse[2]);
        shapeI->setColorDiffuse(GL_LINES,materials[i]._diffuse[0],materials[i]._diffuse[1],materials[i]._diffuse[2]);
        shapeI->setColorDiffuse(GL_TRIANGLES,materials[i]._diffuse[0],materials[i]._diffuse[1],materials[i]._diffuse[2]);
        //specular
        shapeI->setColorSpecular(GL_POINTS,materials[i]._specular[0],materials[i]._specular[1],materials[i]._specular[2]);
        shapeI->setColorSpecular(GL_LINES,materials[i]._specular[0],materials[i]._specular[1],materials[i]._specular[2]);
        shapeI->setColorSpecular(GL_TRIANGLES,materials[i]._specular[0],materials[i]._specular[1],materials[i]._specular[2]);
        //shininess
        shapeI->setShininess(GL_POINTS,materials[i]._shininess);
        shapeI->setShininess(GL_LINES,materials[i]._shininess);
        shapeI->setShininess(GL_TRIANGLES,materials[i]._shininess);
        //texture
        if(!materials[i]._texFile.empty()) {
          std::shared_ptr<DRAWER::Texture> tex=DRAWER::Texture::load(materials[i]._texFile);
          std::cout << "Loaded texture from " << materials[i]._texFile << " w=" <<  tex->width() << " h=" << tex->height() << std::endl;
          shapeI->setTextureDiffuse(tex);
        }
      }
      //assign
      localT.block<3,4>(0,0)=trans[i].template cast<GLfloat>();
      shapeI->setLocalTransform(localT);
      shape->addShape(shapeI);
    }
    return shape;
  } else return NULL;
}
std::shared_ptr<DRAWER::Shape> visualizeEnvironment(std::unordered_set<std::shared_ptr<ShapeExact>> obs,bool wire) {
  std::shared_ptr<DRAWER::CompositeShape> shape(new DRAWER::CompositeShape);
  for(std::shared_ptr<ShapeExact> o:obs)
    shape->addShape(visualizeShapeExact(o,wire));
  return shape;
}
std::shared_ptr<DRAWER::Shape> visualizeEnvironment(std::vector<std::shared_ptr<ShapeExact>> obs,bool wire) {
  std::shared_ptr<DRAWER::CompositeShape> shape(new DRAWER::CompositeShape);
  for(std::shared_ptr<ShapeExact> o:obs)
    shape->addShape(visualizeShapeExact(o,wire));
  return shape;
}
void visualizeContact(std::shared_ptr<DRAWER::CompositeShape> shape,const std::vector<ContactGenerator::ContactManifold>& manifolds,int sz,GLfloat forceCoef) {
  std::shared_ptr<DRAWER::MeshShape> shapeA(new DRAWER::MeshShape);
  std::shared_ptr<DRAWER::MeshShape> shapeB(new DRAWER::MeshShape);
  std::shared_ptr<DRAWER::MeshShape> shapeL(new DRAWER::MeshShape);
  std::shared_ptr<DRAWER::MeshShape> shapeF(new DRAWER::MeshShape);
  int id=0,idF=0;
  for(const auto& m:manifolds)
    for(const auto& p:m._points) {
      shapeA->addVertex(p._ptA.template cast<GLfloat>());
      shapeA->addIndexSingle(id);
      shapeB->addVertex(p._ptB.template cast<GLfloat>());
      shapeB->addIndexSingle(id);
      shapeL->addVertex(p._ptA.template cast<GLfloat>());
      shapeL->addVertex(p._ptB.template cast<GLfloat>());
      shapeL->addIndexSingle(id*2+0);
      shapeL->addIndexSingle(id*2+1);
      if(m._jidA>=0 && forceCoef>0) {
        shapeF->addVertex(p._ptA.template cast<GLfloat>());
        shapeF->addVertex((p._ptA+p._fA*forceCoef).template cast<GLfloat>());
        shapeL->addIndexSingle(idF*2+0);
        shapeL->addIndexSingle(idF*2+1);
        idF++;
      }
      if(m._jidB>=0 && forceCoef>0) {
        shapeF->addVertex(p._ptB.template cast<GLfloat>());
        shapeF->addVertex((p._ptB+p._fB*forceCoef).template cast<GLfloat>());
        shapeL->addIndexSingle(idF*2+0);
        shapeL->addIndexSingle(idF*2+1);
        idF++;
      }
      id++;
    }
  shapeA->setMode(GL_POINTS);
  shapeB->setMode(GL_POINTS);
  shapeL->setMode(GL_LINES);
  shapeF->setMode(GL_LINES);

  shapeA->setColorDiffuse(GL_POINTS,1,0,0);
  shapeB->setColorDiffuse(GL_POINTS,0,1,0);
  shapeL->setColorDiffuse(GL_LINES, 0,0,1);
  shapeF->setColorDiffuse(GL_LINES, 1.,.5,.5);

  shapeA->setPointSize(sz);
  shapeB->setPointSize(sz);
  shapeL->setLineWidth(sz);
  shapeF->setLineWidth(sz);
  //shape
  while(shape->numChildren()>0)
    shape->removeChild(shape->getChild(0));
  shape->addShape(shapeA);
  shape->addShape(shapeB);
  shape->addShape(shapeL);
  shape->addShape(shapeF);
}
std::shared_ptr<DRAWER::Shape> visualizeContact(const std::vector<ContactGenerator::ContactManifold>& manifolds,int sz) {
  std::shared_ptr<DRAWER::CompositeShape> shape(new DRAWER::CompositeShape);
  visualizeContact(shape,manifolds,sz);
  return shape;
}
}
