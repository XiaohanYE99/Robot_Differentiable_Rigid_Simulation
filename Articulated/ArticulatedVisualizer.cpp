#include "ArticulatedVisualizer.h"
#include "ArticulatedLoader.h"
#include <Environment/EnvironmentVisualizer.h>
#include <TinyVisualizer/Bullet3DShape.h>
#include <TinyVisualizer/CellShape.h>
#include <TinyVisualizer/MakeMesh.h>

namespace PHYSICSMOTION {
std::shared_ptr<DRAWER::Shape> visualizeArticulated
(std::shared_ptr<ArticulatedBody> b,
 const Eigen::Matrix<GLfloat,3,1>& colorBody,
 const Eigen::Matrix<GLfloat,3,1>& colorBB,
 const Eigen::Matrix<GLfloat,3,1>& colorIDX,
 bool wire) {
  using namespace DRAWER;
  std::shared_ptr<CompositeShape> shape(new CompositeShape);
  for(int i=0; i<b->nrJ(); i++)
    if(b->joint(i)._mesh) {
      std::shared_ptr<Bullet3DShape> stJ(new Bullet3DShape);
      {
        std::shared_ptr<Shape> s=visualizeShapeExact(b->joint(i)._mesh,wire);
        std::shared_ptr<Bullet3DShape> st(new Bullet3DShape);
        st->setLocalRotate(ROT(b->joint(i)._transMesh).template cast<GLfloat>());
        st->setLocalTranslate(CTR(b->joint(i)._transMesh).template cast<GLfloat>());
        st->addShape(s);
        stJ->addShape(st);
        st->setColorDiffuse(GL_TRIANGLES,colorBody[0],colorBody[1],colorBody[2]);
        st->setColorDiffuse(GL_LINES,colorBody[0],colorBody[1],colorBody[2]);
      }
      shape->addShape(stJ);
    }
  return shape;
}
void highlightJoint
(std::shared_ptr<DRAWER::Shape> root,std::shared_ptr<ArticulatedBody> b,
 const std::vector<char>* selected,
 const std::vector<char>* excluded,
 const std::vector<char>* markFoot,
 const Eigen::Matrix<GLfloat,3,1>& colorBody) {
  std::function<bool(int)> isSelected=[&](int cid)->bool{
    for(; cid>=0; cid=b->joint(cid)._parent)
      if(selected&&selected->at(cid))
        return true;
    return false;
  };
  std::function<bool(int)> isExcluded=[&](int cid)->bool{
    for(; cid>=0; cid=b->joint(cid)._parent)
      if(excluded&&excluded->at(cid))
        return true;
    return false;
  };
  std::function<bool(int)> isFoot=[&](int cid)->bool{
    return markFoot&&markFoot->at(cid);
  };
  std::shared_ptr<DRAWER::CompositeShape> croot=std::dynamic_pointer_cast<DRAWER::CompositeShape>(root);
  for(int i=0,j=0; i<b->nrJ(); i++)
    if(b->joint(i)._mesh) {
      std::shared_ptr<DRAWER::CompositeShape> shapeJ=std::dynamic_pointer_cast<DRAWER::CompositeShape>(croot->getChild(j++));
      std::shared_ptr<DRAWER::Shape> st=shapeJ->getChild(0);
      st->setEnabled(!isExcluded(i));
      st->setColorDiffuse(GL_TRIANGLES,isSelected(i)?1:colorBody[0],colorBody[1],isFoot(i)?1:colorBody[2]);
      st->setColorDiffuse(GL_LINES,isSelected(i)?1:colorBody[0],colorBody[1],isFoot(i)?1:colorBody[2]);
    }
}
int getRayIntersectedJoint(std::shared_ptr<DRAWER::Shape> root,std::shared_ptr<ArticulatedBody> b,const Eigen::Matrix<GLfloat,6,1>& ray) {
  std::shared_ptr<DRAWER::CompositeShape> croot=std::dynamic_pointer_cast<DRAWER::CompositeShape>(root);
  int sjid=-1;
  GLfloat alpha=1e6;
  for(int i=0,j=0; i<b->nrJ(); i++)
    if(b->joint(i)._mesh) {
      std::shared_ptr<DRAWER::CompositeShape> shapeJ=std::dynamic_pointer_cast<DRAWER::CompositeShape>(croot->getChild(j++));
      if(shapeJ->rayIntersect(ray,alpha))
        sjid=i;
    }
  return sjid;
}
void updateArticulatedBody(std::shared_ptr<DRAWER::Shape> shape,std::shared_ptr<ArticulatedBody> b,const ArticulatedBody::Mat3XT& M) {
  using namespace DRAWER;
  std::shared_ptr<CompositeShape> shapeC=std::dynamic_pointer_cast<CompositeShape>(shape);
  for(int i=0,j=0; i<b->nrJ(); i++)
    if(b->joint(i)._mesh) {
      std::shared_ptr<Bullet3DShape> s=std::dynamic_pointer_cast<Bullet3DShape>(shapeC->getChild(j));
      ArticulatedBody::Mat3X4T TI=TRANSI(M,i);
      s->setLocalRotate(ROT(TI).template cast<GLfloat>());
      s->setLocalTranslate(CTR(TI).template cast<GLfloat>());
      j++;
    }
}
void updateArticulatedBody(std::shared_ptr<DRAWER::Shape> shape,std::shared_ptr<ArticulatedBody> b,const ArticulatedBody::Vec& x) {
  updateArticulatedBody(shape,b,b->getT(x));
}
void compareVisualMeshURDF(int argc,char** argv,const std::string& file,bool wire,
                           const Eigen::Matrix<GLfloat,3,1>& color) {
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  //read links
  ArticulatedLoader::Links links,linksRef;
  ArticulatedLoader::readLinks(linksRef,*(pt.RootElement()),std::filesystem::path(file).parent_path(),false);
  ArticulatedLoader::readLinks(links,*(pt.RootElement()),std::filesystem::path(file).parent_path(),true);
  //visualMesh
  using namespace DRAWER;
  Drawer drawer(argc,argv);
  drawer.addCamera3D(90);
  GLfloat deltaSum=0,deltaMax=0;
  std::vector<std::shared_ptr<Shape>> sRefs,ss;
  for(const std::pair<const std::string,Joint>&v:links) {
    std::shared_ptr<Shape> sRef=visualizeShapeExact(linksRef[v.first]._mesh,wire);
    GLfloat offRef=(sRef->getBB().segment<3>(3)-sRef->getBB().segment<3>(0)).norm();
    std::shared_ptr<Shape> s=visualizeShapeExact(links[v.first]._mesh,wire);
    GLfloat off=(s->getBB().segment<3>(3)-s->getBB().segment<3>(0)).norm();
    deltaMax=std::max(deltaMax,off+offRef);
    sRefs.push_back(sRef);
    ss.push_back(s);
  }
  for(int i=0; i<(int)sRefs.size(); i++) {
    std::shared_ptr<Bullet3DShape> sRefT(new Bullet3DShape);
    sRefT->addShape(sRefs[i]);
    sRefT->setColorDiffuse(GL_TRIANGLES,color[0],color[1],color[2]);
    sRefT->setColorDiffuse(GL_LINES,color[0],color[1],color[2]);
    drawer.addShape(sRefT);

    std::shared_ptr<Bullet3DShape> sT(new Bullet3DShape);
    sT->addShape(ss[i]);
    sT->setColorDiffuse(GL_TRIANGLES,color[0],color[1],color[2]);
    sT->setColorDiffuse(GL_LINES,color[0],color[1],color[2]);
    drawer.addShape(sT);

    sRefT->setLocalTranslate(Eigen::Matrix<GLfloat,3,1>::UnitX()*deltaSum);
    sT->setLocalTranslate(Eigen::Matrix<GLfloat,3,1>::UnitX()*deltaSum+Eigen::Matrix<GLfloat,3,1>::UnitZ()*deltaMax);
    deltaSum+=deltaMax;
  }
  //frame
  int off=0;
  std::shared_ptr<MeshShape> frame(new MeshShape);
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>(0,0,0));
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>((int)sRefs.size()*deltaMax,0,0));
  frame->addIndex(Eigen::Matrix<GLuint,2,1>(off+0,off+1));
  off+=2;
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>(0,0,deltaMax));
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>((int)sRefs.size()*deltaMax,0,deltaMax));
  frame->addIndex(Eigen::Matrix<GLuint,2,1>(off+0,off+1));
  off+=2;
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>(0,0,deltaMax*2));
  frame->addVertex(Eigen::Matrix<GLfloat,3,1>((int)sRefs.size()*deltaMax,0,deltaMax*2));
  frame->addIndex(Eigen::Matrix<GLuint,2,1>(off+0,off+1));
  off+=2;
  for(int i=0; i<=(int)sRefs.size(); i++) {
    frame->addVertex(Eigen::Matrix<GLfloat,3,1>(i*deltaMax,0,0));
    frame->addVertex(Eigen::Matrix<GLfloat,3,1>(i*deltaMax,0,deltaMax*2));
    frame->addIndex(Eigen::Matrix<GLuint,2,1>(off+0,off+1));
    off+=2;
  }
  frame->setMode(GL_LINES);
  frame->setColorDiffuse(GL_LINES,0,0,0);
  std::shared_ptr<Bullet3DShape> frameT(new Bullet3DShape);
  frameT->addShape(frame);
  frameT->setLocalTranslate(Eigen::Matrix<GLfloat,3,1>(-deltaMax/2,0,-deltaMax/2));
  drawer.addShape(frameT);
  drawer.mainLoop();
}
}
