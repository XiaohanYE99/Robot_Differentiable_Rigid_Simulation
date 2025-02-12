#include <Environment/SAT.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/EnvironmentVisualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>
#include <Utils/RotationUtils.h>
#include <random>

using namespace DRAWER;
using namespace PHYSICSMOTION;

void visualizeSATFeature(std::shared_ptr<MeshShape> shape,
                         const SAT::Mat3X4T& t,const ShapeExact::Facet& f) {
  shape->clear();
  for(int i=2,id=0; i<(int)f._boundary.size(); i++) {
    shape->addVertex((ROT(t)*f._boundary[0  ]+CTR(t)).template cast<GLfloat>());
    shape->addVertex((ROT(t)*f._boundary[i-1]+CTR(t)).template cast<GLfloat>());
    shape->addVertex((ROT(t)*f._boundary[i  ]+CTR(t)).template cast<GLfloat>());
    shape->addIndexSingle(id++);
    shape->addIndexSingle(id++);
    shape->addIndexSingle(id++);
  }
  shape->setMode(GL_TRIANGLES);
  shape->setColorDiffuse(GL_TRIANGLES,1,0,0);
  shape->computeNormals();
}
void visualizeSATFeature(std::shared_ptr<MeshShape> shape,
                         const SAT::Mat3X4T& t1,const ShapeExact::Edge& e1,
                         const SAT::Mat3X4T& t2,const ShapeExact::Edge& e2) {
  shape->clear();
  shape->addVertex((ROT(t1)*e1._a+CTR(t1)).template cast<GLfloat>());
  shape->addVertex((ROT(t1)*e1._b+CTR(t1)).template cast<GLfloat>());
  shape->addVertex((ROT(t2)*e2._a+CTR(t2)).template cast<GLfloat>());
  shape->addVertex((ROT(t2)*e2._b+CTR(t2)).template cast<GLfloat>());
  shape->addIndexSingle(0);
  shape->addIndexSingle(1);
  shape->addIndexSingle(2);
  shape->addIndexSingle(3);
  shape->setMode(GL_LINES);
  shape->setColorDiffuse(GL_LINES,1,0,0);
  shape->setLineWidth(5);
}
void mainSAT(std::shared_ptr<ShapeExact> s1,
             std::shared_ptr<ShapeExact> s2,
             int argc,char** argv,bool swap) {
  SAT::Vec3T r1=SAT::Vec3T::Identity(),t1=SAT::Vec3T::Zero();
  SAT::Vec3T r2=SAT::Vec3T::Identity(),t2=SAT::Vec3T::Zero();
  std::vector<ShapeExact::Facet> F1=s1->facets(),F2=s2->facets();
  std::vector<ShapeExact::Edge> E1=s1->edges(),E2=s2->edges();
  SAT::T bary=1,speed=0.01,rng=5;

  Drawer drawer(argc,argv);
  std::shared_ptr<Bullet3DShape> shape1(new Bullet3DShape);
  std::shared_ptr<Bullet3DShape> shape2(new Bullet3DShape);
  std::shared_ptr<CompositeShape> manifold(new CompositeShape);
  std::shared_ptr<MeshShape> feature(new MeshShape);
  shape1->addShape(visualizeShapeExact(s1,true));
  shape2->addShape(visualizeShapeExact(s2,true));
  shape1->setLineWidth(1);
  shape2->setLineWidth(1);
  drawer.addShape(shape1);
  drawer.addShape(shape2);
  drawer.addShape(manifold);
  drawer.addShape(feature);

  bool sim=true,intersect;
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.setKeyFunc([&](GLFWwindow* wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.setFrameFunc([&](std::shared_ptr<SceneNode>&) {
    if(sim) {
      if(bary<1)
        bary+=speed;
      else {
        bary=0;
        r1=r2;
        t1=t2;
        r2=SAT::Vec3T::Random()*M_PI;
        t2=SAT::Vec3T::Random()*rng;
      }
    }
    //set transform
    auto r=expWGradV<GLfloat,Eigen::Matrix<GLfloat,3,1>>((r1*(1-bary)+r2*bary).template cast<GLfloat>(),NULL,NULL);
    auto t=(t1*(1-bary)+t2*bary).template cast<GLfloat>();
    shape2->setLocalRotate(r);
    shape2->setLocalTranslate(t);
    //SAT
    SAT::ContactManifold m;
    SAT::Mat3X4T T1=SAT::Mat3X4T::Identity(),T2;
    ROT(T2)=r.template cast<SAT::T>();
    CTR(T2)=t.template cast<SAT::T>();
    SAT::ProjRange rngP=SAT::generateManifold(s1,s2,F1,F2,E1,E2,T1,T2,m,&intersect);
    if(intersect)
      shape2->setColorDiffuse(GL_LINES,1,0,0);
    else shape2->setColorDiffuse(GL_LINES,0,1,0);
    //draw contact
    visualizeContact(manifold, {m});
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
    //draw SAT feature
    feature->clear();
    if(rngP._fidA>=0)
      visualizeSATFeature(feature,T1,F1[rngP._fidA]);
    else if(rngP._fidB>=0)
      visualizeSATFeature(feature,T2,F2[rngP._fidB]);
    else if(rngP._eidA>=0 && rngP._eidB>=0)
      visualizeSATFeature(feature,T1,E1[rngP._eidA],T2,E2[rngP._eidB]);
  });
  drawer.mainLoop();
}
int main(int argc,char** argv) {
  {
    std::shared_ptr<BBoxExact> s1(new BBoxExact(1,2,3));
    std::shared_ptr<BBoxExact> s2(new BBoxExact(.5,.5,.5));
    mainSAT(s1,s2,argc,argv,false);
  }
  {
    std::shared_ptr<BBoxExact> s1(new BBoxExact(.5,.5,.5));
    std::shared_ptr<BBoxExact> s2(new BBoxExact(1,2,3));
    mainSAT(s1,s2,argc,argv,true);
  }
  {
    std::shared_ptr<BBoxExact> s1(new BBoxExact(1,2,3));
    std::shared_ptr<BBoxExact> s2(new SphericalBBoxExact(2,.5));
    mainSAT(s1,s2,argc,argv,false);
  }
  {
    std::shared_ptr<BBoxExact> s1(new SphericalBBoxExact(2,.5));
    std::shared_ptr<BBoxExact> s2(new BBoxExact(1,2,3));
    mainSAT(s1,s2,argc,argv,false);
  }
  return 0;
}
