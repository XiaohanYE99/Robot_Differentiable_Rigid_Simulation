#include "Environment/CompositeShapeExact.h"
#include <Environment/ContactGenerator.h>
#include <Environment/BBoxExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/EnvironmentVisualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>
#include <Utils/RotationUtils.h>
#include <random>

using namespace PHYSICSMOTION;

int off=0,dx=10;
void SS1(std::vector<std::shared_ptr<ShapeExact>>& shapes) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=CTR(t2)=CompositeShapeExact::Vec3T(.1,.1,.1)+CompositeShapeExact::Vec3T::UnitY()*off;
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  off+=dx;
}
void SS2(std::vector<std::shared_ptr<ShapeExact>>& shapes) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(0.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=ROT(t1)*CompositeShapeExact::Vec3T( .1,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(1.1,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  off+=dx;
}
void SC1(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY,bool swap=false) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(.5,1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=ROT(t1)*CompositeShapeExact::Vec3T(.0,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  if(swap) {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  } else {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  }
  off+=dx;
}
void CC1(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(.5,1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(1,.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=ROT(t1)*CompositeShapeExact::Vec3T(.0,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  off+=dx;
}
void CC2(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY,double offZ) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(.5,1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(1,.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1).setIdentity();
  ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>(0,0,M_PI/3),NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=ROT(t1)*CompositeShapeExact::Vec3T(.0,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,offZ)+CompositeShapeExact::Vec3T::UnitY()*off;
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  off+=dx;
}
void CC3(std::vector<std::shared_ptr<ShapeExact>>& shapes,double angle,double offY,bool swap=false) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(.5,1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(.2,.1));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1).setIdentity();
  ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>(0,0,angle),NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=ROT(t1)*CompositeShapeExact::Vec3T(.0,.0,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(.0,offY,.0)+CompositeShapeExact::Vec3T::UnitY()*off;
  if(swap) {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  } else {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  }
  off+=dx;
}
void SB1(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY,double offZ,bool swap=false) {
  std::shared_ptr<BBoxExact> s1(new BBoxExact(1,1.1,1.2));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=ROT(t2)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  CTR(t1)=CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,offZ)+CompositeShapeExact::Vec3T::UnitY()*off;
  if(swap) {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  } else {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  }
  off+=dx;
}
void CB1(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY,double offZ,bool swap=false) {
  std::shared_ptr<BBoxExact> s1(new BBoxExact(1,1.1,1.2));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(1,.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  ROT(t2)=ROT(t1)*eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>(0,M_PI/4,0),NULL,NULL).template cast<CompositeShapeExact::T>();;
  CTR(t1)=CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,offZ)+CompositeShapeExact::Vec3T::UnitY()*off;
  if(swap) {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  } else {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  }
  off+=dx;
}
void CB2(std::vector<std::shared_ptr<ShapeExact>>& shapes,double offX,double offY,double offZ,bool swap=false) {
  std::shared_ptr<BBoxExact> s1(new BBoxExact(1,1.1,1.2));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(1,.5));
  CompositeShapeExact::Mat3X4T t1;
  CompositeShapeExact::Mat3X4T t2;
  ROT(t1)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<CompositeShapeExact::T>();
  ROT(t2)=ROT(t1)*eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>(0,M_PI/4,M_PI/4),NULL,NULL).template cast<CompositeShapeExact::T>();;
  CTR(t1)=CompositeShapeExact::Vec3T::UnitY()*off;
  CTR(t2)=ROT(t2)*CompositeShapeExact::Vec3T(offX,offY,offZ)+CompositeShapeExact::Vec3T::UnitY()*off;
  if(swap) {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
  } else {
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s1}, {t1})));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({s2}, {t2})));
  }
  off+=dx;
}
void sphereSphere(int argc,char** argv) {
  off=0;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  //sphere-sphere
  SS1(shapes);
  //degenerate sphere-sphere
  SS2(shapes);

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
void sphereCapsule(int argc,char** argv) {
  off=0;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  //sphere-capsule
  SC1(shapes,.4,1.1);
  SC1(shapes,.7,1.1);
  SC1(shapes,-.7,1.1);
  SC1(shapes,-.7,1.1,true);
  //degenerate sphere-capsule
  SC1(shapes,.4,0);
  SC1(shapes,.5,0);
  SC1(shapes,-.5,0);

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
void capsuleCapsule(int argc,char** argv) {
  off=0;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  //colinear capsule-capsule
  CC1(shapes,1.9,.5);
  CC1(shapes,-1.9,.5);
  CC1(shapes,1.5,0);
  CC1(shapes,-1.5,0);
  CC1(shapes,-1.1,0);
  CC1(shapes, 1.1,0);
  //non-colinear capsule-capsule
  CC2(shapes,0.,.3,1.2);
  CC2(shapes,1.7,.3,1.2);
  CC3(shapes,M_PI/10,.4);
  CC3(shapes,-M_PI/10,.4);
  CC3(shapes,M_PI/10,.4,true);
  CC3(shapes,-M_PI/10,.4,true);

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
void sphereBox(int argc,char** argv) {
  off=0;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  //sphere box, interior
  SB1(shapes,.9,0,0);
  SB1(shapes,-.9,0,0);
  SB1(shapes,0,.9,0);
  SB1(shapes,0,-.9,0);
  SB1(shapes,0,0,.9);
  SB1(shapes,0,0,-.9);
  //sphere box, face exterior
  SB1(shapes,1.3,0,0);
  SB1(shapes,-1.3,0,0);
  SB1(shapes,1.1,1.2,0);
  SB1(shapes,-1.1,-1.2,0);
  SB1(shapes,1.1,1.2,1.3);
  SB1(shapes,-1.1,-1.2,-1.3);
  SB1(shapes,-1.1,-1.2,-1.3,true);

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
void capsuleBox(int argc,char** argv) {
  off=0;
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  CB1(shapes, 1.4,0,0);
  CB1(shapes,-1.4,0,0);
  CB1(shapes,0, 1.4,0);
  CB1(shapes,0,-1.4,0);
  CB1(shapes, 1.4,0,0,false);
  CB1(shapes,-1.4,0,0,false);
  CB1(shapes,0, 1.4,0,false);
  CB1(shapes,0,-1.4,0,false);
  CB2(shapes, 3,0,0);
  CB2(shapes,-3,0,0);

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
void boxBox(int argc,char** argv) {
  std::vector<std::shared_ptr<ShapeExact>> shapes;
  for(int i=0; i<10; i++) {
    CompositeShapeExact::Mat3X4T t;
    ROT(t)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>(0,0,i*M_PI/10),NULL,NULL).template cast<CompositeShapeExact::T>();
    CTR(t)=CompositeShapeExact::Vec3T::UnitZ()*i*(2-1e-3);
    std::shared_ptr<BBoxExact> bb(new BBoxExact(1,1,1));
    shapes.push_back(std::shared_ptr<CompositeShapeExact>(new CompositeShapeExact({bb}, {t})));
  }

  std::vector<ContactGenerator::ContactManifold> manifolds;
  ContactGenerator cc(std::shared_ptr<ArticulatedBody>(),shapes);
  cc.generateManifolds(0,false,manifolds,ContactGenerator::Mat3XT(),ContactGenerator::STATIC_STATIC);
  for(const auto& m:manifolds) {
    ASSERT(!m._points.empty())
    for(const auto& p:m._points) {
      ASSERT(p.depth()>0)
    }
  }

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  auto shape=visualizeEnvironment(shapes,true);
  shape->setLineWidth(1);
  drawer.addShape(shape);
  drawer.addShape(visualizeContact(manifolds));
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
}
int main(int argc,char** argv) {
  sphereSphere(argc,argv);
  sphereCapsule(argc,argv);
  capsuleCapsule(argc,argv);
  sphereBox(argc,argv);
  capsuleBox(argc,argv);
  boxBox(argc,argv);
  return 0;
}
