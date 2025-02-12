#include <Environment/GJK.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/EnvironmentVisualizer.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>
#include <Utils/RotationUtils.h>
#include <random>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  std::shared_ptr<SphericalBBoxExact> s1(new SphericalBBoxExact(1,2,3,.1));
  std::shared_ptr<SphericalBBoxExact> s2(new SphericalBBoxExact(.5,.5,.5,.1));
  GJK::Vec3T r1=GJK::Vec3T::Identity(),t1=GJK::Vec3T::Zero();
  GJK::Vec3T r2=GJK::Vec3T::Identity(),t2=GJK::Vec3T::Zero();
  GJK::T bary=1,speed=0.01,rng=5;

  using namespace DRAWER;
  Drawer drawer(argc,argv);
  std::shared_ptr<Bullet3DShape> shape1(new Bullet3DShape);
  std::shared_ptr<Bullet3DShape> shape2(new Bullet3DShape);
  std::shared_ptr<MeshShape> contact(new MeshShape);
  shape1->addShape(visualizeShapeExact(s1,true));
  shape2->addShape(visualizeShapeExact(s2,true));
  shape1->setLineWidth(1);
  shape2->setLineWidth(1);
  drawer.addShape(shape1);
  drawer.addShape(shape2);
  drawer.addShape(contact);

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
        r2=GJK::Vec3T::Random()*M_PI;
        t2=GJK::Vec3T::Random()*rng;
      }
    }
    //set transform
    auto r=expWGradV<GLfloat,Eigen::Matrix<GLfloat,3,1>>((r1*(1-bary)+r2*bary).template cast<GLfloat>(),NULL,NULL);
    auto t=(t1*(1-bary)+t2*bary).template cast<GLfloat>();
    shape2->setLocalRotate(r);
    shape2->setLocalTranslate(t);
    //GJK
    GJK::Vec3T pAL,pBL;
    GJK::Mat3X4T T1=GJK::Mat3X4T::Identity(),T2;
    ROT(T2)=r.template cast<GJK::T>();
    CTR(T2)=t.template cast<GJK::T>();
    GJK::runGJK(s1,s2,T1,T2,pAL,pBL,&intersect);
    //draw GJK result
    contact->clear();
    contact->addVertex(pAL.template cast<GLfloat>());
    contact->addVertex(r*pBL.template cast<GLfloat>()+t);
    contact->addIndexSingle(0);
    contact->addIndexSingle(1);
    contact->setMode(GL_LINES);
    if(intersect)
      contact->setColorDiffuse(GL_LINES,1,0,0);
    else contact->setColorDiffuse(GL_LINES,0,1,0);
    contact->setLineWidth(5);
  });
  drawer.mainLoop();
  return 0;
}
