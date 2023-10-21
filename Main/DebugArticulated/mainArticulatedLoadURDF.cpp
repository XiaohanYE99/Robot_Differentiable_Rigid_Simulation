#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Environment/MeshExact.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/TrackballCameraManipulator.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
  //compareVisualMeshURDF(argc,argv,"data/ANYmal/robot_large_limits.urdf");
//  std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",false,true)));
  std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("../data/dragon_ddk/dddk.urdf",false,true)));
  // std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("data/ANYmal/robot_large_limits.urdf",false,true)));
  // std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("../data/BlackBird/urdf/blackbird.urdf",false,true)));
  //std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::readURDF("data/LittleDog/LittleDog.urdf",false,true)));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(*bodyTmp));
  ArticulatedUtils(*body).simplify([&](int,const Joint& J) {
    std::set<std::string> names({"imu_link"});
    return names.find(J._name)!=names.end();
  },ArticulatedUtils::Vec::Zero(body->nrDOF()),10);
  ArticulatedUtils(*body).addBase(3,ArticulatedUtils::Vec3T::Zero());
  ArticulatedUtils(*body).simplify(ArticulatedUtils::Vec::Zero(body->nrDOF()),10);
  std::shared_ptr<Shape> s=visualizeArticulated(body,Eigen::Matrix<GLfloat,3,1>(.7,.7,.7));
  s->setColorAmbient(GL_TRIANGLES,.5,.5,.5);
  s->setColorSpecular(GL_TRIANGLES,1,1,1);
  s->setColor(GL_TRIANGLES,1,1,1);
  updateArticulatedBody(s,body,ArticulatedBody::Vec(ArticulatedBody::Vec::Zero(body->nrDOF())));
  drawer.addShape(s);

#define USE_LIGHT
#ifdef USE_LIGHT
  GLfloat lightRng=0.25f;
  drawer.addLightSystem();
  drawer.getLight()->lightSz(10);
  for(int x=-1; x<=1; x+=2)
    for(int y=-1; y<=1; y+=2)
      for(int z=-1; z<=1; z+=2)
        drawer.getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(lightRng*x,lightRng*y,lightRng*z),
                                    Eigen::Matrix<GLfloat,3,1>(0,0,0),
                                    Eigen::Matrix<GLfloat,3,1>(.2,.2,.2),
                                    Eigen::Matrix<GLfloat,3,1>(.5,.5,.5));
#endif
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,0,1));
  //drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new TrackballCameraManipulator(drawer.getCamera3D())));
  drawer.mainLoop();
  return 0;
}
