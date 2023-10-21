#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedVisualizer.h>
#include <TinyVisualizer/ShadowAndLight.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
  std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::createChain(Joint::ROT_3D_XYZ,0.1,10)));
  //std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::createBird(Joint::ROT_3D_XYZ|Joint::TRANS_3D,true));
  //std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::createSpider(Joint::ROT_3D_XYZ|Joint::TRANS_3D));
  //std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::createBipedal(Joint::ROT_3D_XYZ|Joint::TRANS_3D,true));
  //std::shared_ptr<ArticulatedBody> bodyTmp(new ArticulatedBody(ArticulatedLoader::createArm(1,.1));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(*bodyTmp));
  std::shared_ptr<Shape> s=visualizeArticulated(body,Eigen::Matrix<GLfloat,3,1>(.7,.7,.7));
  updateArticulatedBody(s,body,ArticulatedBody::Vec(ArticulatedBody::Vec::Random(body->nrDOF())));
  drawer.addShape(s);

#define USE_LIGHT
#ifdef USE_LIGHT
  drawer.addLightSystem(0);
  drawer.getLight()->lightSz(10);
  for(int x=-1; x<=1; x+=2)
    for(int y=-1; y<=1; y+=2)
      for(int z=-1; z<=1; z+=2)
        drawer.getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(2*x,2*y,2*z),
                                    Eigen::Matrix<GLfloat,3,1>(0,0,0),
                                    Eigen::Matrix<GLfloat,3,1>(.2,.2,.2),
                                    Eigen::Matrix<GLfloat,3,1>(0,0,0));
#endif
  drawer.addCamera3D(90);
  drawer.mainLoop();
  return 0;
}
