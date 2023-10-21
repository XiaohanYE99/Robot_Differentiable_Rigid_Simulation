#include <Articulated/ArticulatedLoader.h>
#include <Articulated/NEArticulatedGradientInfo.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Articulated/ArticulatedUtils.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  mpfr_float::default_precision(100);
  ArticulatedBody body=ArticulatedLoader::createChain(Joint::TRANS_3D|Joint::ROT_3D_XYZ,0.1,10);
  NEArticulatedGradientInfo<T>::debug(body);
  PBDArticulatedGradientInfo<T>::debug(body);
  //simplify
  ArticulatedUtils(body).simplify([&](int,const Joint& J) {
    std::set<std::string> names({"joint","joint_1_2","joint_1_2_3_4","joint_1_2_3_4_5_6","joint_1_2_3_4_5_6_7_8"});
    return names.find(J._name)!=names.end();
  },ArticulatedBody::Vec::Random(body.nrDOF()),10);
  //scale mass
  std::cout << "mass=" << ArticulatedUtils(body).totalMass();
  ArticulatedUtils(body).scaleMass(2);
  std::cout << " mass*2=" << ArticulatedUtils(body).totalMass() << std::endl;
  //transform mesh
  body.joint(0).debugTransformMesh();
  return 0;
}
