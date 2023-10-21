#include <Articulated/ArticulatedLoader.h>
#include <Articulated/PBCentroidBodyDynamicsGradientInfo.h>
#include <Environment/DeformedEnvironment.h>
#include <Utils/Utils.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  mpfr_float::default_precision(100);
  ArticulatedBody body=ArticulatedLoader::createChain(Joint::TRANS_3D|Joint::ROT_3D_EXP,0.1,10);
  PBCentroidBodyDynamicsGradientInfo<T>::debug(body,4,PBCentroidBodyDynamicsGradientInfo<T>::Vec3T::Random(),0.01);
  PBCentroidBodyDynamicsGradientInfo<T>::debugContinuous(body,4);
  if(exists("DEnv/binary.dat")) {
    std::shared_ptr<DeformedEnvironment<T>> DEnv(new DeformedEnvironment<T>);
    DEnv->readStr("DEnv/binary.dat");
    PBCentroidBodyDynamicsGradientInfo<T>::debug(body,4,PBCentroidBodyDynamicsGradientInfo<T>::Vec3T::Random(),0.01,DEnv);
  }
  return 0;
}
