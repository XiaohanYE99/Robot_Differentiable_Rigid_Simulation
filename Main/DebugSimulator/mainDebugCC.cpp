#include "ConvexHull/ConvexHullDistanceEnergy.h"
#include "ConvexHull/ConvexHullDistanceConvexEnergy.h"
#include "ConvexHull/ConvexHullMeshDistanceEnergy.h"
#include "ConvexHull/ConvexHullDistanceFrictionEnergy.h"
#include "ConvexHull/ConvexHullMeshDistanceFrictionEnergy.h"
#include "Articulated/ArticulatedUtils.h"
#include "Articulated/ArticulatedLoader.h"
#include "Environment/EnvironmentUtils.h"
#include "Environment/ConvexHullExact.h"

#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/MeshBasedPBDSimulator.h>

using namespace PHYSICSMOTION;
template <typename T,typename CCEnergyType>
void debugCCBarrierConvexEnergy() {
  std::shared_ptr<MeshExact> mesh;
  {
    std::vector<Eigen::Matrix<double, 3, 1>> vss;
    std::vector<Eigen::Matrix<int, 3, 1>> iss;
    addBox(vss,iss,Eigen::Matrix<double, 3, 1>(-30,-30,-200),Eigen::Matrix<double, 3, 1>(20,20,-.08));
    mesh.reset(new MeshExact(vss,iss));
  }
  std::vector<Eigen::Matrix<double,3,1>> shape;
  shape.push_back(Eigen::Matrix<double,3,1>(.3,.3,.3));
  shape.push_back(Eigen::Matrix<double,3,1>(.1,.1,.8));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::createPushTask(shape,true)));
  //std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",true,false)));
  ArticulatedUtils(*body).tessellate(true);
  //ArticulatedUtils(*body).BBApproxiate(true);
  //ArticulatedUtils(*body).makeConvex();
  Px barrier;
  barrier._x0=.1;
  T d0=0;
  GJKPolytope<T> p(mesh);
  //CCEnergyType::debugGradient(p,*body,4,barrier._x0,d0);
  CCEnergyType::debugGradient(*body,2,4,barrier._x0,d0);
}
template <typename T,typename CCEnergyType>
void debugCCBarrierEnergy() {
  std::shared_ptr<MeshExact> mesh;
  {
    std::vector<Eigen::Matrix<double, 3, 1>> vss;
    std::vector<Eigen::Matrix<int, 3, 1>> iss;
    addBox(vss,iss,Eigen::Matrix<double, 3, 1>(0,0,0),Eigen::Matrix<double, 3, 1>(.2,.2,.2));
    mesh.reset(new MeshExact(vss,iss));
  }

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",true,false)));
  ArticulatedUtils(*body).tessellate(true);
  ArticulatedUtils(*body).BBApproxiate(true);

  Px barrier;
  barrier._x0=1.0;
  T d0=1e-3;
  GJKPolytope<T> p(mesh);
  CCEnergyType::debugGradient(true,p,*body,7,barrier._x0,d0);
  CCEnergyType::debugGradient(false,p,*body,7,barrier._x0,d0);
  CCEnergyType::debugGradient(true,*body,4,7,barrier._x0,d0);
  CCEnergyType::debugGradient(false,*body,4,7,barrier._x0,d0);
}
template <typename T,typename CCEnergyType>
void debugCCBarrierConvexFrictionEnergy() {
  std::shared_ptr<MeshExact> mesh;
  {
    std::vector<Eigen::Matrix<double, 3, 1>> vss;
    std::vector<Eigen::Matrix<int, 3, 1>> iss;
    addBox(vss,iss,Eigen::Matrix<double, 3, 1>(0,0,0),Eigen::Matrix<double, 3, 1>(.2,.2,.2));
    mesh.reset(new ConvexHullExact(vss));
  }

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",true,false)));
  ArticulatedUtils(*body).tessellate(true);
  //ArticulatedUtils(*body).BBApproxiate(true);
  //ArticulatedUtils(*body).makeConvex();

  Px barrier;
  barrier._x0=1.0;
  T d0=1e-3;
  GJKPolytope<T> p(mesh);
  CCEnergyType::debugGradient(p,*body,7,barrier._x0,d0);
  CCEnergyType::debugGradient(*body,4,7,barrier._x0,d0);
}
template <typename T,typename CCEnergyType>
void debugCCBarrierFrictionEnergy() {
  std::shared_ptr<MeshExact> mesh;
  {
    std::vector<Eigen::Matrix<double, 3, 1>> vss;
    std::vector<Eigen::Matrix<int, 3, 1>> iss;
    addBox(vss,iss,Eigen::Matrix<double, 3, 1>(0,0,0),Eigen::Matrix<double, 3, 1>(.2,.2,.2));
    mesh.reset(new MeshExact(vss,iss));
  }

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::readURDF("../data/kuka_lwr/kuka.urdf",true,false)));
  ArticulatedUtils(*body).tessellate(true);
  ArticulatedUtils(*body).BBApproxiate(true);

  Px barrier;
  barrier._x0=0.2;
  T d0=1e-3;
  GJKPolytope<T> p(mesh);
  CCEnergyType::debugGradient(true,p,*body,7,barrier._x0,d0);
  CCEnergyType::debugGradient(false,p,*body,7,barrier._x0,d0);
  CCEnergyType::debugGradient(true,*body,4,7,barrier._x0,d0);
  CCEnergyType::debugGradient(false,*body,4,7,barrier._x0,d0);
}
int main() {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  //debugCCBarrierEnergy<T,CCBarrierEnergy<T,Px>>();
  //debugCCBarrierEnergy<T,CCBarrierConvexEnergy<T,Px>>();
  debugCCBarrierConvexEnergy<T,CCBarrierMeshEnergy<T,Logx>>();
  //debugCCBarrierFrictionEnergy<T,CCBarrierFrictionEnergy<T,Px>>();
  //debugCCBarrierConvexFrictionEnergy<T,CCBarrierMeshFrictionEnergy<T,Px>>();
  return 0;
}
