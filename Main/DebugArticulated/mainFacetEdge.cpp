#include <Environment/BBoxExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Utils/RotationUtils.h>
#include <Utils/VTKWriter.h>

using namespace PHYSICSMOTION;

#define RANDOM_TRANSFORM
void box() {
  BBoxExact bb(1,2,3);
  BBoxExact::Mat3X4T t;
#ifdef RANDOM_TRANSFORM
  ROT(t)=expWGradV<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<ShapeExact::T>();
  CTR(t).setRandom();
#else
  t.setIdentity();
#endif
  bb.checkPositiveFacets();
  VTKWriter<double> os("box","box.vtk",true);
  bb.writeVTK(os,t);
  VTKWriter<double> osf("facet","boxFacet.vtk",true);
  bb.writeFacetsVTK(osf,t);
  VTKWriter<double> ose("edge","boxEdge.vtk",true);
  bb.writeEdgesVTK(ose,t);
}
void hull() {
  ConvexHullExact hull({Eigen::Matrix<double,3,1>(-1,-2,-3),
                        Eigen::Matrix<double,3,1>( 1,-2,-3),
                        Eigen::Matrix<double,3,1>( 1, 2,-3),
                        Eigen::Matrix<double,3,1>(-1, 2,-3),

                        Eigen::Matrix<double,3,1>(-1,-2, 3),
                        Eigen::Matrix<double,3,1>( 1,-2, 3),
                        Eigen::Matrix<double,3,1>( 1, 2, 3),
                        Eigen::Matrix<double,3,1>(-1, 2, 3),

                        Eigen::Matrix<double,3,1>(-1.3, 0, 0),
                        Eigen::Matrix<double,3,1>( 1.3, 0, 0),
                        Eigen::Matrix<double,3,1>( 0,-2.3, 0),
                        Eigen::Matrix<double,3,1>( 0, 2.3, 0),
                        Eigen::Matrix<double,3,1>( 0, 0,-3.3),
                        Eigen::Matrix<double,3,1>( 0, 0, 3.3),});
  BBoxExact::Mat3X4T t;
#ifdef RANDOM_TRANSFORM
  ROT(t)=expWGradV<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<ShapeExact::T>();
  CTR(t).setRandom();
#else
  t.setIdentity();
#endif
  hull.checkPositiveFacets();
  VTKWriter<double> os("hull","hull.vtk",true);
  hull.writeVTK(os,t);
  VTKWriter<double> osf("facet","hullFacet.vtk",true);
  hull.writeFacetsVTK(osf,t);
  VTKWriter<double> ose("edge","hullEdge.vtk",true);
  hull.writeEdgesVTK(ose,t);
}
void capsule() {
  SphericalBBoxExact c(1,.5);
  BBoxExact::Mat3X4T t;
#ifdef RANDOM_TRANSFORM
  ROT(t)=expWGradV<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<ShapeExact::T>();
  CTR(t).setRandom();
#else
  t.setIdentity();
#endif
  c.checkPositiveFacets();
  VTKWriter<double> os("box","box.vtk",true);
  c.writeVTK(os,t);
  VTKWriter<double> osf("facet","boxFacet.vtk",true);
  c.writeFacetsVTK(osf,t);
  VTKWriter<double> ose("edge","boxEdge.vtk",true);
  c.writeEdgesVTK(ose,t);
}
int main(int argc,char** argv) {
  box();
  hull();
  capsule();
  return 0;
}
