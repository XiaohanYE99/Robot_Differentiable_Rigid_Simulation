#include "ShapeExact.h"
#include <Utils/VTKWriter.h>

namespace PHYSICSMOTION {
//for GJK
ShapeExact::Vec3T ShapeExact::support(const Vec3T& D,int& id) const {
  FUNCTION_NOT_IMPLEMENTED
  return Vec3T();
}
//for SAT
ShapeExact::Vec2T ShapeExact::project(const Vec3T& d) const {
  FUNCTION_NOT_IMPLEMENTED
  return Vec2T();
}
std::vector<ShapeExact::Facet> ShapeExact::facets() const {
  FUNCTION_NOT_IMPLEMENTED
  return std::vector<Facet>();
}
std::vector<ShapeExact::Edge> ShapeExact::edges() const {
  FUNCTION_NOT_IMPLEMENTED
  return std::vector<Edge>();
}
void ShapeExact::checkPositiveFacets() const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> tss;
  for(const auto& f:facets()) {
    //boundary
    Vec3T ctr=Vec3T::Zero();
    for(const auto& v:f._boundary)
      ctr+=v;
    ctr/=(int)f._boundary.size();
    //index
    for(int i=0; i<(int)f._boundary.size(); i++) {
      auto d1=f._boundary[i]-ctr;
      auto d2=f._boundary[(i+1)%(int)f._boundary.size()]-ctr;
      T area=d1.cross(d2).dot(f._n);
      ASSERT_MSG(area>0,"Facet area is negative, this is impossible!")
    }
  }
}
void ShapeExact::writeFacetsVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  std::vector<Eigen::Matrix<double,3,1>> vss,css;
  std::vector<Eigen::Matrix<int,3,1>> tss;
  std::vector<Eigen::Matrix<int,2,1>> lss;
  for(const auto& f:facets()) {
    //boundary
    int id=(int)vss.size();
    Vec3T ctr=Vec3T::Zero();
    Eigen::Matrix<double,3,1> c=(Eigen::Matrix<double,3,1>::Random().array()*0.5+0.5).matrix();
    for(const auto& v:f._boundary) {
      vss.push_back((ROT(trans)*v+CTR(trans)).template cast<double>());
      css.push_back(c);
      ctr+=v;
    }
    ctr/=(int)f._boundary.size();
    //center,normal
    int id2=(int)vss.size();
    vss.push_back((ROT(trans)*ctr+CTR(trans)).template cast<double>());
    css.push_back(c);
    vss.push_back((ROT(trans)*(ctr+f._n)+CTR(trans)).template cast<double>());
    css.push_back(c);
    //index
    for(int i=0; i<(int)f._boundary.size(); i++)
      tss.push_back({id2,id+i,id+(i+1)%(int)f._boundary.size()});
    lss.push_back({id2,id2+1});
  }
  if(vss.empty())
    return;
  os.setRelativeIndex();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(lss.begin(),lss.end(),VTKWriter<double>::LINE,true);
  os.appendCells(tss.begin(),tss.end(),VTKWriter<double>::TRIANGLE,true);
  os.appendCustomPointColorData("color",css.begin(),css.end());
}
void ShapeExact::writeEdgesVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  std::vector<Eigen::Matrix<double,3,1>> vss,css;
  std::vector<Eigen::Matrix<int,2,1>> ess;
  int id=0;
  for(const auto& e:edges()) {
    vss.push_back((ROT(trans)*e._a+CTR(trans)).template cast<double>());
    vss.push_back((ROT(trans)*e._b+CTR(trans)).template cast<double>());
    css.push_back((Eigen::Matrix<double,3,1>::Random().array()*0.5+0.5).matrix());
    css.push_back(css.back());
    ess.push_back({id*2+0,id*2+1});
    id++;
  }
  if(vss.empty())
    return;
  os.setRelativeIndex();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(ess.begin(),ess.end(),VTKWriter<double>::LINE,true);
  os.appendCustomPointColorData("color",css.begin(),css.end());
}
}
