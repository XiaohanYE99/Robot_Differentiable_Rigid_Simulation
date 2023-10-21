#include "ConvexDecomposition.h"
#include "MeshExact.h"
#include <Utils/IO.h>
#include <Utils/Pragma.h>

#include <utility>
#define ENABLE_VHACD_IMPLEMENTATION 1
#include "VHACD.h"

namespace PHYSICSMOTION {
ConvexDecomposition::ConvexDecomposition() {}
ConvexDecomposition::ConvexDecomposition(const std::string& inputFile,ConvexDecomposition::T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHulls) {
  MeshExact mesh(inputFile,false);
  _points.resize(mesh.vss().size()*3);
  for(int i=0; i<(int)mesh.vss().size(); i++)
    for(int d=0; d<3; d++)
      _points[i*3+d]=(float)mesh.vss()[i][d];
  _triangles.resize(mesh.iss().size()*3);
  for(int i=0; i<(int)mesh.iss().size(); i++)
    for(int d=0; d<3; d++)
      _triangles[i*3+d]=mesh.iss()[i][d];
  convexDecomposition(scale,pos,rot,maxConvexHulls);
}
void ConvexDecomposition::convexDecomposition(T scale,const Vec3T& pos,const Mat3T& rot,int maxConvexHulls) {
  VHACD::IVHACD::Parameters p;
  p.m_maxConvexHulls=maxConvexHulls;
  p.m_asyncACD=false;
  VHACD::IVHACD* iface=p.m_asyncACD?VHACD::CreateVHACD_ASYNC():VHACD::CreateVHACD();
  iface->Compute(_points.data(),(int)_points.size()/3,_triangles.data(),(int)_triangles.size()/3,p);
  while(!iface->IsReady())
    std::this_thread::sleep_for(std::chrono::nanoseconds(10000)); //s
  if(iface->GetNConvexHulls()) {
    for(uint32_t i=0; i<iface->GetNConvexHulls(); i++) {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      VHACD::IVHACD::ConvexHull ch;
      iface->GetConvexHull(i,ch);
      for(auto& point:ch.m_points)
        vss.emplace_back(point.mX,point.mY,point.mZ);
      for(auto& m_triangle:ch.m_triangles)
        iss.emplace_back(Eigen::Matrix<int,3,1>((int)m_triangle.mI0,(int)m_triangle.mI1,(int)m_triangle.mI2));
      std::shared_ptr<MeshExact> mesh(new MeshExact(vss,iss,true));
      Mat3X4T trans=Mat3X4T::Identity();
      CTR(trans)=pos;
      ROT(trans)=rot;
      mesh->scale((GEOMETRY_SCALAR)scale);
      mesh->transform(trans.template cast<GEOMETRY_SCALAR>());
      _convexHulls.push_back(mesh);
    }
    iface->Release();
  }
}
bool ConvexDecomposition::read(std::istream& is,IOData* dat) {
  registerType<MeshExact>(dat);
  readBinaryData(_points,is);
  readBinaryData(_triangles,is);
  readBinaryData(_convexHulls,is,dat);
  return is.good();
}
bool ConvexDecomposition::write(std::ostream& os,IOData* dat) const {
  registerType<MeshExact>(dat);
  writeBinaryData(_points,os);
  writeBinaryData(_triangles,os);
  writeBinaryData(_convexHulls,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ConvexDecomposition::copy() const {
  return std::shared_ptr<ConvexDecomposition>(new ConvexDecomposition);
}
std::string ConvexDecomposition::type() const {
  return typeid(ConvexDecomposition).name();
}
int ConvexDecomposition::getNumConvexHull() {
  return _convexHulls.size();
}
const std::vector<std::shared_ptr<MeshExact>>& ConvexDecomposition::getConvexHulls() const {
  return _convexHulls;
}
}
