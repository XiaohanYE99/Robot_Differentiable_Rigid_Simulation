#include "GJKPolytope.h"
#include <Utils/VTKWriter.h>
#include <Utils/CrossSpatialUtils.h>
//CGAL's quadratic programming formulation
#ifdef CGAL_SUPPORT
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf                                       ET;
typedef CGAL::Homogeneous<double>                         K;
typedef K::Point_3                                        Point_3;
// ... and the EXACT traits class based on the inexcat kernel
typedef CGAL::Polytope_distance_d_traits_3<K, ET, double> Traits;
typedef CGAL::Polytope_distance_d<Traits>                 Polytope_distance;
#endif
//our GJK library
#include "GJK.h"

//#define USE_CGAL_GJK
#define USE_NATIVE_GJK

namespace PHYSICSMOTION {
//GJKPolytope
template <typename T>
GJKPolytope<T>::GJKPolytope() {
  _jid=-1;
  _n=0;
}
template <typename T>
GJKPolytope<T>::GJKPolytope(int JID,const ArticulatedBody& body,const CollisionGradInfo<T>& info) {
  _jid=JID;
  _n=info._globalVss[JID].cols();
  _trans=TRANSI(info._info._TM,JID);
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  body.joint(JID)._mesh->getMesh(vss,iss);
  _mesh.reset(new MeshExact(vss,iss,false));
  //_mesh=std::dynamic_pointer_cast<MeshExact>(body.joint(JID)._mesh);
  _coord.resize(_n*3);
  for(int i=0; i<_n; i++) {
    _coord[i*3+0]=(double)info._globalVss[JID](0,i);
    _coord[i*3+1]=(double)info._globalVss[JID](1,i);
    _coord[i*3+2]=(double)info._globalVss[JID](2,i);
  }
}
template <typename T>
GJKPolytope<T>::GJKPolytope(std::shared_ptr<MeshExact> mesh) {
  _jid=-1;
  _n=mesh->vss().size();
  _trans.setIdentity();
  _mesh=mesh;
  _coord.resize(_n*3);
  for(int i=0; i<_n; i++) {
    _coord[i*3+0]=(double)mesh->vss()[i][0];
    _coord[i*3+1]=(double)mesh->vss()[i][1];
    _coord[i*3+2]=(double)mesh->vss()[i][2];
  }
}
template <typename T>
bool GJKPolytope<T>::read(std::istream& is,IOData* dat) {
  readBinaryData(_jid,is);
  readBinaryData(_n,is);
  readBinaryData(_trans,is);
  readBinaryData(_mesh,is,dat);
  readBinaryData(_coord,is);
  return is.good();
}
template <typename T>
bool GJKPolytope<T>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_jid,os);
  writeBinaryData(_n,os);
  writeBinaryData(_trans,os);
  writeBinaryData(_mesh,os,dat);
  writeBinaryData(_coord,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> GJKPolytope<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new GJKPolytope<T>);
}
template <typename T>
std::string GJKPolytope<T>::type() const {
  return typeid(GJKPolytope<T>).name();
}
template <typename T>
int GJKPolytope<T>::jid() const {
  return _jid;
}
template <typename T>
std::shared_ptr<MeshExact> GJKPolytope<T>::mesh() const {
  return _mesh;
}
template <typename T>
void GJKPolytope<T>::writeVTK(const std::string& path,const ArticulatedBody* body) const {
  VTKWriter<T> os("poly",path,true);
  writeVTK(os,body);
}
template <typename T>
void GJKPolytope<T>::writeVTK(VTKWriter<T>& os,const ArticulatedBody* body) const {
//  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Vec3T> vss;
  os.setRelativeIndex();
  if(_mesh) {
    for(int i=0; i<(int)_coord.size()/3; i++)
      vss.push_back(Vec3T(_coord[i*3+0],_coord[i*3+1],_coord[i*3+2]));
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(_mesh->iss().begin(),_mesh->iss().end(),VTKWriter<T>::TRIANGLE,true);
  } else {
    for(int i=0; i<(int)_coord.size()/3; i++)
      vss.push_back(Vec3T(_coord[i*3+0],_coord[i*3+1],_coord[i*3+2]));
    os.appendPoints(vss.begin(),vss.end());
    if(body) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body->joint(_jid)._mesh);
      os.appendCells(local->iss().begin(),local->iss().end(),VTKWriter<T>::TRIANGLE,true);
    } else {
      std::vector<Eigen::Matrix<int,3,1>> iss;
      for(int i=0; i<(int)_coord.size()/3; i++)
        iss.push_back(Eigen::Matrix<int,3,1>::Constant(i));
      os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::POINT,true);
    }
  }
}
template <typename T>
void GJKPolytope<T>::writeConfigVTK(const std::string& path,
                                    const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,
                                    GJKPolytope<T>::Point wp1,GJKPolytope<T>::Point wp2,
                                    const ArticulatedBody* body) {
  VTKWriter<T> os("GJKConfig",path,true);
  p1.writeVTK(os,body);
  p2.writeVTK(os,body);
  //distance
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  vss.push_back(Eigen::Matrix<double,3,1>(wp1[0],wp1[1],wp1[2]));
  vss.push_back(Eigen::Matrix<double,3,1>(wp2[0],wp2[1],wp2[2]));
  iss.push_back(Eigen::Matrix<int,3,1>(0,1,-1));
  os.setRelativeIndex();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
template <typename T>
void GJKPolytope<T>::writeConfig(const std::string& path,
                                 const GJKPolytope<T>& p1,const GJKPolytope<T>& p2) {
  std::ofstream os(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  registerType<MeshExact>(dat.get());
  p1.write(os,dat.get());
  p2.write(os,dat.get());
}
template <typename T>
void GJKPolytope<T>::readConfig(const std::string& path) {
  std::ifstream is(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  registerType<MeshExact>(dat.get());
  GJKPolytope<T> p1,p2;
  p1.read(is,dat.get());
  p2.read(is,dat.get());
  Point pA,pB;
  distance(p1,p2,pA,pB);
  writeConfigVTK("GJKConfigTest.vtk",p1,p2,pA,pB,NULL);
}
template <typename T>
typename GJKPolytope<T>::Vec2T GJKPolytope<T>::project(const GJKPolytope<T>& p,const Vec3T& n) {
  T d;
  Vec2T ret(std::numeric_limits<double>::max(),-std::numeric_limits<double>::max());
  for(int c=0; c<(int)p._coord.size(); c+=3) {
    d=p._coord[c+0]*n[0]+p._coord[c+1]*n[1]+p._coord[c+2]*n[2];
    ret[0]=std::min<T>(ret[0],d);
    ret[1]=std::max<T>(ret[1],d);
  }
  return ret;
}
template <typename T>
T GJKPolytope<T>::distance(const GJKPolytope<T>& A,const GJKPolytope<T>& B,Point pA,Point pB) {
#ifdef USE_CGAL_GJK
  std::vector<Point_3> pAss(A._n),pBss(B._n);
  for(int i=0,j=0; i<A._n; i++,j+=3)
    pAss[i]=Point_3(A._coord[j+0],A._coord[j+1],A._coord[j+2]);
  for(int i=0,j=0; i<B._n; i++,j+=3)
    pBss[i]=Point_3(B._coord[j+0],B._coord[j+1],B._coord[j+2]);
  Polytope_distance pd(pAss.data(),pAss.data()+A._n,pBss.data(),pBss.data()+B._n);
  Polytope_distance::Coordinate_iterator coord_it;
  if(!pd.is_valid())
    return 0;
  if(pA) {
    int it=0;
    double pAH[4];
    for(coord_it=pd.realizing_point_p_coordinates_begin(); coord_it!=pd.realizing_point_p_coordinates_end(); ++coord_it)
      pAH[it++]=CGAL::to_double(*coord_it);
    for(it=0; it<3; it++)
      pA[it]=pAH[it]/pAH[3];
  }
  if(pB) {
    int it=0;
    double pBH[4];
    for(coord_it=pd.realizing_point_q_coordinates_begin(); coord_it!=pd.realizing_point_q_coordinates_end(); ++coord_it)
      pBH[it++]=CGAL::to_double(*coord_it);
    for(it=0; it<3; it++)
      pB[it]=pBH[it]/pBH[3];
  }
  T distSqr=CGAL::to_double(pd.squared_distance_numerator())/CGAL::to_double(pd.squared_distance_denominator());
  return sqrt(distSqr);
#else
  MeshExact::Vec3T pAL,pBL;
  MeshExact::T distSqr;
  OMP_CRITICAL_
  distSqr=GJK::runGJK(A._mesh,B._mesh,
                      A._trans.template cast<MeshExact::T>(),
                      B._trans.template cast<MeshExact::T>(),
                      pAL,pBL);
  if(pA) {
    Vec3T pAG=ROT(A._trans)*pAL.template cast<T>()+CTR(A._trans);
    for(int i=0; i<3; i++)
      pA[i]=(double)pAG[i];
  }
  if(pB) {
    Vec3T pBG=ROT(B._trans)*pBL.template cast<T>()+CTR(B._trans);
    for(int i=0; i<3; i++)
      pB[i]=(double)pBG[i];
  }
  return sqrt((T)distSqr);
#endif
}
//instance
template class GJKPolytope<FLOAT>;
}
