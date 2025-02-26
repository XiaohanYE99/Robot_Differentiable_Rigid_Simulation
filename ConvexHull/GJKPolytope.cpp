#include "GJKPolytope.h"
#include "CollisionGradInfo.h"
#include <Utils/VTKWriter.h>
#include <Utils/CrossSpatialUtils.h>
#include <Environment/GJK.h>

namespace PHYSICSMOTION {
//GJKPolytope
template <typename T>
GJKPolytope<T>::GJKPolytope() {
  _jid=-1;
  _jidP=-1;
}
template <typename T>
GJKPolytope<T>::GJKPolytope(int JID,const ArticulatedBody& body,const CollisionGradInfo<T>& info) {
  _jid=JID;
  _jidP=body.joint(JID)._parent;
  _mesh=std::dynamic_pointer_cast<MeshExact>(body.joint(JID)._mesh);
  setVertexId(info._polytopes[JID].getVertexId());
  resetGlobalVss(&info);
}
template <typename T>
GJKPolytope<T>::GJKPolytope(int JID,std::shared_ptr<MeshExact> mesh,const CollisionGradInfo<T>& info){
  _jid=JID;
  _jidP=-1;
  _mesh=mesh;
  if(_jid>=0) resetGlobalVss(&info);
  else resetGlobalVss(NULL);
}
template <typename T>
GJKPolytope<T>::GJKPolytope(std::shared_ptr<MeshExact> mesh) {
  _jid=-1;
  _jidP=-1;
  _mesh=mesh;
  resetGlobalVss(NULL);
}
template <typename T>
bool GJKPolytope<T>::read(std::istream& is,IOData* dat) {
  readBinaryData(_jid,is);
  readBinaryData(_trans,is);
  readBinaryData(_globalVss,is);
  readBinaryData(_mesh,is,dat);
  return is.good();
}
template <typename T>
bool GJKPolytope<T>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_jid,os);
  writeBinaryData(_trans,os);
  writeBinaryData(_globalVss,os);
  writeBinaryData(_mesh,os,dat);
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
void GJKPolytope<T>::setVertexId(const Eigen::Matrix<int,2,1>& off) {
  _offBeg=off[0];
  _offEnd=off[1];
  return;
}
template <typename T>
Eigen::Matrix<int,2,1> GJKPolytope<T>::getVertexId() const {
  return Eigen::Matrix<int,2,1>(_offBeg,_offEnd);
}
template <typename T>
const typename GJKPolytope<T>::Mat3XT& GJKPolytope<T>::globalVss() const {
  return _globalVss;
}
template <typename T>
const typename GJKPolytope<T>::Mat3X4T& GJKPolytope<T>::trans() const {
  return _trans;
}
template <typename T>
std::shared_ptr<MeshExact> GJKPolytope<T>::mesh() const {
  return _mesh;
}
template <typename T>
const std::vector<Node<int,BBoxExact>>& GJKPolytope<T>::getBVH() const {
  return _bvh;
}
template <typename T>
void GJKPolytope<T>::resetGlobalVss(const CollisionGradInfo<T>* info,bool buildBVH) {
  if(info) _trans=TRANSI(info->_info._TM,_jid);
  else _trans.setIdentity();
  _globalVss.resize(3,_mesh->vss().size());
  for(int i=0; i<(int)_mesh->vss().size(); i++)
    _globalVss.col(i)=ROT(_trans)*_mesh->vss()[i].template cast<T>()+CTR(_trans);
  //update bvh
  if(buildBVH){
    _bvh=_mesh->getBVH();
    for(int i=0; i<(int)_bvh.size(); i++){
      /*if(_bvh[i]._cell>=0) {
        Eigen::Matrix<int,3,1> iss=_mesh->iss()[_bvh[i]._cell];
        _bvh[i]._bb=BBoxExact(_globalVss.col(iss[0]).template cast<GEOMETRY_SCALAR>(),
                              _globalVss.col(iss[1]).template cast<GEOMETRY_SCALAR>(),
                              _globalVss.col(iss[2]).template cast<GEOMETRY_SCALAR>());
      } else if(_bvh[i]._cell<0) {
        _bvh[i]._bb=BBoxExact();
        _bvh[i]._bb.setUnion(_bvh[_bvh[i]._l]._bb);
        _bvh[i]._bb.setUnion(_bvh[_bvh[i]._r]._bb);
      }*/
     _bvh[i]._bb.set_center(ROT(_trans)*_bvh[i]._bb.center()+CTR(_trans));
    }
  }
}
template <typename T>
void GJKPolytope<T>::writeVTK(const std::string& path) const {
  VTKWriter<double> os("poly",path,true);
  writeVTK(os);
}
template <typename T>
void GJKPolytope<T>::writeVTK(VTKWriter<double>& os) const {
  ASSERT_MSG(_mesh!=NULL,"Cannot visualize GJKPolytope<T> without _mesh!")
  _mesh->writeVTK(os,_trans.template cast<GEOMETRY_SCALAR>());
}
template <typename T>
void GJKPolytope<T>::writeConfigVTK(const std::string& path,
                                    const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,
                                    GJKPolytope<T>::Point wp1,GJKPolytope<T>::Point wp2) {
  VTKWriter<double> os("GJKConfig",path,true);
  p1.writeVTK(os);
  p2.writeVTK(os);
  //distance
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  vss.push_back(Eigen::Matrix<double,3,1>(wp1[0],wp1[1],wp1[2]));
  vss.push_back(Eigen::Matrix<double,3,1>(wp2[0],wp2[1],wp2[2]));
  iss.push_back(Eigen::Matrix<int,3,1>(0,1,-1));
  os.setRelativeIndex();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::LINE,true);
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
  writeConfigVTK("GJKConfigTest.vtk",p1,p2,pA,pB);
}
template <typename T>
typename GJKPolytope<T>::Vec2T GJKPolytope<T>::project(const GJKPolytope<T>& p,const Vec3T& n) {
  T d;
  Vec2T ret(std::numeric_limits<double>::max(),-std::numeric_limits<double>::max());
  for(int c=0; c<(int)p._mesh->vss().size(); c++) {
    d=(ROT(p._trans)*p._mesh->vss()[c].template cast<T>()+CTR(p._trans)).dot(n);
    ret[0]=std::min<T>(ret[0],d);
    ret[1]=std::max<T>(ret[1],d);
  }
  return ret;
}
template <typename T>
T GJKPolytope<T>::distance(const GJKPolytope<T>& A,const GJKPolytope<T>& B,Point pA,Point pB) {
  MeshExact::Vec3T pAL,pBL;
  MeshExact::T distSqr;
  OMP_CRITICAL_
  distSqr=GJK::runGJK(A._mesh,B._mesh,
                      A._trans.template cast<MeshExact::T>(),
                      B._trans.template cast<MeshExact::T>(),
                      pAL,pBL,NULL);
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
}
//instance
template class GJKPolytope<FLOAT>;
}
