#include "PointCloudExact.h"
#include "EnvironmentUtils.h"
#include <Utils/VTKWriter.h>
#include <Utils/Interp.h>
#include <Utils/Utils.h>
#include <Utils/IO.h>

namespace PHYSICSMOTION {
PointCloudExact::PointCloudExact() {}
PointCloudExact::PointCloudExact(const MeshExact& mesh,T d0) {
  //subdivide
  bool more=true;
  _vss=mesh.vss();
  std::unordered_map<Eigen::Matrix<int,2,1>,int,EdgeHash> midMap;
  std::vector<Eigen::Matrix<int,3,1>> iss=mesh.iss(),issNew;
  while(more) {
    more=false;
    issNew.clear();
    for(auto i:iss) {
      bool needSubd=false;
      for(int d=0; d<3 && !needSubd; d++)
        if((_vss[i[d]]-_vss[i[(d+1)%3]]).squaredNorm()>d0*d0)
          more=needSubd=true;
      if(!needSubd)
        issNew.push_back(i);
      else {
        //subdivide
        Eigen::Matrix<int,3,1> mid;
        for(int d=0; d<3; d++) {
          Eigen::Matrix<int,2,1> id(i[d],i[(d+1)%3]);
          sort2(id[0],id[1]);
          if(midMap.find(id)!=midMap.end())
            mid[d]=midMap[id];
          else {
            mid[d]=midMap[id]=(int)_vss.size();
            _vss.push_back((_vss[id[0]]+_vss[id[1]])/2);
          }
        }
        issNew.push_back(mid);
        issNew.push_back(Eigen::Matrix<int,3,1>(i[0],mid[0],mid[2]));
        issNew.push_back(Eigen::Matrix<int,3,1>(i[1],mid[0],mid[1]));
        issNew.push_back(Eigen::Matrix<int,3,1>(i[2],mid[1],mid[2]));
      }
    }
    iss.swap(issNew);
  }
  //build BVH
  _bvh.clear();
  for(int i=0; i<(int)_vss.size(); i++) {
    _bvh.emplace_back();
    _bvh.back()._l=-1;
    _bvh.back()._r=-1;
    _bvh.back()._cell=i;
    _bvh.back()._parent=-1;
    _bvh.back()._nrCell=1;
    _bvh.back()._bb=BBoxExact(_vss[i].template cast<BBoxExact::T>(),
                              _vss[i].template cast<BBoxExact::T>());
  }
  Node<int,BBoxExact>::buildBVHVertexBottomUp(_bvh,iss,true);
  //distance field
  initSDF(mesh,d0);
}
bool PointCloudExact::read(std::istream& is,IOData* dat) {
  readBinaryData(_vss,is,dat);
  readBinaryData(_bvh,is,dat);
  //distance field
  readBinaryData(_nrNode,is,dat);
  readBinaryData(_stride,is,dat);
  readBinaryData(_bbSDFMinC,is,dat);
  readBinaryData(_bbSDFMaxC,is,dat);
  readBinaryData(_SDF,is,dat);
  return is.good();
}
bool PointCloudExact::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_vss,os,dat);
  writeBinaryData(_bvh,os,dat);
  //distance field
  writeBinaryData(_nrNode,os,dat);
  writeBinaryData(_stride,os,dat);
  writeBinaryData(_bbSDFMinC,os,dat);
  writeBinaryData(_bbSDFMaxC,os,dat);
  writeBinaryData(_SDF,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> PointCloudExact::copy() const {
  std::shared_ptr<PointCloudExact> ret(new PointCloudExact);
  *ret=*this;
  return ret;
}
std::string PointCloudExact::type() const {
  return typeid(PointCloudExact).name();
}
const BBoxExact& PointCloudExact::getBB() const {
  return _bvh.back()._bb;
}
bool PointCloudExact::empty() const {
  return _vss.empty();
}
void PointCloudExact::getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                              std::vector<Eigen::Matrix<int,3,1>>& iss) const {
  ASSERT_MSG(false,"getMesh not supported for PointCloudExact!")
}
bool PointCloudExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                                   T& rad,Eigen::Matrix<int,2,1>& feat,bool cache,
                                   std::vector<Vec3T>* history) const {
  ASSERT_MSG(false,"closestInner not supported for PointCloudExact!")
  return false;
}
const std::vector<Node<int,BBoxExact>>& PointCloudExact::getBVH() const {
  return _bvh;
}
const std::vector<MeshExact::Vec3T>& PointCloudExact::vss() const {
  return _vss;
}
void PointCloudExact::scale(T coef) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v*=coef;
  init(_bvh.size()>1);
}
void PointCloudExact::translate(const Vec3T& pos) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v+=pos;
  init(_bvh.size()>1);
}
void PointCloudExact::transform(const Mat3X4T& trans) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v=ROT(trans)*v+CTR(trans);
  init(_bvh.size()>1);
}
PointCloudExact::Vec3T PointCloudExact::support(const Vec3T& D,int& id) const {
  id=0;
  T dist,maxDist=_vss[0].dot(D);
  for(int i=1; i<(int)_vss.size(); i++) {
    dist=_vss[i].dot(D);
    if(dist>maxDist) {
      maxDist=dist;
      id=i;
    }
  }
  return _vss[id];
}
void PointCloudExact::writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  os.setRelativeIndex();
  for(auto v:_vss)
    vss.push_back((ROT(trans)*v+CTR(trans)).template cast<double>());
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>(0,0,0),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>((int)vss.size(),0,0),
                 VTKWriter<double>::POINT,true);
}
void PointCloudExact::writeSDFVTK(const std::string& path) const {
  BBoxExact bbSDF(_bbSDFMinC.template cast<T>(),_bbSDFMaxC.template cast<T>());
  VTKWriter<double> os("SDF",path,true,bbSDF,(_nrNode.array()-1).matrix(),false);
  os.appendDatas("value",_SDF.begin(),_SDF.end());
}
//helper
void PointCloudExact::init(bool buildBVH) {
  for(int i=0; i<(int)_bvh.size(); i++)
    if(_bvh[i]._cell>=0) {
      _bvh[i]._bb.minCorner()=_vss[i];
      _bvh[i]._bb.maxCorner()=_vss[i];
    } else {
      _bvh[i]._bb=BBoxExact();
      _bvh[i]._bb.setUnion(_bvh[_bvh[i]._l]._bb);
      _bvh[i]._bb.setUnion(_bvh[_bvh[i]._r]._bb);
    }
}
void PointCloudExact::initSDF(const MeshExact& mesh,T d0,int extend) {
  _bbSDFMinC=mesh.getBB().minCorner().template cast<double>();
  _bbSDFMaxC=mesh.getBB().maxCorner().template cast<double>();
  _bbSDFMinC.array()-=(double)d0*extend;
  _bbSDFMaxC.array()+=(double)d0*extend;
  //size of SDF grid
  Eigen::Matrix<double,3,1> ext=_bbSDFMaxC-_bbSDFMinC;
  _nrNode=((ext/(double)d0).array()+1).matrix().template cast<int>();
  _stride=Eigen::Matrix<int,3,1>(1,_nrNode[0],_nrNode[0]*_nrNode[1]);
  Eigen::Matrix<double,3,1> cellSz=(ext.array()/(_nrNode.array().template cast<double>()-1)).matrix();
  _SDF.resize(_nrNode.prod());
  //build grid
  OMP_PARALLEL_FOR_
  for(int z=0; z<_nrNode[2]; z++)
    for(int y=0; y<_nrNode[1]; y++)
      for(int x=0; x<_nrNode[0]; x++) {
        Eigen::Matrix<int,3,1> id(x,y,z);
        Eigen::Matrix<double,3,3> hessian;
        Eigen::Matrix<double,3,1> pt,n,normal;
        Eigen::Matrix<int,2,1> feat;
        pt=(id.template cast<double>().array()*cellSz.array()).matrix()+_bbSDFMinC;
        _SDF[id.dot(_stride)]=mesh.closest<double>(pt,n,normal,hessian,feat);
      }
}
double PointCloudExact::closestSDFInner(const Eigen::Matrix<double,3,1>& pt,Eigen::Matrix<double,3,1>& normal) const {
  Eigen::Matrix<double,3,1> ext=_bbSDFMaxC-_bbSDFMinC;
  Eigen::Matrix<double,3,1> cellSz=(ext.array()/(_nrNode.array().template cast<double>()-1)).matrix();
  Eigen::Matrix<int,3,1> id=((pt-_bbSDFMinC).array()/cellSz.array()).matrix().template cast<int>();
  id=id.cwiseMax(0).cwiseMin(_nrNode-Eigen::Matrix<int,3,1>::Constant(2));
  double val[8]= {_SDF[(id+Eigen::Matrix<int,3,1>(0,0,0)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(1,0,0)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(0,1,0)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(1,1,0)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(0,0,1)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(1,0,1)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(0,1,1)).dot(_stride)],
                  _SDF[(id+Eigen::Matrix<int,3,1>(1,1,1)).dot(_stride)],
                 };
  //interpolate
  Eigen::Matrix<double,3,1> frac=pt-_bbSDFMinC-(id.array().template cast<double>()*cellSz.array()).matrix();
  frac.array()/=cellSz.array();
  double dist=interp3D<double,double>(val[0],val[1],val[2],val[3],
                                      val[4],val[5],val[6],val[7],
                                      frac[0],frac[1],frac[2]);
  auto g=interp3DGrad<double,double>(val[0],val[1],val[2],val[3],
                                     val[4],val[5],val[6],val[7],
                                     frac[0],frac[1],frac[2]);
  normal=Eigen::Matrix<double,3,1>(std::get<0>(g),std::get<1>(g),std::get<2>(g));
  normal.array()/=cellSz.array();
  //std::cout << "pt=" << pt.transpose().template cast<double>()
  //          << " cellSz=" << cellSz.transpose().template cast<double>()
  //          << " minC=" << _bbSDF.minCorner().transpose().template cast<double>()
  //          << " maxC=" << _bbSDF.maxCorner().transpose().template cast<double>()
  //          << " id=" << id.transpose()
  //          << " frac=" << frac.transpose().template cast<double>()
  //          << " normal=" << normal.template cast<double>().norm() << std::endl;
  return dist;
}
}
