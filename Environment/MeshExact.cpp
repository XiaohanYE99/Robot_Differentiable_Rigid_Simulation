#include "MeshExact.h"
#include "EnvironmentUtils.h"
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include <Utils/IO.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <stack>

namespace PHYSICSMOTION {
MeshExact::MeshExact() {}
MeshExact::MeshExact(const std::string& path,bool buildBVH) {
  Assimp::Importer importer;
  const aiScene* scene=importer.ReadFile(path.c_str(),aiProcess_JoinIdenticalVertices);
  ASSERT_MSGV(scene,"Mesh %s is empty!",path.c_str())
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  init(scene,NULL,vss,iss,buildBVH);
}
MeshExact::MeshExact(const aiScene* scene,bool buildBVH) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  init(scene,NULL,vss,iss,buildBVH);
}
MeshExact::MeshExact(std::vector<Eigen::Matrix<double,3,1>>& vss,
                     std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH) {
  if(buildBVH) {
    makeUniform(iss);
    if(volume(vss,iss)<0)
      makeInsideOut(iss);
  }
  init(vss,iss,buildBVH);
}
MeshExact::MeshExact(std::vector<Eigen::Matrix<double,3,1>>& vss,
                     std::vector<Eigen::Matrix<double,2,1>>& tcss,
                     std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH) {
  if(buildBVH) {
    makeUniform(iss);
    if(volume(vss,iss)<0)
      makeInsideOut(iss);
  }
  init(vss,tcss,iss,buildBVH);
}
void MeshExact::init(const aiScene* scene,const aiNode* node,
                     std::vector<Eigen::Matrix<double,3,1>>& vss,
                     std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH) {
  if(!node) {
    init(scene,scene->mRootNode,vss,iss);
    if(buildBVH) {
      makeUniform(iss);
      if(volume(vss,iss)<0)
        makeInsideOut(iss);
    }
    init<double>(vss,iss,buildBVH);
  } else {
    //add children
    for(int i=0; i<(int)node->mNumChildren; i++)
      init(scene,node->mChildren[i],vss,iss);
    //add mesh
    for(int idm=0; idm<(int)node->mNumMeshes; idm++) {
      const aiMesh* m=scene->mMeshes[node->mMeshes[idm]];
      int off=(int)vss.size();
      for(int i=0; i<(int)m->mNumVertices; i++) {
        const aiVector3D& v=m->mVertices[i];
        vss.push_back(Eigen::Matrix<double,3,1>(v.x,v.y,v.z));
      }
      for(int i=0; i<(int)m->mNumFaces; i++) {
        const aiFace& f=m->mFaces[i];
        for(int j=0; j<(int)f.mNumIndices-2; j++)
          iss.push_back(Eigen::Matrix<int,3,1>(off+f.mIndices[0],off+f.mIndices[j+1],off+f.mIndices[j+2]));
      }
    }
    Eigen::Matrix<double,4,4> T;
    //row0
    T(0,0)=node->mTransformation.a1;
    T(0,1)=node->mTransformation.a2;
    T(0,2)=node->mTransformation.a3;
    T(0,3)=node->mTransformation.a4;
    //row1
    T(1,0)=node->mTransformation.b1;
    T(1,1)=node->mTransformation.b2;
    T(1,2)=node->mTransformation.b3;
    T(1,3)=node->mTransformation.b4;
    //row2
    T(2,0)=node->mTransformation.c1;
    T(2,1)=node->mTransformation.c2;
    T(2,2)=node->mTransformation.c3;
    T(2,3)=node->mTransformation.c4;
    //row3
    T(3,0)=node->mTransformation.d1;
    T(3,1)=node->mTransformation.d2;
    T(3,2)=node->mTransformation.d3;
    T(3,3)=node->mTransformation.d4;
    for(int i=0; i<(int)vss.size(); i++) {
      Eigen::Matrix<double,4,1> v=T*Eigen::Matrix<double,4,1>(vss[i][0],vss[i][1],vss[i][2],1);
      vss[i]=v.template segment<3>(0)/v[3];
    }
  }
}
template <typename T2>
void MeshExact::init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
                     const std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH) {
  //vss
  if(vss.empty() || iss.empty())
    return;
  _vss.resize(vss.size());
  for(int i=0; i<(int)vss.size(); i++)
    _vss[i]=vss[i].template cast<T>();
  //iss
  _iss=iss;
  //bss
  initBss();
  if(buildBVH) {
    //tss
    _tss.resize(_iss.size());
    for(int i=0; i<(int)_tss.size(); i++)
      _tss[i]=TriangleExact(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]);
    //bvh
    _bvh.assign(_iss.size(),Node<int,BBoxExact>());
    for(int i=0; i<(int)_bvh.size(); i++) {
      Node<int,BBoxExact>& n=_bvh[i];
      n._bb=BBoxExact(_tss[i].getBB());
      n._nrCell=1;
      n._cell=i;
    }
    Node<int,BBoxExact>::buildBVHTriangleBottomUp(_bvh,iss,true);
  } else {
    //tss
    _tss.clear();
    //bvh
    _bvh.resize(1);
    ASSERT_MSG(_vss.size()>=3,"We do not accept mesh having <= 3 vertices!")
    _bvh[0]._bb=BBoxExact(_vss[0],_vss[1],_vss[2]);
    for(int i=3; i<(int)_vss.size(); i++)
      _bvh[0]._bb.setUnion(_vss[i]);
  }
}
template <typename T2>
void MeshExact::init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
                     const std::vector<Eigen::Matrix<T2,2,1>>& tcss,
                     const std::vector<Eigen::Matrix<int,3,1>>& iss,bool buildBVH) {
  init(vss,iss,buildBVH);
  for(int i=0; i<(int)tcss.size(); i++)
    _tcss.push_back(tcss[i].template cast<T>());
  ASSERT_MSGV(tcss.size()==vss.size() || tcss.empty(),"Texture coordinate array size(%d)!=number of vertices(%d)!",(int)tcss.size(),(int)vss.size())
}
bool MeshExact::read(std::istream& is,IOData* dat) {
  readBinaryData(_bss,is);
  readBinaryData(_vss,is,dat);
  readBinaryData(_tcss,is,dat);
  readBinaryData(_iss,is,dat);
  readBinaryData(_tss,is,dat);
  readBinaryData(_bvh,is,dat);
  return is.good();
}
bool MeshExact::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_bss,os);
  writeBinaryData(_vss,os,dat);
  writeBinaryData(_tcss,os,dat);
  writeBinaryData(_iss,os,dat);
  writeBinaryData(_tss,os,dat);
  writeBinaryData(_bvh,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> MeshExact::copy() const {
  std::shared_ptr<MeshExact> ret(new MeshExact);
  *ret=*this;
  return ret;
}
std::string MeshExact::type() const {
  return typeid(MeshExact).name();
}
const BBoxExact& MeshExact::getBB() const {
  return _bvh.back()._bb;
}
bool MeshExact::empty() const {
  return _vss.empty();
}
const std::vector<Node<int,BBoxExact>>& MeshExact::getBVH() const {
  return _bvh;
}
const std::vector<char>& MeshExact::bss() const {
  return _bss;
}
const std::vector<MeshExact::Vec3T>& MeshExact::vss() const {
  return _vss;
}
const std::vector<MeshExact::Vec2T>& MeshExact::tcss() const {
  return _tcss;
}
const std::vector<Eigen::Matrix<int,3,1>>& MeshExact::iss() const {
  return _iss;
}
void MeshExact::getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                        std::vector<Eigen::Matrix<int,3,1>>& iss) const {
  vss.clear();
  for(const Vec3T& v:_vss)
    vss.push_back(v.template cast<double>());
  iss=_iss;
}
bool MeshExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                             T&,Eigen::Matrix<int,2,1>& feat,bool,
                             std::vector<Vec3T>*) const {
  Eigen::Matrix<int,2,1> featTmp;
  Vec3T cp,cpTmp,bTmp;
  feat=Eigen::Matrix<int,2,1>(-1,-1);
  T distSqr=std::numeric_limits<double>::max(),distSqrTmp;
  cp.setConstant(distSqr);
  //main loop
  std::stack<std::pair<T,int>> ss;
  ss.push(std::make_pair(_bvh.back()._bb.distToSqr(pt),(int)_bvh.size()-1));
  while(!ss.empty()) {
    T distSqrToBB=ss.top().first;
    const Node<int,BBoxExact>& node=_bvh[ss.top().second];
    ss.pop();
    if(distSqrToBB>distSqr)
      continue;
    else if(node._cell>=0) {
      const TriangleExact& te=_tss[node._cell];
      te.calcPointDist(pt,distSqrTmp,cpTmp,bTmp,featTmp);
      //get feature id
      for(int d=0; d<2; d++)
        if(featTmp[d]>=0)
          featTmp[d]=_iss[node._cell][featTmp[d]];
      //make feature id Unique
      if(featTmp[1]>=0 && featTmp[0]>featTmp[1])
        std::swap(featTmp[0],featTmp[1]);
      else if(featTmp[0]==-1)
        featTmp[1]=node._cell;
      //update distance
      if(featTmp[0]>=0 && feat==featTmp) {
        //update due to same feature
#ifdef USE_RATIONAL
        ASSERT(distSqr==distSqrTmp)
#endif
        T align=abs((cp-pt).dot(normal));
        T alignCurr=abs((cpTmp-pt).dot(te.normal()));
        if(alignCurr>align) {
          distSqr=distSqrTmp;
          normal=te.normal();
          cp=cpTmp;
          n=cp-pt;
        }
      } else if(distSqrTmp<distSqr) {
        //update due to new distance
        //if((cp-pt).dot(te.normal())==0)
        //  continue;
        feat=featTmp;
        distSqr=distSqrTmp;
        normal=te.normal();
        cp=cpTmp;
        n=cp-pt;
      }
    } else {
      const Node<int,BBoxExact>& nl=_bvh[node._l];
      const Node<int,BBoxExact>& nr=_bvh[node._r];
      T distSqrToBBL=nl._bb.distToSqr(pt);
      T distSqrToBBR=nr._bb.distToSqr(pt);
      if(distSqrToBBL<distSqrToBBR) {
        ss.push(std::make_pair(distSqrToBBR,node._r));
        ss.push(std::make_pair(distSqrToBBL,node._l));
      } else {
        ss.push(std::make_pair(distSqrToBBL,node._l));
        ss.push(std::make_pair(distSqrToBBR,node._r));
      }
    }
  }
  //adjust hessian
  if(distSqr==0) {
    //surface
    hessian.setZero();
  } else if(feat[0]==-1) {
    //surface
    hessian.setZero();
  } else if(feat[1]==-1) {
    //vertex
    ASSERT(feat[0]>=0)
    hessian=n*n.transpose();
    hessian/=hessian.trace();
    hessian-=Mat3T::Identity();
  } else {
    //edge
    ASSERT(feat[0]>=0 && feat[1]>=0)
    Vec3T e=_vss[feat[0]]-_vss[feat[1]];
    T eDotE=e.dot(e);
    Vec3T d=n-n.dot(e)*e/eDotE;
    T dDotD=d.dot(d);
    Vec3T nd=e.dot(d)*e/eDotE/dDotD;
    hessian=e*e.transpose()/eDotE-Mat3T::Identity();
    hessian+=d*d.transpose()/dDotD;
    hessian-=d*nd.transpose();
  }
  return n.dot(normal)>0;
}
void MeshExact::scale(T coef) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v*=coef;
  bool buildBVH=_bvh.size()>1;
  init(vss,_iss,buildBVH);
}
void MeshExact::translate(const Vec3T& pos) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v+=pos;
  bool buildBVH=_bvh.size()>1;
  init(vss,_iss,buildBVH);
}
void MeshExact::transform(const Mat3X4T& trans) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v=ROT(trans)*v+CTR(trans);
  bool buildBVH=_bvh.size()>1;
  init(vss,_iss,buildBVH);
}
MeshExact MeshExact::addMesh(const MeshExact& mesh) {
  std::vector<Vec3T> vss=_vss;
  std::vector<Eigen::Matrix<int,3,1>> iss=_iss;
  std::vector<Eigen::Matrix<int,3,1>> meshIss=mesh.iss();
  std::vector<Node<int,BBoxExact>> bvh=_bvh;
  for(int i=0; i<(int)meshIss.size(); ++i)
    meshIss[i]+=Eigen::Matrix<int,3,1>::Constant(vss.size());
  vss.insert(vss.end(),mesh.vss().begin(),mesh.vss().end());
  iss.insert(iss.end(),meshIss.begin(),meshIss.end());
  init(vss,iss,false);
  if(!bvh.empty()) {
    int id=(int)bvh.size()-1;
    std::vector<Node<int,BBoxExact>> meshBVH=mesh.getBVH();
    for(auto& n:meshBVH) {
      n._l=n._l>=0? n._l+(int)bvh.size(): -1;
      n._r=n._r>=0? n._r+(int)bvh.size(): -1;
      n._parent=n._parent>=0? n._parent+(int)bvh.size(): -1;
    }
    bvh.insert(bvh.end(),meshBVH.begin(),meshBVH.end());
    int id1=(int)bvh.size()-1;
    std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
    Eigen::Matrix<int,2,1> pair(id,id1);
    edgeMap[pair]=std::make_pair(pair[0],pair[1]);
    Node<int,BBoxExact>::buildBVHBottomUp(bvh,edgeMap);
  } else bvh.insert(bvh.end(),mesh.getBVH().begin(),mesh.getBVH().end());
  _bvh=bvh;
  while(Node<int,BBoxExact>::optimizeBVH(_bvh))
    Node<int,BBoxExact>::parityCheck(_bvh);
  Node<int,BBoxExact>::compact(_bvh);
  return *this;
}
MeshExact::Vec2T MeshExact::project(const Vec3T& d) const {
  Vec2T rng;
  rng[0]=std::numeric_limits<double>::max();
  rng[1]=-rng[0];
  for(const auto& v:_vss) {
    rng[0]=std::min<T>(rng[0],v.dot(d));
    rng[1]=std::max<T>(rng[1],v.dot(d));
  }
  return rng;
}
MeshExact::Vec3T MeshExact::support(const Vec3T& D,int& id) const {
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
void MeshExact::writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  os.setRelativeIndex();
  for(auto v:_vss)
    vss.push_back((ROT(trans)*v+CTR(trans)).template cast<double>());
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(_iss.begin(),_iss.end(),VTKWriter<double>::TRIANGLE,true);
}
void MeshExact::initBss() {
  //bss structure: char 00-edge[321]-vertices[321]
  _bss.resize(_iss.size());
  std::unordered_set<int> verticesStat;
  std::unordered_set<Eigen::Matrix<int,2,1>,EdgeHash> edgeStat;
  for(int i=0; i<(int)_iss.size(); ++i) {
    int cnt=0;
    _bss[i]=0;
    //vertices
    for(int j=0; j<3; ++j) {
      if(verticesStat.find(_iss[i][j])==verticesStat.end()) {
        verticesStat.insert(_iss[i][j]);
        _bss[i]|=1<<cnt;
      }
      cnt++;
    }
    for(int j=0; j<3; ++j) {
      Eigen::Matrix<int,2,1> edge(_iss[i][j],_iss[i][(j+1)%3]);
      sort2(edge[0],edge[1]);
      if(edgeStat.find(edge)==edgeStat.end()) {
        edgeStat.insert(edge);
        _bss[i]|=1<<cnt;
      }
      cnt++;
    }
  }
}
}
