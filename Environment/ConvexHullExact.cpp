#include "ConvexHullExact.h"
#include "EnvironmentUtils.h"
#include <Utils/DebugGradient.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include <Utils/IO.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <iomanip>
#include <stack>

namespace PHYSICSMOTION {
ConvexHullExact::ConvexHullExact() {}
ConvexHullExact::ConvexHullExact(const std::string& path) {
  Assimp::Importer importer;
  const aiScene *scene=importer.ReadFile(path.c_str(),aiProcess_JoinIdenticalVertices);
  ASSERT_MSGV(scene,"Mesh %s is empty!",path.c_str())
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  init(scene,NULL,vss);
}
ConvexHullExact::ConvexHullExact(const aiScene* scene) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  init(scene,NULL,vss);
}
ConvexHullExact::ConvexHullExact(const MeshExact& m) {
  init<T>(m.vss());
}
ConvexHullExact::ConvexHullExact(const std::vector<Eigen::Matrix<double,3,1>>& vss) {
  init<double>(vss);
}
void ConvexHullExact::init(const aiScene* scene,const aiNode* node,
                           std::vector<Eigen::Matrix<double,3,1>>& vss) {
  if(!node) {
    init(scene,scene->mRootNode,vss);
    init<double>(vss);
  } else {
    //add children
    for(int i=0; i<(int)node->mNumChildren; i++)
      init(scene,node->mChildren[i],vss);
    //add mesh
    for(int idm=0; idm<(int)node->mNumMeshes; idm++) {
      const aiMesh* m=scene->mMeshes[node->mMeshes[idm]];
      for(int i=0; i<(int)m->mNumVertices; i++) {
        const aiVector3D& v=m->mVertices[i];
        vss.push_back(Eigen::Matrix<double,3,1>(v.x,v.y,v.z));
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
void ConvexHullExact::init(const std::vector<Eigen::Matrix<T2,3,1>>& pss) {
  if(pss.empty())
    return;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  vss.resize(pss.size());
  for(int i=0; i<(int)pss.size(); i++)
    vss[i]=pss[i].template cast<double>();
  std::vector<Eigen::Matrix<int,3,1>> iss;
  makeConvex(vss,iss);
  init<double>(vss,iss);
}
template <typename T2>
void ConvexHullExact::init(const std::vector<Eigen::Matrix<T2,3,1>>& vss,
                           const std::vector<Eigen::Matrix<int,3,1>>& iss) {
  //vss
  if(vss.empty() || iss.empty())
    return;
  _vss.resize(vss.size());
  for(int i=0; i<(int)vss.size(); i++)
    _vss[i]=vss[i].template cast<T>();
  //iss
  _iss=iss;
  //tss
  _tss.resize(_iss.size());
  for(int i=0; i<(int)_tss.size(); i++) {
    _tss[i]=TriangleExact(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]);
    _tss[i]._vid=_iss[i];
  }
  //bvh
  if(!_vss.empty()) {
    _bvh.assign(_iss.size(),Node<int,BBoxExact>());
    for(int i=0; i<(int)_bvh.size(); i++) {
      Node<int,BBoxExact>& n=_bvh[i];
      //n._bb=BBoxExact(_tss[i].getBB());
      Vec4T SBB=_tss[i].getSBB();
      n._bb=BBoxExact(SBB.segment<3>(0),SBB[3]);
      n._nrCell=1;
      n._cell=i;
    }
    Node<int,BBoxExact>::buildBVHTriangleBottomUp(_bvh,iss,true);
  }
  //edge
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> ess;
  buildEdge(iss,ess);
  //ess
  int eid=0;
  _ess.resize(ess.size());
  _eNss.resize(_vss.size());
  for(const std::pair<const Eigen::Matrix<int,2,1>,std::pair<int,int>>& E:ess) {
    _ess[eid]=EdgeExact(_vss[E.first[0]],_vss[E.first[1]]);
    _ess[eid]._tNId=Eigen::Matrix<int,2,1>(E.second.first,E.second.second);
    _ess[eid]._vid=E.first;
    for(int d=0; d<2; d++) {
      int tid=_ess[eid]._tNId[d];
      TriangleExact& T=_tss[tid];
      bool found=false;
      for(int d2=0; d2<3; d2++)
        if(T._vid[d2]!=E.first[0] && T._vid[d2]!=E.first[1]) {
          T._eNId[d2]=eid;
          found=true;
          break;
        }
      ASSERT_MSGV(found,"Cannot find Edge (%d,%d) in Triangle %d",E.first[0],E.first[1],tid)
    }
    _eNss[E.first[0]].push_back(eid);
    _eNss[E.first[1]].push_back(eid);
    eid++;
  }
  parityCheck();
  //bss
  initBss();
}
bool ConvexHullExact::read(std::istream& is,IOData* dat) {
  MeshExact::read(is,dat);
  readBinaryData(_ess,is,dat);
  readBinaryData(_eNss,is,dat);
  return is.good();
}
bool ConvexHullExact::write(std::ostream& os,IOData* dat) const {
  MeshExact::write(os,dat);
  writeBinaryData(_ess,os,dat);
  writeBinaryData(_eNss,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ConvexHullExact::copy() const {
  std::shared_ptr<ConvexHullExact> ret(new ConvexHullExact);
  *ret=*this;
  return ret;
}
std::string ConvexHullExact::type() const {
  return typeid(ConvexHullExact).name();
}
void ConvexHullExact::parityCheck() const {
  for(int tid=0; tid<(int)_tss.size(); tid++)
    for(int d=0; d<3; d++) {
      int eid=_tss[tid]._eNId[d];
      ASSERT_MSGV(eid>=0,"Triangle %d's Edge %d is -1",tid,d)
      const EdgeExact& e=_ess[eid];
      ASSERT_MSGV(e._tNId[0]==tid || e._tNId[1]==tid,"Edge %d does not contain Traingle %d",eid,tid)
    }
  for(int eid=0; eid<(int)_ess.size(); eid++)
    for(int d=0; d<2; d++) {
      int tid=_ess[eid]._tNId[d];
      ASSERT_MSGV(tid>=0,"Edge %d's Triangle %d is -1",eid,d)
      const TriangleExact& t=_tss[tid];
      ASSERT_MSGV(t._eNId[0]==eid || t._eNId[1]==eid || t._eNId[2]==eid,"Triangle %d does not contain Edge %d",tid,eid)
      int vid=_ess[eid]._vid[d];
      ASSERT_MSGV(std::find(_eNss[vid].begin(),_eNss[vid].end(),eid)!=_eNss[vid].end(),"Cannot find Edge %d as neighbor of Vertex %d",eid,vid)
    }
  for(int vid=0; vid<(int)_eNss.size(); vid++)
    for(int eid:_eNss[vid]) {
      const EdgeExact& e=_ess[eid];
      ASSERT_MSGV(e._vid[0]==vid || e._vid[1]==vid,"Edge %d does not contain Vertex %d",eid,vid)
    }
}
bool ConvexHullExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                                   T&,Eigen::Matrix<int,2,1>& feat,bool cache,
                                   std::vector<Vec3T>* history) const {
  //initialize
  int bestVid=-1;
  T distSqr,bestDistSqr;
  if(cache) {
    if(feat[0]>=0)
      bestVid=feat[0];
    else bestVid=_tss[feat[1]]._vid[0];
  } else {
    for(int vid=0; vid<(int)_vss.size(); vid++) {
      distSqr=(_vss[vid]-pt).squaredNorm();
      if(bestVid==-1 || distSqr<bestDistSqr) {
        bestDistSqr=distSqr;
        bestVid=vid;
      }
    }
  }
  //local optimization
  Vec3T p=_vss[bestVid];
  if(history)
    history->push_back(p);
  int currentVid=bestVid;
  int currentEid=-1;
  int currentTid=-1;
  while(true)
    if(currentTid>=0) {
      ASSERT(currentVid==-1 && currentEid>=0)
      const EdgeExact& E=_ess[currentEid];
      const TriangleExact& T=_tss[currentTid];
      std::pair<int,Vec3T> pT=T.moveFromEdge(E._vid,p,pt-p);
      p=pT.second;
      if(history)
        history->push_back(p);
      if(pT.first==-1) {
        normal=T.normal();
        feat=Eigen::Matrix<int,2,1>(-1,currentTid);
        break;
      } else if(T._vid[pT.first]!=E._vid[0] && T._vid[pT.first]!=E._vid[1]) {
        currentVid=T._vid[pT.first];
        currentEid=-1;
        currentTid=-1;
      } else {
        currentVid=-1;
        currentEid=T._eNId[pT.first];
        currentTid=-1;
      }
    } else if(currentVid>=0) {
      ASSERT(currentEid==-1 && currentTid==-1)
      //check edge
      T bestDirGrad=0,dirGrad;
      for(int eNId:_eNss[currentVid]) {
        dirGrad=_ess[eNId].dirGrad(currentVid,pt-p);
        if(dirGrad>bestDirGrad) {
          currentEid=eNId;
          bestDirGrad=dirGrad;
        }
      }
      if(bestDirGrad>0)
        currentVid=-1;
      //no edge, return
      if(currentVid>=0) {
        int eNId=_eNss[currentVid][0];
        normal=_tss[_ess[eNId]._tNId[0]].normal();
        feat=Eigen::Matrix<int,2,1>(currentVid,-1);
        break;
      }
    } else if(currentEid>=0) {
      ASSERT(currentVid==-1 && currentTid==-1)
      const EdgeExact& E=_ess[currentEid];
      std::pair<int,Vec3T> pE=E.moveFromVertex(p,pt-p);
      p=pE.second;
      if(history)
        history->push_back(p);
      if(pE.first>=0) {
        currentVid=pE.first;
        currentEid=-1;
      } else {
        //check triangle
        for(int d=0; d<2; d++) {
          const TriangleExact& T=_tss[E._tNId[d]];
          if(T.dirGrad(E._vid,pt-p)>0) {
            currentTid=E._tNId[d];
            break;
          }
        }
        //no triangle, return
        if(currentTid==-1) {
          normal=_tss[E._tNId[0]].normal();
          feat=E._vid;
          break;
        }
      }
    } else {
      ASSERT_MSG(false,"Invalid configuration with vid=eid=tid=-1")
    }
  n=p-pt;
  distSqr=n.squaredNorm();
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
void ConvexHullExact::scale(T coef) {
  std::vector<Vec3T> vss=_vss;
  for(Vec3T& v:vss)
    v*=coef;
  init(vss,_iss);
}
const std::vector<EdgeExact>& ConvexHullExact::ess() const {
  return _ess;
}
Eigen::Matrix<ConvexHullExact::T,4,1> ConvexHullExact::plane(int i) const {
  Eigen::Matrix<T,4,1> ret;
  Eigen::Matrix<T,3,1> a=_vss[_iss[i][0]];
  Eigen::Matrix<T,3,1> b=_vss[_iss[i][1]];
  Eigen::Matrix<T,3,1> c=_vss[_iss[i][2]];
  Eigen::Matrix<T,3,1> normal=(b-a).cross(c-a);
  ret.template segment<3>(0)=normal;
  ret[3]=-a.dot(normal);
  while(ret.squaredNorm()<1)
    ret*=2;
  return ret;
}
int ConvexHullExact::nrPlane() const {
  return (int)_iss.size();
}
//for GJK
ConvexHullExact::Vec3T ConvexHullExact::support(const Vec3T& D,int& id) const {
  id=0;
  T dist,maxDist=_vss[0].dot(D);
  while(true) {
    bool updated=false;
    for(int eid:_eNss[id]) {
      const EdgeExact& e=_ess[eid];
      int oid=e._vid[0]==id?e._vid[1]:e._vid[0];
      dist=_vss[oid].dot(D);
      if(dist>maxDist) {
        updated=true;
        maxDist=dist;
        id=oid;
        break;
      }
    }
    if(!updated)
      break;
  }
  return _vss[id];
}
//for SAT
std::vector<ConvexHullExact::Facet> ConvexHullExact::facets() const {
  //build edge map
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> emap;
  for(int i=0; i<(int)_iss.size(); i++)
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,2,1> e(_iss[i][d],_iss[i][(d+1)%3]);
      sort2(e[0],e[1]);
      if(emap.find(e)==emap.end())
        emap[e]=std::make_pair(i,-1);
      else emap[e].second=i;
    }
  //merge face
  std::vector<int> parent(_iss.size(),-1);
  for(const auto& e:emap) {
    const Vec3T& na=_tss[e.second.first].normal();
    const Vec3T& nb=_tss[e.second.second].normal();
    if(na.template cast<double>().normalized().dot(nb.template cast<double>().normalized())>=1-_epsEdge) {
      int pa=e.second.first;
      int pb=e.second.second;
      while(parent[pa]>=0)
        pa=parent[pa];
      while(parent[pb]>=0)
        pb=parent[pb];
      if(pa!=pb)
        parent[pa]=pb;
    }
  }
  std::unordered_map<int,std::unordered_set<int>> facetMap;
  for(int i=0; i<(int)_iss.size(); i++) {
    int p=i;
    while(parent[p]>=0)
      p=parent[p];
    for(int d=0; d<3; d++)
      facetMap[p].insert(_iss[i][d]);
  }
  //use project convex hull to compute facets
  std::vector<Facet> facets;
  for(const auto& f:facetMap) {
    Facet facet;
    std::vector<Eigen::Matrix<double,3,1>> vss;
    for(const auto& vid:f.second)
      vss.push_back(_vss[vid].template cast<double>());
    //safe normal computation
    T maxNorm=0;
    Eigen::Matrix<double,3,1> n,maxN;
    for(int i=0; i<(int)vss.size(); i++)
      for(int j=i+1; j<(int)vss.size(); j++)
        for(int k=j+1; k<(int)vss.size(); k++) {
          n=(vss[j]-vss[i]).cross(vss[k]-vss[i]);
          if(n.norm()>maxNorm) {
            maxNorm=n.norm();
            maxN=n;
          }
        }
    //normal direction
    facet._n=maxN.normalized().template cast<T>();
    Vec2T rng=project(facet._n);
    if((rng[0]+rng[1])/2>vss[0].template cast<T>().dot(facet._n))
      facet._n*=-1;
    //projected convex hull
    Vec3T::Index id;
    facet._n.cwiseAbs().minCoeff(&id);
    T d=vss[0].template cast<T>().dot(facet._n);
    Vec3T t1=Vec3T::Unit(id).cross(facet._n).template cast<double>().normalized().template cast<T>();
    Vec3T t2=facet._n.cross(t1);
    for(auto& v:vss)
      v=Eigen::Matrix<double,3,1>(t1.template cast<double>().dot(v),t2.template cast<double>().dot(v),0);
    makeConvexProject(vss);
    for(const auto& v:vss)
      facet._boundary.push_back(t1*v[0]+t2*v[1]+facet._n*d);
    facets.push_back(facet);
  }
  return facets;
}
std::vector<ConvexHullExact::Edge> ConvexHullExact::edges() const {
  std::vector<Edge> edges;
  //build edge map
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> emap;
  for(int i=0; i<(int)_iss.size(); i++)
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,2,1> e(_iss[i][d],_iss[i][(d+1)%3]);
      sort2(e[0],e[1]);
      if(emap.find(e)==emap.end())
        emap[e]=std::make_pair(i,-1);
      else emap[e].second=i;
    }
  //detect feature edge
  for(const auto& e:emap) {
    const Vec3T& na=_tss[e.second.first].normal();
    const Vec3T& nb=_tss[e.second.second].normal();
    if(na.template cast<double>().normalized().dot(nb.template cast<double>().normalized())<1-_epsEdge)
      edges.push_back(Edge({_vss[e.first[0]],_vss[e.first[1]]}));
  }
  return edges;
}
void ConvexHullExact::writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> ess;
  os.setRelativeIndex();
  for(auto v:_vss)
    vss.push_back((ROT(trans)*v+CTR(trans)).template cast<double>());
  for(auto e:_ess)
    ess.push_back(e._vid);
  os.setRelativeIndex();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(ess.begin(),ess.end(),VTKWriter<double>::LINE,true);
}
ConvexHullExact::T ConvexHullExact::_epsEdge=1e-3;
}
