#include "BBoxExact.h"
#include <Utils/IO.h>
#include <Utils/VTKWriter.h>

namespace PHYSICSMOTION {
BBoxExact::BBoxExact()
  :_minC(Eigen::Matrix<double,3,1>::Constant(std::numeric_limits<double>::max()).template cast<T>()),
   _maxC(Eigen::Matrix<double,3,1>::Constant(-std::numeric_limits<double>::max()).template cast<T>()) {}
BBoxExact::BBoxExact(T halfX,T halfY,T halfZ):_minC(-halfX,-halfY,-halfZ),_maxC(halfX,halfY,halfZ) {}
BBoxExact::BBoxExact(const Vec3T& a,const Vec3T& b):_minC(a.cwiseMin(b)),_maxC(a.cwiseMax(b)) {}
BBoxExact::BBoxExact(const Vec3T& a,const Vec3T& b,const Vec3T& c):_minC(a.cwiseMin(b).cwiseMin(c)),_maxC(a.cwiseMax(b).cwiseMax(c)) {}
bool BBoxExact::read(std::istream& is,IOData*) {
  readBinaryData(_minC,is);
  readBinaryData(_maxC,is);
  return is.good();
}
bool BBoxExact::write(std::ostream& os,IOData*) const {
  writeBinaryData(_minC,os);
  writeBinaryData(_maxC,os);
  return os.good();
}
std::shared_ptr<SerializableBase> BBoxExact::copy() const {
  return std::shared_ptr<SerializableBase>(new BBoxExact(*this));
}
std::string BBoxExact::type() const {
  return typeid(BBoxExact).name();
}
const BBoxExact& BBoxExact::getBB() const {
  return *this;
}
bool BBoxExact::empty() const {
  return (_minC.array()>_maxC.array()).any();
}
void BBoxExact::getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                        std::vector<Eigen::Matrix<int,3,1>>& iss) const {
  makeBox(vss,iss,1,(_maxC-_minC).template cast<double>()/2);
  for(int i=0; i<(int)vss.size(); i++)
    vss[i]+=(_maxC+_minC).template cast<double>()/2;
}
bool BBoxExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                             T&,Eigen::Matrix<int,2,1>& feat,bool,
                             std::vector<Vec3T>*) const {
  if(BBoxExact::contain(pt)) {
    T minDist=(_maxC-_minC).maxCoeff(),dist;
    for(int d=0; d<3; d++) {
      dist=pt[d]-_minC[d];
      if(dist<minDist) {
        minDist=dist;
        n=-Vec3T::Unit(d)*dist;
        normal=-Vec3T::Unit(d);
        feat=Eigen::Matrix<int,2,1>(-1,d*2+0);
      }
      dist=_maxC[d]-pt[d];
      if(dist<minDist) {
        minDist=dist;
        n=Vec3T::Unit(d)*dist;
        normal=Vec3T::Unit(d);
        feat=Eigen::Matrix<int,2,1>(-1,d*2+1);
      }
    }
    hessian.setZero();
    return true;
  } else {
    int nrBlocked=0;
    Vec3T dist(pt);
    normal.setZero();
    for(int i=0; i<3; i++) {
      if (pt[i]<_minC[i]) {
        dist[i]=_minC[i];
        normal[i]=-1;
        nrBlocked++;
      } else if(pt[i]>_maxC[i]) {
        dist[i]=_maxC[i];
        normal[i]=1;
        nrBlocked++;
      }
    }
    n=dist-pt;
    //calculate feat
    if(nrBlocked==3)
      feat=Eigen::Matrix<int,2,1>(normalToVid(normal),-1);
    else if(nrBlocked==2)
      feat=Eigen::Matrix<int,2,1>(normalToEid0(normal),normalToEid1(normal));
    else if(nrBlocked==1)
      feat=Eigen::Matrix<int,2,1>(-1,normalToFid(normal));
    else {
      ASSERT_MSG(false,"Impossible configuration between point and BBoxExact!")
    }
    //adjust hessian
    if(n.isZero()) {
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
      Vec3T e=vertex(feat[0])-vertex(feat[1]);
      T eDotE=e.dot(e);
      Vec3T d=n-n.dot(e)*e/eDotE;
      T dDotD=d.dot(d);
      Vec3T nd=e.dot(d)*e/eDotE/dDotD;
      hessian=e*e.transpose()/eDotE-Mat3T::Identity();
      hessian+=d*d.transpose()/dDotD;
      hessian-=d*nd.transpose();
    }
    return false;
  }
}
void BBoxExact::scale(T coef) {
  _minC*=coef;
  _maxC*=coef;
}
void BBoxExact::setUnion(const BBoxExact& other) {
  for(int d=0; d<3; d++) {
    _minC[d]=std::min<T>(_minC[d],other._minC[d]);
    _maxC[d]=std::max<T>(_maxC[d],other._maxC[d]);
  }
}
void BBoxExact::setUnion(const Vec3T& other) {
  for(int d=0; d<3; d++) {
    _minC[d]=std::min<T>(_minC[d],other[d]);
    _maxC[d]=std::max<T>(_maxC[d],other[d]);
  }
}
void BBoxExact::extendUnion(const T x0){
  for(int d=0;d<3;d++){
    _minC[d]-=x0/2.0;
    _maxC[d]+=x0/2.0;
  }
}
const BBoxExact::Vec3T& BBoxExact::minCorner() const {
  return _minC;
}
const BBoxExact::Vec3T& BBoxExact::maxCorner() const {
  return _maxC;
}
BBoxExact::Vec3T& BBoxExact::minCorner() {
  return _minC;
}
BBoxExact::Vec3T& BBoxExact::maxCorner() {
  return _maxC;
}
bool BBoxExact::intersect(const BBoxExact& other) const {
  for(int i=0; i<3; i++)
    if(_maxC[i]<other._minC[i] || other._maxC[i]<_minC[i])
      return false;
  return true;
}
bool BBoxExact::contain(const Vec3T& pt) const {
  return (_minC.array()<=pt.array()).all() && (_maxC.array()>=pt.array()).all();
}
BBoxExact BBoxExact::enlargedEps(T eps) const {
  return enlarged((maxCorner()-minCorner())*eps);
}
BBoxExact BBoxExact::enlarged(const Vec3T& ext) const {
  return BBoxExact(_minC-ext,_maxC+ext);
}
BBoxExact::T BBoxExact::distToSqr(const BBoxExact& other) const {
  Vec3T dist=Vec3T::Zero();
  for(int i=0; i<3; i++) {
    if (other._maxC[i]<_minC[i])
      dist[i]=other._maxC[i]-_minC[i];
    else if(other._minC[i]>_maxC[i])
      dist[i]=other._minC[i]-_maxC[i];
  }
  return dist.squaredNorm();
}
BBoxExact::T BBoxExact::distToSqr(const Vec3T& pt) const {
  Vec3T dist=Vec3T::Zero();
  for(int i=0; i<3; i++) {
    if (pt[i]<_minC[i])
      dist[i]=pt[i]-_minC[i];
    else if(pt[i]>_maxC[i])
      dist[i]=pt[i]-_maxC[i];
  }
  return dist.squaredNorm();
}
bool BBoxExact::operator==(BBoxExact &bbox) {
  return (_minC == bbox._minC) && (_maxC == bbox._maxC);
}
//for GJK
BBoxExact::Vec3T BBoxExact::support(const Vec3T& D,int& id) const {
  id=-1;
  Vec3T v,maxV;
  T dist,maxDist=-std::numeric_limits<double>::max();
  for(int i=0; i<8; i++) {
    v[0]=(i&1)?_maxC[0]:_minC[0];
    v[1]=(i&2)?_maxC[1]:_minC[1];
    v[2]=(i&4)?_maxC[2]:_minC[2];
    dist=v.dot(D);
    if(dist>maxDist) {
      maxDist=dist;
      maxV=v;
      id=i;
    }
  }
  return maxV;
}
//for SAT
BBoxExact::Vec2T BBoxExact::project(const Vec3T& d) const {
  Vec3T ctr=(_maxC+_minC)/2,ext=(_maxC-_minC)/2;
  T ctrD=d.dot(ctr),delta=(d.array()*ext.array()).abs().sum();
  return Vec2T(ctrD-delta,ctrD+delta);
}
std::vector<ShapeExact::Facet> BBoxExact::facets() const {
  Vec3T pos;
  std::vector<Facet> facets;
  for(int a=0; a<3; a++) {
    int a2=(a+1)%3;
    int a3=(a+2)%3;
    for(int d=0; d<2; d++) {
      Facet f;
      //np
      if(d==0) {
        f._n=Vec3T::Unit(a);
        pos[a]=_maxC[a];
      } else {
        f._n=-Vec3T::Unit(a);
        pos[a]=_minC[a];
        std::swap(a2,a3);
      }
      //v0
      pos[a2]=_minC[a2];
      pos[a3]=_minC[a3];
      f._boundary.push_back(pos);
      //v1
      pos[a2]=_maxC[a2];
      pos[a3]=_minC[a3];
      f._boundary.push_back(pos);
      //v2
      pos[a2]=_maxC[a2];
      pos[a3]=_maxC[a3];
      f._boundary.push_back(pos);
      //v3
      pos[a2]=_minC[a2];
      pos[a3]=_maxC[a3];
      f._boundary.push_back(pos);
      //insert
      facets.push_back(f);
    }
  }
  return facets;
}
std::vector<ShapeExact::Edge> BBoxExact::edges() const {
  Edge e;
  std::vector<Edge> edges;
  for(int a=0; a<3; a++) {
    int a2=(a+1)%3;
    int a3=(a+2)%3;
    //e00
    e._a[a]=_minC[a];
    e._a[a2]=_minC[a2];
    e._a[a3]=_minC[a3];
    e._b=e._a;
    e._b[a]=_maxC[a];
    edges.push_back(e);
    //e10
    e._a[a]=_minC[a];
    e._a[a2]=_maxC[a2];
    e._a[a3]=_minC[a3];
    e._b=e._a;
    e._b[a]=_maxC[a];
    edges.push_back(e);
    //e11
    e._a[a]=_minC[a];
    e._a[a2]=_maxC[a2];
    e._a[a3]=_maxC[a3];
    e._b=e._a;
    e._b[a]=_maxC[a];
    edges.push_back(e);
    //e01
    e._a[a]=_minC[a];
    e._a[a2]=_minC[a2];
    e._a[a3]=_maxC[a3];
    e._b=e._a;
    e._b[a]=_maxC[a];
    edges.push_back(e);
  }
  return edges;
}
void BBoxExact::writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const {
  Eigen::Matrix<double,3,1> halfExt=(_maxC-_minC).template cast<double>()/2;
  Eigen::Matrix<double,3,1> ctr=(_minC+_maxC).template cast<double>()/2;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  makeBox(vss,iss,2,halfExt);
  os.setRelativeIndex();
  for(auto& v:vss)
    v=(ROT(trans)*(v+ctr).template cast<T>()+CTR(trans)).template cast<double>();
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::TRIANGLE,true);
}
//helper
void BBoxExact::makeGrid(std::vector<Eigen::Matrix<double,3,1>>& vss,std::vector<Eigen::Matrix<int,3,1>>& iss,
                         int RESX,int RESY,const Eigen::Matrix<double,3,1>& ctr,const Eigen::Matrix<double,3,1>& d0,const Eigen::Matrix<double,3,1>& d1) {
#define ID(X,Y) (X)*(RESY+1)+(Y)+off
  //skip degenerate case
  if(d0.isZero() || d1.isZero())
    return;
  int off=(int)vss.size();
  for(int x=0; x<=RESX; x++)
    for(int y=0; y<=RESY; y++) {
      vss.push_back(ctr+d0*(2*x/(double)RESX-1)+d1*(2*y/(double)RESY-1));
      if(x<RESX && y<RESY) {
        iss.push_back(Eigen::Matrix<int,3,1>(ID(x,y),ID(x+1,y),ID(x+1,y+1)));
        iss.push_back(Eigen::Matrix<int,3,1>(ID(x,y),ID(x+1,y+1),ID(x,y+1)));
      }
    }
#undef ID
}
void BBoxExact::makeBox(std::vector<Eigen::Matrix<double,3,1>>& vss,std::vector<Eigen::Matrix<int,3,1>>& iss,
                        int RES,const Eigen::Matrix<double,3,1>& halfExt) {
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>(-halfExt[0],0,0),Eigen::Matrix<double,3,1>(0,0,halfExt[2]),Eigen::Matrix<double,3,1>(0,halfExt[1],0));
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>( halfExt[0],0,0),Eigen::Matrix<double,3,1>(0,halfExt[1],0),Eigen::Matrix<double,3,1>(0,0,halfExt[2]));
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>(0,-halfExt[1],0),Eigen::Matrix<double,3,1>(halfExt[0],0,0),Eigen::Matrix<double,3,1>(0,0,halfExt[2]));
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>(0, halfExt[1],0),Eigen::Matrix<double,3,1>(0,0,halfExt[2]),Eigen::Matrix<double,3,1>(halfExt[0],0,0));
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>(0,0,-halfExt[2]),Eigen::Matrix<double,3,1>(0,halfExt[1],0),Eigen::Matrix<double,3,1>(halfExt[0],0,0));
  makeGrid(vss,iss,RES,RES,Eigen::Matrix<double,3,1>(0,0, halfExt[2]),Eigen::Matrix<double,3,1>(halfExt[0],0,0),Eigen::Matrix<double,3,1>(0,halfExt[1],0));
}
int BBoxExact::normalToVid(const Vec3T& normal) {
  return ((normal.cast<int>()+Eigen::Matrix<int,3,1>(1,1,1))/2).dot(Eigen::Matrix<int,3,1>(1,2,4));
}
int BBoxExact::normalToEid0(const Vec3T& normal) {
  for(int i=0; i<3; i++)
    if(normal[i]==0)
      return normalToVid(normal-Vec3T::Unit(i));
  return normalToVid(normal);
}
int BBoxExact::normalToEid1(const Vec3T& normal) {
  for(int i=0; i<3; i++)
    if(normal[i]==0)
      return normalToVid(normal+Vec3T::Unit(i));
  return normalToVid(normal);
}
int BBoxExact::normalToFid(const Vec3T& normal) {
  for(int i=0; i<3; i++)
    if(normal[i]!=0)
      return normal[i]==-1?(i*2+0):(i*2+1);
  return -1;
}
BBoxExact::Vec3T BBoxExact::vertex(int id) const {
  return Vec3T((id&1)?_maxC[0]:_minC[0],
               (id&2)?_maxC[1]:_minC[1],
               (id&4)?_maxC[2]:_minC[2]);
}
}
