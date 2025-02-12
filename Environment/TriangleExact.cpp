#include "TriangleExact.h"
#include <Utils/Pragma.h>
#include <Utils/IO.h>

namespace PHYSICSMOTION {
TriangleExact::TriangleExact():_vid(-1,-1,-1),_eNId(-1,-1,-1) {}
TriangleExact::TriangleExact(const Vec3T& a,const Vec3T& b,const Vec3T& c):_vid(-1,-1,-1),_eNId(-1,-1,-1) {
  _v[0]=a;
  _v[1]=b;
  _v[2]=c;
  reset();
}
bool TriangleExact::read(std::istream& is,IOData*) {
  readBinaryData(_vid,is);
  readBinaryData(_eNId,is);
  for(int d=0; d<3; d++) {
    readBinaryData(_v[d],is);
    readBinaryData(_n[d],is);
    readBinaryData(_nO[d],is);
  }
  readBinaryData(_nOSqr,is);
  readBinaryData(_nt,is);
  readBinaryData(_nTnN,is);
  readBinaryData(_invM,is);
  return is.good();
}
bool TriangleExact::write(std::ostream& os,IOData*) const {
  writeBinaryData(_vid,os);
  writeBinaryData(_eNId,os);
  for(int d=0; d<3; d++) {
    writeBinaryData(_v[d],os);
    writeBinaryData(_n[d],os);
    writeBinaryData(_nO[d],os);
  }
  writeBinaryData(_nOSqr,os);
  writeBinaryData(_nt,os);
  writeBinaryData(_nTnN,os);
  writeBinaryData(_invM,os);
  return os.good();
}
std::shared_ptr<SerializableBase> TriangleExact::copy() const {
  return std::shared_ptr<SerializableBase>(new TriangleExact);
}
std::string TriangleExact::type() const {
  return typeid(TriangleExact).name();
}
const TriangleExact::Vec3T& TriangleExact::v(int i) const {
  return _v[i];
}
const TriangleExact::Vec3T& TriangleExact::normal() const {
  return _nt;
}
BBoxExact TriangleExact::getBB() const {
  return BBoxExact(_v[0],_v[1],_v[2]);
}
TriangleExact::Vec4T TriangleExact::getSBB() const {
  Vec3T center=(_v[0]+_v[1]+_v[2])/3.0;
  T rad=std::max(std::max((center-_v[0]).norm(),(center-_v[1]).norm()),(center-_v[2]).norm());
  return Vec4T(center[0],center[1],center[2],rad);
}
TriangleExact::Vec3T TriangleExact::bary(const Vec3T& p) const {
  Mat3X2T LHS;
  LHS.col(0)=_v[1]-_v[0];
  LHS.col(1)=_v[2]-_v[0];
  Mat2T LHS2=LHS.transpose()*LHS;
  Vec2T ab=LHS2.inverse()*(LHS.transpose()*(p-_v[0]));
  return Vec3T(1-ab.sum(),ab[0],ab[1]);
}
TriangleExact::Vec3T TriangleExact::interp(const Vec3T& bary) const {
  return _v[0]*bary[0]+_v[1]*bary[1]+_v[2]*bary[2];
}
TriangleExact::Vec2T TriangleExact::project(const Vec3T& d) const {
  T a=_v[0].dot(d);
  T b=_v[1].dot(d);
  T c=_v[2].dot(d);
  return Vec2T(std::min(a,std::min(b,c)),std::max(a,std::max(b,c)));
}
bool TriangleExact::intersect(const BBoxExact& bb) const {
  std::function<bool(const Vec2T&,const Vec2T&)> testSegment=[&](const Vec2T& a,const Vec2T& b) {
    return a.x() > b.y() || b.x() > a.y();
  };
  std::function<bool(const Vec3T&,const TriangleExact&,const BBoxExact&)> sep=[&](const Vec3T& d,const TriangleExact& A,const BBoxExact& B) {
    return testSegment(A.project(d),B.project(d));
  };
  Vec3T d=(_v[1]-_v[0]).cross(_v[2]-_v[0]);
  if(testSegment(bb.project(d),Vec2T::Constant(_v[0].dot(d))))
    return false;
  for(int i=0; i<3; i++)
    if(sep(Vec3T::Unit(i),*this,bb))
      return false;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      if(sep(Vec3T::Unit(i).cross(_v[j]-_v[(j+1)%3]),*this,bb))
        return false;
  return true;
}
TriangleExact::T TriangleExact::dirGrad(const Eigen::Matrix<int,2,1>& evid,const Vec3T& D) const {
  for(int d=0; d<3; d++)
    if(_vid[d]!=evid[0] && _vid[d]!=evid[1]) {
      ASSERT((_vid[(d+1)%3]==evid[0] && _vid[(d+2)%3]==evid[1]) || (_vid[(d+2)%3]==evid[0] && _vid[(d+1)%3]==evid[1]))
      return -_n[d].dot(D);
    }
  ASSERT_MSGV(false,"Cannot find vertex opposite edge (%d,%d)",evid[0],evid[1])
  return 0;
}
std::pair<int,TriangleExact::Vec3T> TriangleExact::moveFromEdge(const Eigen::Matrix<int,2,1>& evid,const Vec3T& p0,Vec3T D) const {
  T t[2];
  int vid=-1;
  for(int d=0; d<3; d++)
    if(_vid[d]!=evid[0] && _vid[d]!=evid[1]) {
      ASSERT((_vid[(d+1)%3]==evid[0] && _vid[(d+2)%3]==evid[1]) || (_vid[(d+2)%3]==evid[0] && _vid[(d+1)%3]==evid[1]))
      vid=d;
      break;
    }
  D-=_nTnN*D;
  for(int off=0; off<2; off++) {
    const Vec3T& tp0=_v[(vid+off+2)%3];
    const Vec3T& td=_nO[(vid+off+1)%3];
    Vec2T RHS2;
    Mat2T LHS2;
    Mat3X2T LHS;
    LHS.col(0)=D;
    LHS.col(1)=-td;
    LHS2=LHS.transpose()*LHS;
    RHS2=LHS.transpose()*(tp0-p0);
    //invert
    T detLHS2=LHS2.determinant();
    if(detLHS2==0)
      t[off]=2;
    else {
      std::swap(LHS2(0,0),LHS2(1,1));
      LHS2(0,1)=-LHS2(0,1);
      LHS2(1,0)=-LHS2(1,0);
      LHS2/=detLHS2;
      t[off]=(LHS2*RHS2)[0];
    }
  }
  ASSERT(t[0]>0 || t[1]>0)
  if(t[0]<0) {
    if(t[1]<1)
      return std::make_pair((vid+2)%3,p0+D*t[1]);
    else return std::make_pair(-1,p0+D);
  } else if(t[1]<0) {
    if(t[0]<1)
      return std::make_pair((vid+1)%3,p0+D*t[0]);
    else return std::make_pair(-1,p0+D);
  } else if(t[0]>1 && t[1]>1)
    return std::make_pair(-1,p0+D);
  else if(t[0]<t[1])
    return std::make_pair((vid+1)%3,p0+D*t[0]);
  else if(t[1]<t[0])
    return std::make_pair((vid+2)%3,p0+D*t[1]);
  else
    return std::make_pair(vid,p0+D*t[0]);
}
void TriangleExact::calcPointDist(const Vec3T& pt,T& sqrDistance,Vec3T& cp,Vec3T& b,Eigen::Matrix<int,2,1>& feat) const {
  feat=Eigen::Matrix<int,2,1>(-1,-1);
  bool outside=false;
  for(int d=0; d<3; d++) {
    if((pt-_v[(d+1)%3]).dot(_n[d])>=0) {
      outside=true;
      T t=(pt-_v[(d+1)%3]).dot(_nO[d])/_nOSqr[d];
      if(t>0&&t<1) {
        //Edge's dirichlet region
        b.setZero();
        b[(d+1)%3]=1-t;
        b[(d+2)%3]=t;
        feat=Eigen::Matrix<int,2,1>((d+1)%3,(d+2)%3);
        cp=b[0]*_v[0]+b[1]*_v[1]+b[2]*_v[2];
        sqrDistance=(pt-cp).squaredNorm();
        break;
      }
    }
  }
  if(!outside) {
    //Triangle's dirichlet region
    Vec2T RHS;
    RHS[0]=(pt-_v[0]).dot(_nO[1]);
    RHS[1]=(pt-_v[0]).dot(_nO[2]);
    b.segment<2>(1)=_invM*RHS;
    b[0]=1-b[1]-b[2];
    feat=Eigen::Matrix<int,2,1>(-1,-1);
    cp=b[0]*_v[0]+b[1]*_v[1]+b[2]*_v[2];
    sqrDistance=(pt-cp).squaredNorm();
  } else if(feat[0]==-1) {
    //Vertex's dirichlet region
    for(int d=0; d<3; d++) {
      T sqrDistanceNew=(pt-_v[d]).squaredNorm();
      if(d==0||sqrDistanceNew<sqrDistance) {
        b=Vec3T::Unit(d);
        feat=Eigen::Matrix<int,2,1>(d,-1);
        cp=_v[d];
        sqrDistance=sqrDistanceNew;
      }
    }
  }
}
void TriangleExact::reset() {
  for(int d=0; d<3; d++) {
    _nO[d]=_v[(d+2)%3]-_v[(d+1)%3];
    _nOSqr[d]=_nO[d].dot(_nO[d]);
    T t=(_v[d]-_v[(d+1)%3]).dot(_nO[d])/_nOSqr[d];
    _n[d]=_v[(d+1)%3]+_nO[d]*t-_v[d];
  }
  _nt=(_v[1]-_v[0]).cross(_v[2]-_v[0]);
  _nTnN=(_nt*_nt.transpose())/_nt.squaredNorm();
  Mat2T M;
  M(0,0)= _nO[1].dot(_nO[2]);
  M(1,0)= _nO[2].dot(_nO[2]);
  M(0,1)=-_nO[1].dot(_nO[1]);
  M(1,1)=-_nO[2].dot(_nO[1]);
  _invM=M.inverse();
}
}
