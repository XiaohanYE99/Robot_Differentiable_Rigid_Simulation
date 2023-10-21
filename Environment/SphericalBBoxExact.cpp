#include "SphericalBBoxExact.h"
#include <Utils/IO.h>

namespace PHYSICSMOTION {
SphericalBBoxExact::SphericalBBoxExact() {}
SphericalBBoxExact::SphericalBBoxExact(T rad):BBoxExact(Vec3T(0,0,0),Vec3T(0,0,0)),_rad(rad),_radSqr(rad*rad) {
  _bb=enlarged(Vec3T(_rad,_rad,_rad));
}
SphericalBBoxExact::SphericalBBoxExact(T a,T rad):BBoxExact(Vec3T(-a,0,0),Vec3T(a,0,0)),_rad(rad),_radSqr(rad*rad) {
  _bb=enlarged(Vec3T(_rad,_rad,_rad));
}
SphericalBBoxExact::SphericalBBoxExact(T a,T b,T rad):BBoxExact(Vec3T(-a,-b,0),Vec3T(a,b,0)),_rad(rad),_radSqr(rad*rad) {
  _bb=enlarged(Vec3T(_rad,_rad,_rad));
}
SphericalBBoxExact::SphericalBBoxExact(T a,T b,T c,T rad):BBoxExact(Vec3T(-a,-b,-c),Vec3T(a,b,c)),_rad(rad),_radSqr(rad*rad) {
  _bb=enlarged(Vec3T(_rad,_rad,_rad));
}
bool SphericalBBoxExact::read(std::istream& is,IOData*) {
  BBoxExact::read(is,NULL);
  readBinaryData(_bb,is);
  readBinaryData(_rad,is);
  readBinaryData(_radSqr,is);
  return is.good();
}
bool SphericalBBoxExact::write(std::ostream& os,IOData*) const {
  BBoxExact::write(os,NULL);
  writeBinaryData(_bb,os);
  writeBinaryData(_rad,os);
  writeBinaryData(_radSqr,os);
  return os.good();
}
std::shared_ptr<SerializableBase> SphericalBBoxExact::copy() const {
  std::shared_ptr<SphericalBBoxExact> ret(new SphericalBBoxExact());
  *ret=*this;
  return ret;
}
std::string SphericalBBoxExact::type() const {
  return typeid(SphericalBBoxExact).name();
}
const BBoxExact& SphericalBBoxExact::getBB() const {
  return _bb;
}
bool SphericalBBoxExact::empty() const {
  return BBoxExact::empty() && _rad<=0;
}
void SphericalBBoxExact::getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                 std::vector<Eigen::Matrix<int,3,1>>& iss) const {
  makeBox(vss,iss,17,Eigen::Matrix<double,3,1>::Ones());
  Eigen::Matrix<double,3,1> halfExt=(_maxC-_minC).template cast<double>()/2;
  for(int i=0; i<(int)vss.size(); i++) {
    Eigen::Matrix<double,3,1>& v=vss[i];
    v=v.normalized()*(double)_rad;
    for(int d=0; d<3; d++)
      if(v[d]<0)
        v[d]-=halfExt[d];
      else v[d]+=halfExt[d];
    v+=(_maxC+_minC).template cast<double>()/2;
  }
}
bool SphericalBBoxExact::closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                                      T& rad,Eigen::Matrix<int,2,1>& feat,bool cache,
                                      std::vector<Vec3T>* history) const {
  rad=_rad;
  return BBoxExact::closestInner(pt,n,normal,hessian,
                                 rad,feat,cache,
                                 history);
}
void SphericalBBoxExact::scale(T coef) {
  BBoxExact::scale(coef);
  _rad*=coef;
}
bool SphericalBBoxExact::contain(const Vec3T& pt) const {
  return BBoxExact::distToSqr(pt)<=_radSqr;
}
bool SphericalBBoxExact::isRoundCube() const {
  return (_minC.array()==_maxC.array()).template cast<unsigned char>().sum()==0;
}
bool SphericalBBoxExact::isRoundBoard() const {
  return (_minC.array()==_maxC.array()).template cast<unsigned char>().sum()==1;
}
bool SphericalBBoxExact::isCapsule() const {
  return (_minC.array()==_maxC.array()).template cast<unsigned char>().sum()==2;
}
bool SphericalBBoxExact::isSphere() const {
  return (_minC.array()==_maxC.array()).template cast<unsigned char>().sum()==3;
}
SphericalBBoxExact::T SphericalBBoxExact::radius() const {
  return _rad;
}
//for SAT
std::vector<ShapeExact::Facet> SphericalBBoxExact::facets() const {
  ASSERT_MSG(isCapsule(),"We only support SAT for capsule!")
  return std::vector<Facet>();
}
std::vector<ShapeExact::Edge> SphericalBBoxExact::edges() const {
  ASSERT_MSG(isCapsule(),"We only support SAT for capsule!")
  return std::vector<Edge>({Edge({_minC,_maxC})});
}
}
