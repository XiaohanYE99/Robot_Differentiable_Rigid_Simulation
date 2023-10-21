#include "EdgeExact.h"
#include <Utils/IO.h>

namespace PHYSICSMOTION {
EdgeExact::EdgeExact():_vid(-1,-1),_tNId(-1,-1) {}
EdgeExact::EdgeExact(const Vec3T& a,const Vec3T& b):_vid(-1,-1),_tNId(-1,-1) {
  _v[0]=a;
  _v[1]=b;
  reset();
}
bool EdgeExact::read(std::istream& is,IOData*) {
  readBinaryData(_vid,is);
  readBinaryData(_tNId,is);
  readBinaryData(_v[0],is);
  readBinaryData(_v[1],is);
  readBinaryData(_d,is);
  readBinaryData(_dTdN,is);
  readBinaryData(_lenSqr,is);
  return is.good();
}
bool EdgeExact::write(std::ostream& os,IOData*) const {
  writeBinaryData(_vid,os);
  writeBinaryData(_tNId,os);
  writeBinaryData(_v[0],os);
  writeBinaryData(_v[1],os);
  writeBinaryData(_d,os);
  writeBinaryData(_dTdN,os);
  writeBinaryData(_lenSqr,os);
  return os.good();
}
std::shared_ptr<SerializableBase> EdgeExact::copy() const {
  return std::shared_ptr<SerializableBase>(new EdgeExact);
}
std::string EdgeExact::type() const {
  return typeid(EdgeExact).name();
}
std::pair<int,EdgeExact::Vec3T> EdgeExact::moveFromVertex(const Vec3T& d0,const Vec3T& D) const {
  Vec3T dir=_dTdN*D;
  T sgn=dir.dot(_d);
  if(sgn>0) {
    T lenSqr=(_v[1]-d0).squaredNorm();
    if(lenSqr>D.dot(dir))
      return std::make_pair(-1,d0+dir);
    else return std::make_pair(_vid[1],_v[1]);
  } else if(sgn<0) {
    T lenSqr=(_v[0]-d0).squaredNorm();
    if(lenSqr>D.dot(dir))
      return std::make_pair(-1,d0+dir);
    else return std::make_pair(_vid[0],_v[0]);
  } else {
    return std::make_pair(-1,d0+dir);
  }
}
EdgeExact::T EdgeExact::dirGrad(int vid,const Vec3T& D) const {
  if(_vid[0]==vid)
    return D.dot(_d);
  else {
    ASSERT(_vid[1]==vid)
    return -D.dot(_d);
  }
}
void EdgeExact::reset() {
  _d=_v[1]-_v[0];
  _dTdN=_d*_d.transpose()/_d.dot(_d);
  _lenSqr=_d.dot(_d);
}
}
