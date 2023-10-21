#ifndef CROSS_SPATIAL_UTILS_H
#define CROSS_SPATIAL_UTILS_H

#include "Pragma.h"
#include <Eigen/Dense>

namespace PHYSICSMOTION {
//cross
template <typename T>
Eigen::Matrix<T,3,3> cross(const Eigen::Matrix<T,3,1>& v) {
  Eigen::Matrix<T,3,3> ret;
  ret <<
      0,-v[2],v[1],
      v[2],0,-v[0],
      -v[1],v[0],0;
  return ret;
}
template <typename T>
Eigen::Matrix<T,3,1> invCross(const Eigen::Matrix<T,3,3>& wCross) {
  return Eigen::Matrix<T,3,1>(wCross(2,1),wCross(0,2),wCross(1,0));
}
template <typename T>   //trace([w]*m) = w.dot(invCrossMatTrace(m))
Eigen::Matrix<T,3,1> invCrossMatTrace(const Eigen::Matrix<T,3,3>& m) {
  return Eigen::Matrix<T,3,1>(m(1,2)-m(2,1),m(2,0)-m(0,2),m(0,1)-m(1,0));
}
template <typename T>   //trace([wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(m)*wB)
Eigen::Matrix<T,3,3> invDoubleCrossMatTrace(const Eigen::Matrix<T,3,3>& m) {
  Eigen::Matrix<T,3,3> ret;
  ret <<
      -m(1,1)-m(2,2),m(1,0),m(2,0),
      m(0,1),-m(0,0)-m(2,2),m(2,1),
      m(0,2),m(1,2),-m(0,0)-m(1,1);
  return ret;
}
template <typename T>   //trace([l]*[wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(l,m)*wB)
Eigen::Matrix<T,3,3> invDoubleCrossMatTrace(const Eigen::Matrix<T,3,1>& l,const Eigen::Matrix<T,3,3>& m) {
  Eigen::Matrix<T,3,3> ret;
  ret<<
     l[0]*m(2,1)-l[0]*m(1,2),
     -l[2]*m(2,2)-l[0]*m(2,0)-l[1]*m(1,2),
     l[2]*m(2,1)+l[1]*m(1,1)+l[0]*m(1,0),

     l[2]*m(2,2)+l[1]*m(2,1)+l[0]*m(0,2),
     m(0,2)*l[1]-l[1]*m(2,0),
     -l[2]*m(2,0)-m(0,1)*l[1]-l[0]*m(0,0),

     -m(1,2)*l[2]-l[1]*m(1,1)-l[0]*m(0,1),
     m(0,2)*l[2]+l[1]*m(1,0)+l[0]*m(0,0),
     m(1,0)*l[2]-m(0,1)*l[2];
  return ret;
}
//spatial
template <typename T>
Eigen::Matrix<T,6,6> toSpatial(const Eigen::Matrix<T,3,4>& t) {
  Eigen::Matrix<T,6,6> ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=t.template block<3,3>(0,0);
  ret.template block<3,3>(3,3)=t.template block<3,3>(0,0);
  ret.template block<3,3>(3,0)=cross<T>(t.template block<3,1>(0,3))*t.template block<3,3>(0,0);
  return ret;
}
template <typename T>
Eigen::Matrix<T,3,4> fromSpatial(const Eigen::Matrix<T,6,6>& t) {
  Eigen::Matrix<T,3,4> ret;
  ret.template block<3,3>(0,0)=t.template block<3,3>(0,0);
  ret.template block<3,1>(0,3)=invCross<T>(t.template block<3,3>(3,0)*t.template block<3,3>(0,0).transpose());
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,6> spatialCross(const Eigen::Matrix<T,6,1>& v) {
  Eigen::Matrix<T,6,6> ret;
  ret.setZero();
  ret.template block<3,3>(3,0)=cross<T>(v.template segment<3>(3));
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=cross<T>(v.template segment<3>(0));
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,1> spatialCross(const Eigen::Matrix<T,6,1>& v,const Eigen::Matrix<T,6,1>& w) {
  Eigen::Matrix<T,6,1> ret;
  ret.template segment<3>(0)=v.template segment<3>(0).cross(w.template segment<3>(0));
  ret.template segment<3>(3)=v.template segment<3>(0).cross(w.template segment<3>(3))+v.template segment<3>(3).cross(w.template segment<3>(0));
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,6> spatialCrossStar(const Eigen::Matrix<T,6,1>& v) {
  Eigen::Matrix<T,6,6> ret;
  ret.setZero();
  ret.template block<3,3>(0,3)=cross<T>(v.template segment<3>(3));
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=cross<T>(v.template segment<3>(0));
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,1> spatialCrossStar(const Eigen::Matrix<T,6,1>& v,const Eigen::Matrix<T,6,1>& w) {
  Eigen::Matrix<T,6,1> ret;
  ret.template segment<3>(0)=v.template segment<3>(0).cross(w.template segment<3>(0))+v.template segment<3>(3).cross(w.template segment<3>(3));
  ret.template segment<3>(3)=v.template segment<3>(0).cross(w.template segment<3>(3));
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,6> spatialXStar(const Eigen::Matrix<T,6,6>& X) {
  Eigen::Matrix<T,6,6> ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=X.template block<3,3>(0,0);
  ret.template block<3,3>(0,3)=-X.template block<3,3>(0,0)*X.template block<3,3>(3,0).transpose()*X.template block<3,3>(0,0);
  return ret;
}
template <typename T>
Eigen::Matrix<T,6,6> spatialInv(const Eigen::Matrix<T,6,6>& X) {
  Eigen::Matrix<T,6,6> ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=X.template block<3,3>(0,0).transpose();
  ret.template block<3,3>(3,0)=-X.template block<3,3>(0,0).transpose()*X.template block<3,3>(3,0)*X.template block<3,3>(0,0).transpose();
  return ret;
}
template <typename T>
Eigen::Matrix<T,3,1> spatialVel(const Eigen::Matrix<T,6,1>& V,const Eigen::Matrix<T,3,1>& v) {
  return V.template segment<3>(0).cross(v)+V.template segment<3>(3);
}
template <typename T>
Eigen::Matrix<T,3,1> applyTrans(const Eigen::Matrix<T,3,4>& X,const Eigen::Matrix<T,3,1>& v) {
  return X.template block<3,3>(0,0)*v+X.template block<3,1>(0,3);
}
template <typename T>
Eigen::Matrix<T,3,1> applyTransInv(const Eigen::Matrix<T,3,4>& X,const Eigen::Matrix<T,3,1>& v) {
  return X.template block<3,3>(0,0).transpose()*(v-X.template block<3,1>(0,3));
}
template <typename T>
Eigen::Matrix<T,3,1> spatialApplyTrans(const Eigen::Matrix<T,6,6>& X,const Eigen::Matrix<T,3,1>& v) {
  return applyTrans<T>(fromSpatial<T>(X),v);
}
template <typename T>
Eigen::Matrix<T,3,1> spatialApplyTransInv(const Eigen::Matrix<T,6,6>& X,const Eigen::Matrix<T,3,1>& v) {
  return applyTrans<T>(fromSpatial<T>(spatialInv<T>(X)),v);
}
}

#endif
