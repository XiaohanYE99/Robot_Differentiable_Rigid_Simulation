#ifndef ROTATION_UTILS_H
#define ROTATION_UTILS_H

#include "Epsilon.h"
#include "CrossSpatialUtils.h"

namespace PHYSICSMOTION {
template <typename T>
struct SinCosTraits {
  static void sincos(T theta,T* a,T* b) {
    *a=sin(theta);
    *b=cos(theta);
  }
};
template <>
struct SinCosTraits<double> {
  static void sincos(double theta,double* a,double* b) {
#ifdef __APPLE__
    __sincos(theta,a,b);
#elif defined(WIN32)
    *a=sin(theta);
    *b=cos(theta);
#else
    ::sincos(theta,a,b);
#endif
  }
};
template <>
struct SinCosTraits<float> {
  static void sincos(float theta,float* a,float* b) {
#if defined(__APPLE__)
    __sincosf(theta,a,b);
#elif defined(WIN32)
    *a=sinf(theta);
    *b=cosf(theta);
#else
    ::sincosf(theta,a,b);
#endif
  }
};
//rotation: expWZ
template <typename T,typename T2>   //hinge joint
Eigen::Matrix<T,3,3> expWZ(T z,T2* DRDZ) {
  T sz,cz;
  SinCosTraits<T>::sincos(z,&sz,&cz);

  Eigen::Matrix<T,3,3> ret;
  ret <<
      cz,-sz,0,
      sz, cz,0,
      0, 0,1;

  if(DRDZ)
    *DRDZ <<   0, 0,1;
  return ret;
}
//rotation: expWYZ
template <typename T,typename T2>   //ball joint
Eigen::Matrix<T,3,3> expWYZ(T y,T z,T2* DRDY,T2* DRDZ,T2* DRDYDZ) {
  T sy,cy,sz,cz;
  SinCosTraits<T>::sincos(y,&sy,&cy);
  SinCosTraits<T>::sincos(z,&sz,&cz);

  Eigen::Matrix<T,3,3> ret;
  ret <<
      cy*cz,-sz,sy*cz,
      cy*sz, cz,sy*sz,
      -sy,  0,   cy;

  if(DRDY)
    *DRDY << -sz,cz,0;
  if(DRDZ)
    *DRDZ <<   0, 0,1;

  if(DRDYDZ)
    *DRDYDZ << -cz,-sz,0;
  return ret;
}
template <typename T,typename T2>
void expWYZLambda(T z,T lambdaZ,T2& DRDYDZ) {
  T sz,cz;
  SinCosTraits<T>::sincos(z,&sz,&cz);
  DRDYDZ << sz*lambdaZ,-cz*lambdaZ,0;
}
//rotation: eulerXYZ
template <typename T>
Eigen::Matrix<T,3,1> invEulerX1Y3Z2(const Eigen::Matrix<T,3,3>& R) {
  Eigen::Matrix<T,3,1> ret;
  T cz=sqrt(R(1,1)*R(1,1)+R(1,2)*R(1,2));
  if(cz<1e-6) {
    ret[2]=M_PI/2;
    ret[0]=ret[1]=atan2(R(0,2),-R(0,1))/2;
  } else {
    ret[2]=atan2(R(1,0),cz);
    ret[1]=atan2(-R(2,0),R(0,0));
    ret[0]=atan2(-R(1,2),R(1,1));
  }
  return ret;
}
template <typename T,typename T2>   //rotation 3D euler angle
Eigen::Matrix<T,3,3> eulerX1Y3Z2(const Eigen::Matrix<T,3,1>& w,T2 diffV[3],T2 ddiffV[9]) {
  //this is a special type of Euler-Angle in favor of our definition in RootInvariantMLP
  T sx,cx,sy,cy,sz,cz;
  SinCosTraits<T>::sincos(w[0],&sx,&cx);
  SinCosTraits<T>::sincos(w[1],&sy,&cy);
  SinCosTraits<T>::sincos(w[2],&sz,&cz);

  Eigen::Matrix<T,3,3> ret;
  ret(0,0)=cy*cz;
  ret(0,1)=sx*sy-cx*cy*sz;  //-cos(x+y) if sz=1
  ret(0,2)=sx*cy*sz+cx*sy;  //sin(x+y) if sz=1
  ret(1,0)=sz;
  ret(1,1)=cx*cz;
  ret(1,2)=-sx*cz;
  ret(2,0)=-sy*cz;
  ret(2,1)=cx*sy*sz+sx*cy;
  ret(2,2)=cx*cy-sx*sy*sz;

  if(diffV) {
    diffV[0] << cy*cz,sz,-sy*cz;
    diffV[1] <<     0, 1,     0;
    diffV[2] <<    sy, 0,    cy;
  }

  if(ddiffV) {
    ddiffV[0] <<      0, 0,     0;
    ddiffV[1] << -sy*cz, 0,-cy*cz;
    ddiffV[2] << -cy*sz,cz, sy*sz;

    ddiffV[3] << 0,0,0;
    ddiffV[4] << 0,0,0;
    ddiffV[5] << 0,0,0;

    ddiffV[6] <<  0,0,  0;
    ddiffV[7] << cy,0,-sy;
    ddiffV[8] <<  0,0,  0;
  }
  return ret;
}
template <typename T,typename T2>
void eulerX1Y3Z2Lambda(const Eigen::Matrix<T,3,1>& w,const Eigen::Matrix<T,3,1>& lambda,T2 ddiffV[9]) {
  T sx,cx,sy,cy,sz,cz;
  SinCosTraits<T>::sincos(w[0],&sx,&cx);
  SinCosTraits<T>::sincos(w[1],&sy,&cy);
  SinCosTraits<T>::sincos(w[2],&sz,&cz);

  ddiffV[0] <<                                0,            0,                              0;
  ddiffV[1] << -cy*cz*lambda[1]+sy*sz*lambda[2],            0,sy*cz*lambda[1]+cy*sz*lambda[2];
  ddiffV[2] <<  sy*sz*lambda[1]-cy*cz*lambda[2],-sz*lambda[2],cy*sz*lambda[1]+sy*cz*lambda[2];

  ddiffV[3] << 0,0,0;
  ddiffV[4] << 0,0,0;
  ddiffV[5] << 0,0,0;

  ddiffV[6] <<             0,0,            0;
  ddiffV[7] << -sy*lambda[1],0,-cy*lambda[1];
  ddiffV[8] <<             0,0,            0;
}
template <typename T>   //rotation 3D rodrigues, helper for diffV
T f_tMst_t3(T t,T st,T ct,T eps) {
  T t2=t*t,t3=t*t2;
  if(fabs(t)<eps)
    return (1037836800+t2*(t2*(1235520+t2*(t2*(156-t2)-17160))-51891840))/6227020800;
  else return (t-st)/t3;
}
template <typename T>
T f_1Mct_t2(T t,T st,T ct,T eps) {
  T t2=t*t;
  if(fabs(t)<eps)
    return (239500800+t2*(t2*(665280+t2*(t2*(132-t2)-11880))-19958400))/479001600;
  else return (1-ct)/t2;
}
template <typename T>
T f_st_t(T t,T st,T ct,T eps) {
  T t2=t*t;
  if(fabs(t)<eps)
    return (39916800+t2*(t2*(332640+t2*(t2*(110-t2)-7920))-6652800))/39916800;
  else return st/t;
}
template <typename T>   //helper for ddiffV
T f_3stM2tMtct_t5(T t,T st,T ct,T eps) {
  T t2=t*t,t5=t2*t2*t;
  if(fabs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-175)+21840)-1801800)+86486400)-1816214400)/108972864000;
  else return (3*st-2*t-t*ct)/t5;
}
template <typename T>
T f_tstP2ctM2_t4(T t,T st,T ct,T eps) {
  T t2=t*t,t4=t2*t2;
  if(fabs(t)<eps)
    return (t2*(t2*(t2*(t2*(3*t2-455)+48048)-3243240)+121080960)-1816214400)/21794572800;
  else return (-2+2*ct+t*st)/t4;
}
template <typename T>
T f_stMtct_t3(T t,T st,T ct,T eps) {
  T t2=t*t,t3=t*t2;
  if(fabs(t)<eps)
    return (172972800+t2*(t2*(617760+t2*(t2*(130-t2)-11440))-17297280))/518918400;
  else return (st-t*ct)/t3;
}
template <typename T>   //helper for dddiffVLambda
T fD_tMst_t3(T t,T st,T ct,T eps) {
  T t2=t*t,t5=t2*t2*t;
  if(fabs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-175)+21840)-1801800)+86486400)-1816214400)/108972864000;
  else return (-2*t-t*ct+3*st)/t5;
}
template <typename T>
T fD_1Mct_t2(T t,T st,T ct,T eps) {
  T t2=t*t,t4=t2*t2;
  if(fabs(t)<eps)
    return (t2*(t2*(t2*(t2*(3*t2-455)+48048)-3243240)+121080960)-1816214400)/21794572800;
  else return (-2+2*ct+t*st)/t4;
}
template <typename T>
T fD_3stM2tMtct_t5(T t,T st,T ct,T eps) {
  T t2=t*t,t7=t2*t2*t2*t;
  if(fabs(t)<eps)
    return (23524300800+t2*(t2*(17821440+t2*(t2*(1360-7*t2)-190400))-980179200))/14820309504000;
  else return (8*t+7*t*ct+(t2-15)*st)/t7;
}
template <typename T>
T fD_tstP2ctM2_t4(T t,T st,T ct,T eps) {
  T t2=t*t,t6=t2*t2*t2;
  if(fabs(t)<eps)
    return (9686476800+t2*(t2*(11531520+t2*(t2*(1200-7*t2)-145600))-518918400))/871782912000;
  else return -(-8+(8-t2)*ct+5*t*st)/t6;
}
template <typename T>
T fD_stMtct_t3(T t,T st,T ct,T eps) {
  T t2=t*t,t5=t2*t2*t;
  if(fabs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-150)+15600)-1029600)+37065600)-518918400)/7783776000;
  else return (3*t*ct+(t2-3)*st)/t5;
}
//rotation: expWGradV
template <typename T,typename T2>
Eigen::Matrix<T,3,3> expWGradV(const Eigen::Matrix<T,3,1>& w,T2* diffV=NULL,T2* ddiffV=NULL) {
  typedef Eigen::Matrix<T,3,3> MAT3;
  typedef Eigen::Matrix<T,3,1> VEC3;
  T epsT=1E-1f;
  T t=w.norm(),st,ct;
  SinCosTraits<T>::sincos(t,&st,&ct);
  MAT3 cw=cross<T>(w);
  //W
  T _tMst_t3=f_tMst_t3<T>(t,st,ct,epsT);
  T _1Mct_t2=f_1Mct_t2<T>(t,st,ct,epsT);
  T _st_t=f_st_t<T>(t,st,ct,epsT);
  if(diffV) {
#define DIFFV(I) diffV[I]=_tMst_t3*w[I]*w+_1Mct_t2*w.cross(VEC3::Unit(I))+_st_t*VEC3::Unit(I);
    DIFFV(0)
    DIFFV(1)
    DIFFV(2)
#undef DIFFV
  }
  //DW
  if(ddiffV) {
    T _3stM2tMtct_t5=f_3stM2tMtct_t5<T>(t,st,ct,epsT);
    T _tstP2ctM2_t4=f_tstP2ctM2_t4<T>(t,st,ct,epsT);
    T _stMtct_t3=f_stMtct_t3<T>(t,st,ct,epsT);
#define DDIFFV(I,J) ddiffV[I*3+J]=  \
    _tMst_t3*(VEC3::Unit(J)*w[I]+(I==J?1:0)*w)+ \
    _3stM2tMtct_t5*(w*w[I]*w[J])+   \
    _1Mct_t2*(VEC3::Unit(J)).cross(VEC3::Unit(I))+  \
    _tstP2ctM2_t4*w.cross(VEC3::Unit(I))*w[J]-  \
    _stMtct_t3*VEC3::Unit(I)*w[J];
    DDIFFV(0,0)
    DDIFFV(0,1)
    DDIFFV(0,2)
    DDIFFV(1,0)
    DDIFFV(1,1)
    DDIFFV(1,2)
    DDIFFV(2,0)
    DDIFFV(2,1)
    DDIFFV(2,2)
#undef DDIFFV
  }
  //R
  return MAT3::Identity()+cw*_st_t+cw*cw*_1Mct_t2;
}
template <typename T,typename T2>
void expWGradVLambda(const Eigen::Matrix<T,3,1>& w,const Eigen::Matrix<T,3,1>& lambda,T2* ddiffV) {
  typedef Eigen::Matrix<T,3,1> VEC3;
  T epsT=1E-1f;
  T t=w.norm(),st,ct,dwl=w.dot(lambda);
  SinCosTraits<T>::sincos(t,&st,&ct);
  //W
  T _tMst_t3=f_tMst_t3<T>(t,st,ct,epsT);
  T _3stM2tMtct_t5=f_3stM2tMtct_t5<T>(t,st,ct,epsT);
  T _tstP2ctM2_t4=f_tstP2ctM2_t4<T>(t,st,ct,epsT);
  T _stMtct_t3=f_stMtct_t3<T>(t,st,ct,epsT);

  T D_tMst_t3=fD_tMst_t3<T>(t,st,ct,epsT)*dwl;
  T D_1Mct_t2=fD_1Mct_t2<T>(t,st,ct,epsT)*dwl;
  T D_3stM2tMtct_t5=fD_3stM2tMtct_t5<T>(t,st,ct,epsT)*dwl;
  T D_tstP2ctM2_t4=fD_tstP2ctM2_t4<T>(t,st,ct,epsT)*dwl;
  T D_stMtct_t3=fD_stMtct_t3<T>(t,st,ct,epsT)*dwl;
#define DDIFFV(I,J) ddiffV[I*3+J]=  \
  _tMst_t3*(VEC3::Unit(J)*lambda[I]+(I==J?1:0)*lambda)+ \
  _3stM2tMtct_t5*(lambda*w[I]*w[J]+w*lambda[I]*w[J]+w*w[I]*lambda[J])+   \
  _tstP2ctM2_t4*(lambda.cross(VEC3::Unit(I))*w[J]+w.cross(VEC3::Unit(I))*lambda[J])-  \
  _stMtct_t3*VEC3::Unit(I)*lambda[J]+   \
    \
  D_tMst_t3*(VEC3::Unit(J)*w[I]+(I==J?1:0)*w)+ \
  D_3stM2tMtct_t5*(w*w[I]*w[J])+   \
  D_1Mct_t2*(VEC3::Unit(J)).cross(VEC3::Unit(I))+  \
  D_tstP2ctM2_t4*w.cross(VEC3::Unit(I))*w[J]-  \
  D_stMtct_t3*VEC3::Unit(I)*w[J];
  DDIFFV(0,0)
  DDIFFV(0,1)
  DDIFFV(0,2)
  DDIFFV(1,0)
  DDIFFV(1,1)
  DDIFFV(1,2)
  DDIFFV(2,0)
  DDIFFV(2,1)
  DDIFFV(2,2)
#undef DDIFFV
}
template <typename T>
Eigen::Matrix<T,3,1> invExpW(const Eigen::Matrix<T,3,3>& R) {
  T cosTheta=(R.trace()-1)/2,coef;
  if(cosTheta < Epsilon<T>::rotationEps()-1) {
    int a,b,c;
    R.diagonal().maxCoeff(&a);
    b=(a+1)%3;
    c=(a+2)%3;

    T s=sqrt(R(a,a)-R(b,b)-R(c,c)+1);
    Eigen::Matrix<T,3,1> ret;
    ret[a]=s/2;
    ret[b]=(R(a,b)+R(b,a))/(2*s);
    ret[c]=(R(a,c)+R(c,a))/(2*s);
    return ret/sqrt(ret.squaredNorm())*M_PI;
  } else if(cosTheta > 1-Epsilon<T>::rotationEps()) {
    T theta=acos(std::min<T>(cosTheta,1));
    coef=0.5f+theta*theta/12.0f;  //using taylor expansion
  } else {
    T theta=acos(cosTheta),sinTheta=sin(theta);
    coef=theta/(2*sinTheta);
  }
  return Eigen::Matrix<T,3,1>(R(2,1)-R(1,2),R(0,2)-R(2,0),R(1,0)-R(0,1))*coef;
}
template <typename T>
T getAngle3D(const Eigen::Matrix<T,3,1>& a,const Eigen::Matrix<T,3,1>& b) {
  T denom=std::max<T>(a.norm()*b.norm(),Epsilon<T>::rotationEps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-Epsilon<T>::rotationEps())),(T)(-1.0f+Epsilon<T>::rotationEps())));
}
template <typename T>
T getAngle2D(const Eigen::Matrix<T,2,1>& a,const Eigen::Matrix<T,2,1>& b) {
  T denom=std::max<T>(a.norm()*b.norm(),Epsilon<T>::rotationEps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-Epsilon<T>::rotationEps())),(T)(-1.0f+Epsilon<T>::rotationEps())));
}
}

#endif
