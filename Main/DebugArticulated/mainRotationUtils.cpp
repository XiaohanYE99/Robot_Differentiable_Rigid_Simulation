#include <Utils/CrossSpatialUtils.h>
#include <Utils/SpatialRotationUtils.h>
#include <Utils/RotationUtils.h>
#include <Utils/DebugGradient.h>
#include <random>

using namespace PHYSICSMOTION;

template <typename T>
void debugCross() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  std::cout << "-------------------------------------------------------------DebugCross" << std::endl;
  Vec3T l=Vec3T::Random();
  Vec3T w=Vec3T::Random();
  Vec3T wA=Vec3T::Random();
  Vec3T wB=Vec3T::Random();
  Mat3T m=Mat3T::Random();
  DEBUG_GRADIENT("cross",w.norm(),(w-invCross<T>(cross<T>(w))).norm())
  DEBUG_GRADIENT("invCrossMatTrace",sqrt((cross<T>(w)*m).trace()),sqrt((cross<T>(w)*m).trace()-w.dot(invCrossMatTrace<T>(m))))
  DEBUG_GRADIENT("invDoubleCrossMatTrace",(cross<T>(wA)*m*cross<T>(wB)).trace(),(cross<T>(wA)*m*cross<T>(wB)).trace()-wA.dot(invDoubleCrossMatTrace<T>(m)*wB))
  DEBUG_GRADIENT("invDoubleCrossMatTrace2",(cross<T>(l)*cross<T>(wA)*m*cross<T>(wB)).trace(),(cross<T>(l)*cross<T>(wA)*m*cross<T>(wB)).trace()-wA.dot(invDoubleCrossMatTrace<T>(l,m)*wB))
}
template <typename T>
void debugSpatial() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  std::cout << "-------------------------------------------------------------DebugSpatial" << std::endl;
  Vec4T q;
  Vec6T a,a2;
  Mat3X4T t,t2;
  q.setRandom();
  q.normalize();
  a.setRandom();
  a2.setRandom();
  ROT(t)=QuatT(q[0],q[1],q[2],q[3]).toRotationMatrix();
  CTR(t)=Vec3T::Random();
  Mat6T ST=toSpatial<T>(t);
  t2=fromSpatial<T>(ST);
  Mat6T invST=spatialInv<T>(ST);
  Mat6T invTST=spatialXStar<T>(ST);
  DEBUG_GRADIENT("toFromSpatial",t.norm(),(t-t2).norm())
  DEBUG_GRADIENT("invSpatial",ST.norm(),(ST*invST-Mat6T::Identity()).norm())
  DEBUG_GRADIENT("spatialXStar",invST.norm(),(invST-invTST.transpose()).norm())
  DEBUG_GRADIENT("spatialCross",spatialCross<T>(a,a2).norm(),(spatialCross<T>(a,a2)-spatialCross<T>(a)*a2).norm())
  DEBUG_GRADIENT("spatialCrossStar",spatialCrossStar<T>(a,a2).norm(),(spatialCrossStar<T>(a,a2)-spatialCrossStar<T>(a)*a2).norm())
}
template <typename T>
void debugRotation() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  std::cout << "-------------------------------------------------------------debugRotation" << std::endl;
  Vec3T lambda,DRDY,DRDY2,DRDZ,DRDYDZ,DRDYDZ2,DRDYDZLambda;
  Vec3T w,w2,dw,diffV[3],diffV2[3],ddiffV[9],ddiffV2[9],dddiffVLambda[9],diffVw,ddiffVw;
  Mat3T R,R2;

  //expWZ
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1,1);
  T y=dis(gen);
  T z=dis(gen);
  R=expWZ<T,Vec3T>(z,&DRDZ);
  R2=expWZ<T,Vec3T>(z+DELTA,NULL);
  DEBUG_GRADIENT("expWZ-DZ",(cross<T>(DRDZ)*R).norm(),(cross<T>(DRDZ)*R-(R2-R)/DELTA).norm())

  //expWYZ
  R=expWYZ<T,Vec3T>(y,z,&DRDY,&DRDZ,&DRDYDZ);
  R2=expWYZ<T,Vec3T>(y+DELTA,z,&DRDY2,NULL,NULL);
  DEBUG_GRADIENT("expWYZ-DY",(cross<T>(DRDY)*R).norm(),(cross<T>(DRDY)*R-(R2-R)/DELTA).norm())
  R2=expWYZ<T,Vec3T>(y,z+DELTA,&DRDY2,NULL,&DRDYDZ2);
  DEBUG_GRADIENT("expWYZ-DZ",(cross<T>(DRDZ)*R).norm(),(cross<T>(DRDZ)*R-(R2-R)/DELTA).norm())
  DEBUG_GRADIENT("expWYZ-DYDZ",(DRDYDZ).norm(),((DRDYDZ-(DRDY2-DRDY)/DELTA).norm()))
  expWYZLambda<T,Vec3T>(z,1,DRDYDZLambda);
  DEBUG_GRADIENT("expWYZ-DYDZLambda",DRDYDZLambda.norm(),(DRDYDZLambda-(DRDYDZ2-DRDYDZ)/DELTA).norm())

  //eulerX1Y3Z2
  w=Vec3T::Random();
  R=eulerX1Y3Z2<T,Vec3T>(w,NULL,NULL);
  w2=invEulerX1Y3Z2<T>(R);
  DEBUG_GRADIENT("invEulerXYZ",w.norm(),(w-w2).norm())
  w[0]=w[1];
  w[2]=M_PI/2;
  R=eulerX1Y3Z2<T,Vec3T>(w,NULL,NULL);
  w2=invEulerX1Y3Z2<T>(R);
  DEBUG_GRADIENT("invEulerXYZ",w.norm(),(w-w2).norm())

  //expWXYZ
  w=Vec3T::Random();
  dw=Vec3T::Random();
  lambda=Vec3T::Random();
  R=eulerX1Y3Z2<T,Vec3T>(w,diffV,ddiffV);
  R2=eulerX1Y3Z2<T,Vec3T>(w+dw*DELTA,diffV2,ddiffV2);
  diffVw=diffV[0]*dw[0]+diffV[1]*dw[1]+diffV[2]*dw[2];
  DEBUG_GRADIENT("eulerXYZ-diffV",(cross<T>(diffVw)*R).norm(),(cross<T>(diffVw)*R-(R2-R)/DELTA).norm())
  for(int d=0; d<3; d++) {
    ddiffVw=ddiffV[d*3+0]*dw[0]+ddiffV[d*3+1]*dw[1]+ddiffV[d*3+2]*dw[2];
    DEBUG_GRADIENT("eulerXYZ-ddiffV"+std::to_string(d),ddiffVw.norm(),((diffV2[d]-diffV[d])/DELTA-ddiffVw).norm())
  }
  //expWXYZ-Lambda
  eulerX1Y3Z2Lambda<T,Vec3T>(w,dw,dddiffVLambda);
  for(int d=0; d<9; d++) {
    DEBUG_GRADIENT("eulerXYZ-dddiffVLambda",dddiffVLambda[d].norm(),(dddiffVLambda[d]-(ddiffV2[d]-ddiffV[d])/DELTA).norm())
  }

  //expW
  R=expWGradV<T,Vec3T>(w,diffV,ddiffV);
  w2=invExpW<T>(R);
  DEBUG_GRADIENT("expW-inv",w.norm(),(w-w2).norm())
  R2=expWGradV<T,Vec3T>(w+dw*DELTA,diffV2,ddiffV2);
  diffVw=diffV[0]*dw[0]+diffV[1]*dw[1]+diffV[2]*dw[2];
  DEBUG_GRADIENT("expW-diffV",(cross<T>(diffVw)*R).norm(),(cross<T>(diffVw)*R-(R2-R)/DELTA).norm())
  for(int d=0; d<3; d++) {
    ddiffVw=ddiffV[d*3+0]*dw[0]+ddiffV[d*3+1]*dw[1]+ddiffV[d*3+2]*dw[2];
    DEBUG_GRADIENT("expW-ddiffV"+std::to_string(d),ddiffVw.norm(),((diffV2[d]-diffV[d])/DELTA-ddiffVw).norm())
  }
  //expW-Lambda
  expWGradVLambda<T,Vec3T>(w,dw,dddiffVLambda);
  for(int d=0; d<9; d++) {
    DEBUG_GRADIENT("expW-dddiffVLambda",dddiffVLambda[d].norm(),(dddiffVLambda[d]-(ddiffV2[d]-ddiffV[d])/DELTA).norm())
  }
}
template <typename T>
void debugSpatialRotation() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  std::cout << "-------------------------------------------------------------debugSpatialRotation" << std::endl;
  Mat6XT S,S2,DvJDq,DdvJDq,DdvJDdq,DSTDqa;
  Vec6T vJ,vJ2,dvJ,dvJ2,f,Sa,Sa2;
  Mat6T X,X2;
  Mat3T R;
  f.setRandom();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1,1);
  for(int pass=0; pass<2; pass++) {
    {
      T q,dq,dq2,ddq,deltaq,STf,STf2,a=dis(gen);
      S.resize(6,1);
      S2.resize(6,1);
      DvJDq.resize(6,1);
      DdvJDq.resize(6,1);
      DdvJDdq.resize(6,1);
      q=dis(gen);
      dq=dis(gen);
      ddq=dis(gen);
      deltaq=dis(gen);
      if(pass==1)
        ddq=0;
      R=spatialExpWZ<T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&dq,pass==0?&ddq:NULL,&X,&S,&vJ,&dvJ,&DvJDq,&DdvJDq,&DdvJDdq);
      STf=S.col(0).dot(f);
      Sa=S.col(0)*a;
      DEBUG_GRADIENT("spatialExpWZ",sqrt(R.squaredNorm()),sqrt((R-expWZ<T,Vec3T>(q,NULL)).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWZ-S",sqrt((S*dq).squaredNorm()),sqrt((S*dq-vJ).squaredNorm()))
      spatialExpWZ<T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+dq*DELTA,&(dq2=dq+ddq*DELTA),pass==0?&ddq:NULL,&X2,&S2,&vJ2);
      STf2=S2.col(0).dot(f);
      Sa2=S2.col(0)*a;
      DEBUG_GRADIENT("spatialExpWZ-STf",0,(STf2-STf)/DELTA)
      DEBUG_GRADIENT("spatialExpWZ-Sa",0,sqrt(((Sa2-Sa)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWZ-dvJ",sqrt((X*dvJ).squaredNorm()),sqrt((X*dvJ-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      spatialExpWZ<T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+deltaq*DELTA,&dq,pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpWZ-DvJDq",sqrt((X*DvJDq*deltaq).squaredNorm()),sqrt((X*DvJDq*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWZ-DdvJDq",sqrt((X*DdvJDq*deltaq).squaredNorm()),sqrt((X*DdvJDq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
      spatialExpWZ<T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&(dq2=dq+deltaq*DELTA),pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpWZ-DvJDdq",sqrt((X*S*deltaq).squaredNorm()),sqrt((X*S*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWZ-DdvJDdq",sqrt((X*DdvJDdq*deltaq).squaredNorm()),sqrt((X*DdvJDdq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
    }
    {
      Vec2T q,dq,dq2,ddq,deltaq,STf,STf2,a;
      Mat2T DSTDqTf=Mat2T::Zero();
      S.resize(6,2);
      S2.resize(6,2);
      DvJDq.resize(6,2);
      DdvJDq.resize(6,2);
      DdvJDdq.resize(6,2);
      DSTDqa.setZero(6,2);
      q.setRandom();
      dq.setRandom();
      ddq.setRandom();
      deltaq.setRandom();
      a.setRandom();
      if(pass==1)
        ddq.setZero();
      R=spatialExpWYZ<T,Vec2T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&dq,pass==0?&ddq:NULL,&X,&S,&vJ,&dvJ,&DvJDq,&DdvJDq,&DdvJDdq);
      spatialExpWYZDSTDqf<T,Vec2T,Mat2T,Vec6T>(q,DSTDqTf,f);
      spatialExpWYZDSDqf<T,Vec2T,Mat6XT,Vec2T>(q,DSTDqa,a);
      STf=S.transpose()*f;
      Sa=S*a;
      DEBUG_GRADIENT("spatialExpWYZ",sqrt(R.squaredNorm()),sqrt((R-expWYZ<T,Vec3T>(q[0],q[1],NULL,NULL,NULL)).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWYZ-S",sqrt((S*dq).squaredNorm()),sqrt((S*dq-vJ).squaredNorm()))
      spatialExpWYZ<T,Vec2T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+dq*DELTA,&(dq2=dq+ddq*DELTA),pass==0?&ddq:NULL,&X2,&S2,&vJ2);
      STf2=S2.transpose()*f;
      Sa2=S2*a;
      DEBUG_GRADIENT("spatialExpWYZ-STf",sqrt((DSTDqTf*dq).squaredNorm()),sqrt((DSTDqTf*dq-(STf2-STf)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWYZ-Sa",sqrt((DSTDqa*dq).squaredNorm()),sqrt((DSTDqa*dq-(Sa2-Sa)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWYZ-dvJ",sqrt((X*dvJ).squaredNorm()),sqrt((X*dvJ-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      spatialExpWYZ<T,Vec2T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+deltaq*DELTA,&dq,pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpWYZ-DvJDq",sqrt((X*DvJDq*deltaq).squaredNorm()),sqrt((X*DvJDq*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWYZ-DdvJDq",sqrt((X*DdvJDq*deltaq).squaredNorm()),sqrt((X*DdvJDq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
      spatialExpWYZ<T,Vec2T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&(dq2=dq+deltaq*DELTA),pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpWYZ-DvJDdq",sqrt((X*S*deltaq).squaredNorm()),sqrt((X*S*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpWYZ-DdvJDdq",sqrt((X*DdvJDdq*deltaq).squaredNorm()),sqrt((X*DdvJDdq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
    }
    {
      Vec3T q,dq,dq2,ddq,deltaq,STf,STf2,a;
      Mat3T DSTDqf=Mat3T::Zero();
      S.resize(6,3);
      S2.resize(6,3);
      DvJDq.resize(6,3);
      DdvJDq.resize(6,3);
      DdvJDdq.resize(6,3);
      DSTDqa.setZero(6,3);
      q.setRandom();
      dq.setRandom();
      ddq.setRandom();
      deltaq.setRandom();
      a.setRandom();
      if(pass==1)
        ddq.setZero();
      R=spatialEulerX1Y3Z2<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&dq,pass==0?&ddq:NULL,&X,&S,&vJ,&dvJ,&DvJDq,&DdvJDq,&DdvJDdq);
      spatialEulerX1Y3Z2DSTDqf<T,Vec3T,Mat3T,Vec6T>(q,DSTDqf,f);
      spatialEulerX1Y3Z2DSDqf<T,Vec3T,Mat6XT,Vec3T>(q,DSTDqa,a);
      STf=S.transpose()*f;
      Sa=S*a;
      DEBUG_GRADIENT("spatialEulerXYZ",sqrt(R.squaredNorm()),sqrt((R-eulerX1Y3Z2<T,Vec3T>(q,NULL,NULL)).squaredNorm()))
      DEBUG_GRADIENT("spatialEulerXYZ-S",sqrt((S*dq).squaredNorm()),sqrt((S*dq-vJ).squaredNorm()))
      spatialEulerX1Y3Z2<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+dq*DELTA,&(dq2=dq+ddq*DELTA),pass==0?&ddq:NULL,&X2,&S2,&vJ2);
      STf2=S2.transpose()*f;
      Sa2=S2*a;
      DEBUG_GRADIENT("spatialEulerXYZ-STf",sqrt((DSTDqf*dq).squaredNorm()),sqrt((DSTDqf*dq-(STf2-STf)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialEulerXYZ-Sa",sqrt((DSTDqa*dq).squaredNorm()),sqrt((DSTDqa*dq-(Sa2-Sa)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialEulerXYZ-dvJ",sqrt((X*dvJ).squaredNorm()),sqrt((X*dvJ-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      spatialEulerX1Y3Z2<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+deltaq*DELTA,&dq,pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialEulerXYZ-DvJDq",sqrt((X*DvJDq*deltaq).squaredNorm()),sqrt((X*DvJDq*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialEulerXYZ-DdvJDq",sqrt((X*DdvJDq*deltaq).squaredNorm()),sqrt((X*DdvJDq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
      spatialEulerX1Y3Z2<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&(dq2=dq+deltaq*DELTA),pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialEulerXYZ-DvJDdq",sqrt((X*S*deltaq).squaredNorm()),sqrt((X*S*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialEulerXYZ-DdvJDdq",sqrt((X*DdvJDdq*deltaq).squaredNorm()),sqrt((X*DdvJDdq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
    }
    {
      Vec3T q,dq,dq2,ddq,deltaq,STf,STf2,a;
      Mat3T DSTDqf=Mat3T::Zero();
      S.resize(6,3);
      S2.resize(6,3);
      DvJDq.resize(6,3);
      DdvJDq.resize(6,3);
      DdvJDdq.resize(6,3);
      DSTDqa.setZero(6,3);
      q.setRandom();
      dq.setRandom();
      ddq.setRandom();
      deltaq.setRandom();
      a.setRandom();
      if(pass==1)
        ddq.setZero();
      R=spatialExpW<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&dq,pass==0?&ddq:NULL,&X,&S,&vJ,&dvJ,&DvJDq,&DdvJDq,&DdvJDdq);
      spatialExpWDSTDqf<T,Vec3T,Mat3T,Vec6T>(q,DSTDqf,f);
      spatialExpWDSDqf<T,Vec3T,Mat6XT,Vec3T>(q,DSTDqa,a);
      STf=S.transpose()*f;
      Sa=S*a;
      DEBUG_GRADIENT("spatialExpW",sqrt(R.squaredNorm()),sqrt((R-expWGradV<T,Vec3T>(q,NULL,NULL)).squaredNorm()))
      DEBUG_GRADIENT("spatialExpW-S",sqrt((S*dq).squaredNorm()),sqrt((S*dq-vJ).squaredNorm()))
      spatialExpW<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+dq*DELTA,&(dq2=dq+ddq*DELTA),pass==0?&ddq:NULL,&X2,&S2,&vJ2);
      STf2=S2.transpose()*f;
      Sa2=S2*a;
      DEBUG_GRADIENT("spatialExpW-STf",sqrt((DSTDqf*dq).squaredNorm()),sqrt((DSTDqf*dq-(STf2-STf)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpW-Sa",sqrt((DSTDqa*dq).squaredNorm()),sqrt((DSTDqa*dq-(Sa2-Sa)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpW-dvJ",sqrt((X*dvJ).squaredNorm()),sqrt((X*dvJ-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      spatialExpW<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q+deltaq*DELTA,&dq,pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpW-DvJDq",sqrt((X*DvJDq*deltaq).squaredNorm()),sqrt((X*DvJDq*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpW-DdvJDq",sqrt((X*DdvJDq*deltaq).squaredNorm()),sqrt((X*DdvJDq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
      spatialExpW<T,Vec3T,Vec6T,Mat6T,Mat6XT,Mat6XT>(q,&(dq2=dq+deltaq*DELTA),pass==0?&ddq:NULL,&X2,NULL,&vJ2,&dvJ2);
      DEBUG_GRADIENT("spatialExpW-DvJDdq",sqrt((X*S*deltaq).squaredNorm()),sqrt((X*S*deltaq-(X2*vJ2-X*vJ)/DELTA).squaredNorm()))
      DEBUG_GRADIENT("spatialExpW-DdvJDdq",sqrt((X*DdvJDdq*deltaq).squaredNorm()),sqrt((X*DdvJDdq*deltaq-(X2*dvJ2-X*dvJ)/DELTA).squaredNorm()))
    }
  }
}
int main(int argc,char** argv) {
  mpfr_float::default_precision(100);
  debugCross<float>();
  debugCross<double>();
  debugCross<float128>();
  debugCross<mpfr_float>();

  debugSpatial<float>();
  debugSpatial<double>();
  debugSpatial<float128>();
  debugSpatial<mpfr_float>();

  debugRotation<float>();
  debugRotation<double>();
  debugRotation<float128>();
  debugRotation<mpfr_float>();

  debugSpatialRotation<float>();
  debugSpatialRotation<double>();
  debugSpatialRotation<float128>();
  debugSpatialRotation<mpfr_float>();
  return 0;
}
