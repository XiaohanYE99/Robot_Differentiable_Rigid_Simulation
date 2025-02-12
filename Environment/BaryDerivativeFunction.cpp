#include "BaryDerivativeFunction.h"

namespace PHYSICSMOTION {
template <typename T>
void debugBaryDerivativePointLine() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  char feat;
  Eigen::Matrix<T,2,3> dbdpt;
  Eigen::Matrix<T,2,3> dbdv[2];
  Eigen::Matrix<T,3,1> d,pt,pt2;
  Eigen::Matrix<T,3,1> v[2],v2[2];
  Eigen::Matrix<T,2,1> bary,bary2;
  Eigen::Matrix<T,3,1> cp;
  for(int pass=0; pass<3; pass++)
    while(true) {
      pt.setRandom();
      v[0].setRandom();
      v[1].setRandom();
      distToSqrLineSegment(pt,v,Eigen::Map<Eigen::Matrix<T,2,1>>(bary.data()),cp,&feat);
      distToSqrLineSegment(dbdpt,dbdv,pt,v,bary,feat);
      if( (pass==0 && feat==-1) ||
          (pass==1 && feat==0) ||
          (pass==2 && feat==1)) {
        d.setRandom();
        pt2=pt+d*DELTA;
        distToSqrLineSegment(pt2,v,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdpt",(dbdpt*d).norm(),(dbdpt*d-(bary2-bary)/DELTA).norm())
        v2[0]=v[0]+d*DELTA;
        v2[1]=v[1];
        distToSqrLineSegment(pt,v2,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdv[0]",(dbdv[0]*d).norm(),(dbdv[0]*d-(bary2-bary)/DELTA).norm())
        v2[0]=v[0];
        v2[1]=v[1]+d*DELTA;
        distToSqrLineSegment(pt,v2,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdv[1]",(dbdv[1]*d).norm(),(dbdv[1]*d-(bary2-bary)/DELTA).norm())
        break;
      }
    }
}
template <typename T>
void debugBaryDerivativeLineLine() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  Eigen::Matrix<char,2,1> feat;
  Eigen::Matrix<T,2,3> dbda[2],dbdb[2];
  Eigen::Matrix<T,3,1> d,a[2],b[2],a2[2],b2[2];
  Eigen::Matrix<T,2,1> bary,bary2;
  Eigen::Matrix<T,3,1> cpa,cpb;
  for(int pass=0; pass<4; pass++)
    while(true) {
      a[0].setRandom();
      a[1].setRandom();
      b[0].setRandom();
      b[1].setRandom();
      distToSqrLineSegment(a,b,Eigen::Map<Eigen::Matrix<T,2,1>>(bary.data()),cpa,cpb,&feat);
      distToSqrLineSegment(dbda,dbdb,a,b,bary,feat);
      if( (pass==0 && feat[0]==-1 && feat[1]==-1) ||
          (pass==1 && feat[0]!=-1 && feat[1]==-1) ||
          (pass==2 && feat[0]==-1 && feat[1]!=-1) ||
          (pass==3 && feat[0]!=-1 && feat[1]!=-1)) {
        d.setRandom();
        a2[0]=a[0]+d*DELTA;
        a2[1]=a[1];
        distToSqrLineSegment(a2,b,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cpa,cpb,&feat);
        DEBUG_GRADIENT("dbda[0]",(dbda[0]*d).norm(),(dbda[0]*d-(bary2-bary)/DELTA).norm())
        a2[0]=a[0];
        a2[1]=a[1]+d*DELTA;
        distToSqrLineSegment(a2,b,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cpa,cpb,&feat);
        DEBUG_GRADIENT("dbda[1]",(dbda[1]*d).norm(),(dbda[1]*d-(bary2-bary)/DELTA).norm())
        b2[0]=b[0]+d*DELTA;
        b2[1]=b[1];
        distToSqrLineSegment(a,b2,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cpa,cpb,&feat);
        DEBUG_GRADIENT("dbdb[0]",(dbdb[0]*d).norm(),(dbdb[0]*d-(bary2-bary)/DELTA).norm())
        b2[0]=b[0];
        b2[1]=b[1]+d*DELTA;
        distToSqrLineSegment(a,b2,Eigen::Map<Eigen::Matrix<T,2,1>>(bary2.data()),cpa,cpb,&feat);
        DEBUG_GRADIENT("dbdb[1]",(dbdb[1]*d).norm(),(dbdb[1]*d-(bary2-bary)/DELTA).norm())
        break;
      }
    }
}
template <typename T>
void debugBaryDerivativePointTriangle() {
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  Eigen::Matrix<char,2,1> feat;
  Eigen::Matrix<T,3,3> dbdpt;
  Eigen::Matrix<T,3,3> dbdv[3];
  Eigen::Matrix<T,3,1> d,pt,pt2;
  Eigen::Matrix<T,3,1> v[3],v2[3];
  Eigen::Matrix<T,3,1> bary,bary2;
  Eigen::Matrix<T,3,1> cp;
  for(int pass=0; pass<3; pass++)
    while(true) {
      pt.setRandom();
      v[0].setRandom();
      v[1].setRandom();
      v[2].setRandom();
      distToSqrTriangle(pt,v,Eigen::Map<Eigen::Matrix<T,3,1>>(bary.data()),cp,&feat);
      distToSqrTriangle(dbdpt,dbdv,pt,v,bary,feat);
      if( (pass==0 && feat[0]==-1 && feat[1]==-1) ||
          (pass==1 && feat[0]!=-1 && feat[1]==-1) ||
          (pass==2 && feat[0]!=-1 && feat[1]!=-1)) {
        d.setRandom();
        pt2=pt+d*DELTA;
        distToSqrTriangle(pt2,v,Eigen::Map<Eigen::Matrix<T,3,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdpt",(dbdpt*d).norm(),(dbdpt*d-(bary2-bary)/DELTA).norm())
        v2[0]=v[0]+d*DELTA;
        v2[1]=v[1];
        v2[2]=v[2];
        distToSqrTriangle(pt,v2,Eigen::Map<Eigen::Matrix<T,3,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdv[0]",(dbdv[0]*d).norm(),(dbdv[0]*d-(bary2-bary)/DELTA).norm())
        v2[0]=v[0];
        v2[1]=v[1]+d*DELTA;
        v2[2]=v[2];
        distToSqrTriangle(pt,v2,Eigen::Map<Eigen::Matrix<T,3,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdv[1]",(dbdv[1]*d).norm(),(dbdv[1]*d-(bary2-bary)/DELTA).norm())
        v2[0]=v[0];
        v2[1]=v[1];
        v2[2]=v[2]+d*DELTA;
        distToSqrTriangle(pt,v2,Eigen::Map<Eigen::Matrix<T,3,1>>(bary2.data()),cp,&feat);
        DEBUG_GRADIENT("dbdv[2]",(dbdv[2]*d).norm(),(dbdv[2]*d-(bary2-bary)/DELTA).norm())
        break;
      }
    }
}
void debugBaryDerivativeAll() {
  debugBaryDerivativePointLine<FLOAT>();
  debugBaryDerivativeLineLine<FLOAT>();
  debugBaryDerivativePointTriangle<FLOAT>();
}
}
