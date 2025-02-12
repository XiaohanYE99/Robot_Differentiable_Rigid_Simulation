#ifndef BARY_DERIVATIVE_FUNCTION_H
#define BARY_DERIVATIVE_FUNCTION_H

#include "DistanceFunction.h"

namespace PHYSICSMOTION {
template <typename T>
void distToSqrLineSegment(Eigen::Matrix<T,2,3>& dbdpt,
                          Eigen::Matrix<T,2,3> dbdv[2],
                          const Eigen::Matrix<T,3,1>& pt,
                          const Eigen::Matrix<T,3,1> v[2],
                          const Eigen::Matrix<T,2,1>& bary,
                          char feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  if(feat==0 || feat==1)
    dbdpt=dbdv[0]=dbdv[1]=Eigen::Matrix<T,2,3>::Zero();
  else {
    //feat==-1
    Vec3T dAlpha;
    T alpha=bary[1];
    Vec3T LHSE=v[1]-v[0],RHSE=pt-v[0];
    //pt
    dAlpha=LHSE/LHSE.dot(LHSE);
    dbdpt.row(1).array()=dAlpha.array();
    dbdpt.row(0).array()=-dAlpha.array();
    //dbdv[0]
    dAlpha=(-RHSE-LHSE+2*alpha*LHSE)/LHSE.dot(LHSE);
    dbdv[0].row(1).array()=dAlpha.array();
    dbdv[0].row(0).array()=-dAlpha.array();
    //dbdv[1]
    dAlpha=(RHSE-2*alpha*LHSE)/LHSE.dot(LHSE);
    dbdv[1].row(1).array()=dAlpha.array();
    dbdv[1].row(0).array()=-dAlpha.array();
  }
}
template <typename T>
void distToSqrLineSegment(Eigen::Matrix<T,2,3> dbda[2],
                          Eigen::Matrix<T,2,3> dbdb[2],
                          const Eigen::Matrix<T,3,1> a[2],
                          const Eigen::Matrix<T,3,1> b[2],
                          const Eigen::Matrix<T,2,1>& bary,
                          const Eigen::Matrix<char,2,1>& feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  if(feat[0]!=-1 && feat[1]!=-1)
    dbda[0]=dbda[1]=dbdb[0]=dbdb[1]=Eigen::Matrix<T,2,3>::Zero();
  else if(feat[1]!=-1) {
    Eigen::Matrix<T,2,3> dbdpt;
    dbdb[0]=dbdb[1]=Eigen::Matrix<T,2,3>::Zero();
    distToSqrLineSegment(dbdpt,dbda,b[(int)feat[1]],a,Eigen::Matrix<T,2,1>(1-bary[0],bary[0]),-1);
    dbdb[(int)feat[1]].row(0)=-dbdpt.row(0);
    dbda[0]*=-1;
    dbda[0].row(1).setZero();
    dbda[1]*=-1;
    dbda[1].row(1).setZero();
  } else if(feat[0]!=-1) {
    Eigen::Matrix<T,2,3> dbdpt;
    dbda[0]=dbda[1]=Eigen::Matrix<T,2,3>::Zero();
    distToSqrLineSegment(dbdpt,dbdb,a[(int)feat[0]],b,Eigen::Matrix<T,2,1>(1-bary[1],bary[1]),-1);
    dbda[(int)feat[0]].row(1)=-dbdpt.row(0);
    dbdb[0].row(0).setZero();
    dbdb[1].row(0).setZero();
  } else {
    //feat==-1,-1
    Mat3X2T LHS,dLHS;
    LHS.col(0)= (a[1]-a[0]);
    LHS.col(1)=-(b[1]-b[0]);
    Vec3T RHS=b[0]-a[0],dRHS;
    Mat2T LTL=LHS.transpose()*LHS;
    Mat2T invLTL=LTL.inverse();
    Vec2T alpha=invLTL*(LHS.transpose()*RHS);
    //a[0]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,0)=-1;
      dRHS.setZero();
      dRHS[d]=-1;
      dbda[0].col(d)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS+LHS.transpose()*dRHS);
    }
    //a[1]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,0)=1;
      dbda[1].col(d)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS);
    }
    //b[0]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,1)=1;
      dRHS.setZero();
      dRHS[d]=1;
      dbdb[0].col(d)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS+LHS.transpose()*dRHS);
    }
    //b[1]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,1)=-1;
      dbdb[1].col(d)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS);
    }
  }
}
template <typename T>
void distToSqrTriangle(Eigen::Matrix<T,3,3>& dbdpt,
                       Eigen::Matrix<T,3,3> dbdv[3],
                       const Eigen::Matrix<T,3,1>& pt,
                       const Eigen::Matrix<T,3,1> v[3],
                       const Eigen::Matrix<T,3,1>& bary,
                       const Eigen::Matrix<char,2,1>& feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  if(feat[0]!=-1 && feat[1]==-1) {
    //vertex
    dbdpt=dbdv[0]=dbdv[1]=dbdv[2]=Eigen::Matrix<T,3,3>::Zero();
  } else if(feat[0]!=-1 && feat[1]!=-1) {
    //edge
    Eigen::Matrix<T,2,3> dbdpt2;
    Eigen::Matrix<T,2,3> dbde[2];
    Eigen::Matrix<T,3,1> e[2]= {v[(int)feat[0]],v[(int)feat[1]]};
    dbdpt=dbdv[0]=dbdv[1]=dbdv[2]=Eigen::Matrix<T,3,3>::Zero();
    distToSqrLineSegment(dbdpt2,dbde,pt,e,Eigen::Matrix<T,2,1>(bary[(int)feat[0]],bary[(int)feat[1]]),-1);
    dbdpt.row((int)feat[0])=dbdpt2.row(0);
    dbdpt.row((int)feat[1])=dbdpt2.row(1);
    dbdv[(int)feat[0]].row((int)feat[0])=dbde[0].row(0);
    dbdv[(int)feat[0]].row((int)feat[1])=dbde[0].row(1);
    dbdv[(int)feat[1]].row((int)feat[0])=dbde[1].row(0);
    dbdv[(int)feat[1]].row((int)feat[1])=dbde[1].row(1);
  } else {
    //interior
    Mat3X2T LHS,dLHS;
    LHS.col(0)=v[1]-v[0];
    LHS.col(1)=v[2]-v[0];
    Vec3T RHS=pt-v[0],dRHS;
    Mat2T LTL=LHS.transpose()*LHS;
    Mat2T invLTL=LTL.inverse();
    Vec2T alpha=invLTL*(LHS.transpose()*RHS);
    //v[0]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,0)=-1;
      dLHS(d,1)=-1;
      dRHS.setZero();
      dRHS[d]=-1;
      dbdv[0].col(d).template segment<2>(1)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS+LHS.transpose()*dRHS);
      dbdv[0].col(d)[0]=-dbdv[0].col(d).template segment<2>(1).sum();
    }
    //v[1]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,0)=1;
      dbdv[1].col(d).template segment<2>(1)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS);
      dbdv[1].col(d)[0]=-dbdv[1].col(d).template segment<2>(1).sum();
    }
    //v[2]
    for(int d=0; d<3; d++) {
      dLHS.setZero();
      dLHS(d,1)=1;
      dbdv[2].col(d).template segment<2>(1)=-invLTL*(dLHS.transpose()*LHS+LHS.transpose()*dLHS)*alpha+invLTL*(dLHS.transpose()*RHS);
      dbdv[2].col(d)[0]=-dbdv[2].col(d).template segment<2>(1).sum();
    }
    //pt
    for(int d=0; d<3; d++) {
      dRHS.setZero();
      dRHS[d]=1;
      dbdpt.col(d).template segment<2>(1)=invLTL*(LHS.transpose()*dRHS);
      dbdpt.col(d)[0]=-dbdpt.col(d).template segment<2>(1).sum();
    }
  }
}
extern void debugBaryDerivativeAll();
}
#endif
