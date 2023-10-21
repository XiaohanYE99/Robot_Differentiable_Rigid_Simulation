#ifndef SPATIAL_ROTATION_UTILS_H
#define SPATIAL_ROTATION_UTILS_H

#include "RotationUtils.h"
#include "DebugGradient.h"
#include <random>

namespace PHYSICSMOTION {
//spatial rotation: expWZ
template <typename T,typename TV,typename TX,typename TS,typename TD>   //rotation 3D euler angle
typename Eigen::Matrix<T,3,3> spatialExpWZ(T q,const T* dq,const T* ddq,TX* X=NULL,TS* S=NULL,TV* vJ=NULL,TV* dvJ=NULL,TD* DvJDq=NULL,TD* DdvJDq=NULL,TD* DdvJDdq=NULL) {
  //this is a special type of Euler-Angle in favor of our definition in RootInvariantMLP
  T sz,cz;
  SinCosTraits<T>::sincos(q,&sz,&cz);

  typename Eigen::Matrix<T,3,3> R;
  R(0,0)=cz;
  R(0,1)=-sz;
  R(0,2)=0;
  R(1,0)=sz;
  R(1,1)=cz;
  R(1,2)=0;
  R(2,0)=0;
  R(2,1)=0;
  R(2,2)=1;
  if(X) {
    X->setZero();
    X->template block<3,3>(0,0)=X->template block<3,3>(3,3)=R;
  }

  if(S) {
    S->setZero();
    (*S)(2,0)=1;
  }

  if(vJ) {
    vJ->setZero();
    (*vJ)[2]=dq[0];
  }

  if(dvJ) {
    dvJ->setZero();
    if(ddq)
      (*dvJ)[2]=*ddq;
  }

  if(DvJDq)
    DvJDq->setZero();

  if(DdvJDq)
    DdvJDq->setZero();

  if(DdvJDdq)
    DdvJDdq->setZero();
  return R;
}
//spatial rotation: expWYZ
template <typename T,typename TCV,typename TV,typename TX,typename TS,typename TD>   //rotation 3D euler angle
typename Eigen::Matrix<T,3,3> spatialExpWYZ(const TCV& q,const TCV* dq,const TCV* ddq,TX* X=NULL,TS* S=NULL,TV* vJ=NULL,TV* dvJ=NULL,TD* DvJDq=NULL,TD* DdvJDq=NULL,TD* DdvJDdq=NULL) {
  //this is a special type of Euler-Angle in favor of our definition in RootInvariantMLP
  T sy,cy,sz,cz;
  SinCosTraits<T>::sincos(q[0],&sy,&cy);
  SinCosTraits<T>::sincos(q[1],&sz,&cz);

  typename Eigen::Matrix<T,3,3> R;
  R(0,0)=cy*cz;
  R(0,1)=-sz;
  R(0,2)=sy*cz;
  R(1,0)=cy*sz;
  R(1,1)=cz;
  R(1,2)=sy*sz;
  R(2,0)=-sy;
  R(2,1)=0;
  R(2,2)=cy;
  if(X) {
    X->setZero();
    X->template block<3,3>(0,0)=X->template block<3,3>(3,3)=R;
  }

  if(S) {
    S->setZero();
    (*S)(1,0)=1;
    (*S)(0,1)=-sy;
    (*S)(2,1)=cy;
  }

  if(vJ) {
    vJ->setZero();
    (*vJ)[0]-=(*dq)[1]*sy;
    (*vJ)[1]=(*dq)[0];
    (*vJ)[2]=(*dq)[1]*cy;
  }

  if(dvJ) {
    dvJ->setZero();
    if(dq) {
      (*dvJ)[0]-=(*dq)[0]*(*dq)[1]*cy;
      (*dvJ)[2]-=(*dq)[0]*(*dq)[1]*sy;
    }
    if(ddq) {
      (*dvJ)[0]-=(*ddq)[1]*sy;
      (*dvJ)[1]+=(*ddq)[0];
      (*dvJ)[2]+=(*ddq)[1]*cy;
    }
  }

  if(DvJDq) {
    DvJDq->setZero();
    (*DvJDq)(0,1)-=(*dq)[0]*cy;
    (*DvJDq)(2,1)-=(*dq)[0]*sy;
  }

  if(DdvJDq) {
    DdvJDq->setZero();
    if(dq)
      (*DdvJDq)(1,1)-=(*dq)[0]*(*dq)[1];
    if(ddq) {
      (*DdvJDq)(0,1)-=(*ddq)[0]*cy;
      (*DdvJDq)(2,1)-=(*ddq)[0]*sy;
    }
  }

  if(DdvJDdq) {
    DdvJDdq->setZero();
    (*DdvJDdq)(0,0)-=(*dq)[1]*cy;
    (*DdvJDdq)(0,1)-=(*dq)[0]*cy;
    (*DdvJDdq)(2,0)-=(*dq)[1]*sy;
    (*DdvJDdq)(2,1)-=(*dq)[0]*sy;
  }
  return R;
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialExpWYZDSTDqf(const TCV& q,TTAU& DSTDqTf,const TF& f) {
  T sy,cy;
  SinCosTraits<T>::sincos(q[0],&sy,&cy);
  DSTDqTf(1,0)-=f[2]*sy+f[0]*cy;
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialExpWYZDSDqf(const TCV& q,TTAU& DSTDqf,const TF& f) {
  T sy,cy;
  SinCosTraits<T>::sincos(q[0],&sy,&cy);
  DSTDqf(0,0)-=f[1]*cy;
  DSTDqf(2,0)-=f[1]*sy;
}
//spatial rotation: eulerXYZ
template <typename T,typename TCV,typename TV,typename TX,typename TS,typename TD>   //rotation 3D euler angle
typename Eigen::Matrix<T,3,3> spatialEulerX1Y3Z2(const TCV& q,const TCV* dq,const TCV* ddq,TX* X=NULL,TS* S=NULL,TV* vJ=NULL,TV* dvJ=NULL,TD* DvJDq=NULL,TD* DdvJDq=NULL,TD* DdvJDdq=NULL) {
  //this is a special type of Euler-Angle in favor of our definition in RootInvariantMLP
  T sx,cx,sy,cy,sz,cz;
  SinCosTraits<T>::sincos(q[0],&sx,&cx);
  SinCosTraits<T>::sincos(q[1],&sy,&cy);
  SinCosTraits<T>::sincos(q[2],&sz,&cz);

  typename Eigen::Matrix<T,3,3> R;
  R(0,0)=cy*cz;
  R(0,1)=sx*sy-cx*cy*sz;
  R(0,2)=sx*cy*sz+cx*sy;
  R(1,0)=sz;
  R(1,1)=cx*cz;
  R(1,2)=-sx*cz;
  R(2,0)=-sy*cz;
  R(2,1)=cx*sy*sz+sx*cy;
  R(2,2)=cx*cy-sx*sy*sz;
  if(X) {
    X->setZero();
    X->template block<3,3>(0,0)=X->template block<3,3>(3,3)=R;
  }

  if(S) {
    S->setZero();
    (*S)(0,0)=1;
    (*S)(0,1)=sz;
    (*S)(1,1)=cx*cz;
    (*S)(2,1)=-sx*cz;
    (*S)(1,2)=sx;
    (*S)(2,2)=cx;
  }

  if(vJ) {
    vJ->setZero();
    (*vJ)[0]=(*dq)[1]*sz+(*dq)[0];
    (*vJ)[1]=(*dq)[1]*cx*cz+(*dq)[2]*sx;
    (*vJ)[2]=(*dq)[2]*cx-(*dq)[1]*sx*cz;
  }

  if(dvJ) {
    dvJ->setZero();
    if(dq) {
      (*dvJ)[0]=(*dq)[1]*(*dq)[2]*cz;
      (*dvJ)[1]=-(*dq)[1]*(*dq)[2]*cx*sz-(*dq)[0]*(*dq)[1]*sx*cz+(*dq)[0]*(*dq)[2]*cx;
      (*dvJ)[2]=(*dq)[1]*(*dq)[2]*sx*sz-(*dq)[0]*(*dq)[1]*cx*cz-(*dq)[0]*(*dq)[2]*sx;
    }
    if(ddq) {
      (*dvJ)[0]+=(*ddq)[1]*sz+(*ddq)[0];
      (*dvJ)[1]+=(*ddq)[1]*cx*cz+(*ddq)[2]*sx;
      (*dvJ)[2]+=-(*ddq)[1]*sx*cz+(*ddq)[2]*cx;
    }
  }

  if(DvJDq) {
    DvJDq->setZero();
    (*DvJDq)(0,1)=(*dq)[2]*cz;
    (*DvJDq)(1,1)=-(*dq)[2]*cx*sz-(*dq)[0]*sx*cz;
    (*DvJDq)(2,1)=(*dq)[2]*sx*sz-(*dq)[0]*cx*cz;
    (*DvJDq)(1,2)=(*dq)[0]*cx;
    (*DvJDq)(2,2)=-(*dq)[0]*sx;
  }

  if(DdvJDq) {
    DdvJDq->setZero();
    if(dq) {
      (*DdvJDq)(0,1)=-(*dq)[0]*(*dq)[1]*cz*cz;
      (*DdvJDq)(0,2)=-(*dq)[0]*(*dq)[2];
      (*DdvJDq)(1,1)=((*dq)[0]*(*dq)[2]*sx+(*dq)[0]*(*dq)[1]*cx*cz)*sz-(*dq)[1]*(*dq)[2]*sx;
      (*DdvJDq)(1,2)=(*dq)[0]*(*dq)[1]*sx*sz;
      (*DdvJDq)(2,1)=((*dq)[0]*(*dq)[2]*cx-(*dq)[0]*(*dq)[1]*sx*cz)*sz-(*dq)[1]*(*dq)[2]*cx;
      (*DdvJDq)(2,2)=(*dq)[0]*(*dq)[1]*cx*sz;
    }
    if(ddq) {
      (*DdvJDq)(0,1)+=(*ddq)[2]*cz;
      (*DdvJDq)(1,1)+=-(*ddq)[2]*cx*sz-(*ddq)[0]*sx*cz;
      (*DdvJDq)(1,2)+=(*ddq)[0]*cx;
      (*DdvJDq)(2,1)+=(*ddq)[2]*sx*sz-(*ddq)[0]*cx*cz;
      (*DdvJDq)(2,2)+=-(*ddq)[0]*sx;
    }
  }

  if(DdvJDdq) {
    DdvJDdq->setZero();
    (*DdvJDdq)(0,1)=(*dq)[2]*cz;
    (*DdvJDdq)(0,2)=(*dq)[1]*cz;
    (*DdvJDdq)(1,0)=(*dq)[2]*cx-(*dq)[1]*sx*cz;
    (*DdvJDdq)(1,1)=-(*dq)[2]*cx*sz-(*dq)[0]*sx*cz;
    (*DdvJDdq)(1,2)=(*dq)[0]*cx-(*dq)[1]*cx*sz;
    (*DdvJDdq)(2,0)=-(*dq)[1]*cx*cz-(*dq)[2]*sx;
    (*DdvJDdq)(2,1)=(*dq)[2]*sx*sz-(*dq)[0]*cx*cz;
    (*DdvJDdq)(2,2)=(*dq)[1]*sx*sz-(*dq)[0]*sx;
  }
  return R;
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialEulerX1Y3Z2DSTDqf(const TCV& q,TTAU& DSTDqTf,const TF& f) {
  T sx,cx,sz,cz;
  SinCosTraits<T>::sincos(q[0],&sx,&cx);
  SinCosTraits<T>::sincos(q[2],&sz,&cz);
  DSTDqTf(1,0)-=f[1]*sx*cz+f[2]*cx*cz;
  DSTDqTf(2,0)+=f[1]*cx-f[2]*sx;
  DSTDqTf(1,2)+=f[2]*sx*sz-f[1]*cx*sz+f[0]*cz;
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialEulerX1Y3Z2DSDqf(const TCV& q,TTAU& DSTDqf,const TF& f) {
  T sx,cx,sz,cz;
  SinCosTraits<T>::sincos(q[0],&sx,&cx);
  SinCosTraits<T>::sincos(q[2],&sz,&cz);
  DSTDqf(1,0)+=f[2]*cx-f[1]*sx*cz;
  DSTDqf(2,0)-=f[1]*cx*cz+f[2]*sx;
  DSTDqf(0,2)+=f[1]*cz;
  DSTDqf(1,2)-=f[1]*cx*sz;
  DSTDqf(2,2)+=f[1]*sx*sz;
}
//spatial rotation: expW
template <typename T,typename TCV,typename TV,typename TX,typename TS,typename TD>   //rotation 3D euler angle
typename Eigen::Matrix<T,3,3> spatialExpW(const TCV& q,const TCV* dq,const TCV* ddq,TX* X=NULL,TS* S=NULL,TV* vJ=NULL,TV* dvJ=NULL,TD* DvJDq=NULL,TD* DdvJDq=NULL,TD* DdvJDdq=NULL) {
  typedef typename Eigen::Matrix<T,3,1> Vec3T;
  typedef typename Eigen::Matrix<T,3,3> Mat3T;
  Vec3T diffV[3],ddiffV[9],ddiffVLambda[9];
  Mat3T R=expWGradV<T,Vec3T>(q,diffV,ddiffV);
  if(X) {
    X->setZero();
    X->template block<3,3>(0,0)=X->template block<3,3>(3,3)=R;
  }

  if(S) {
    S->setZero();
    S->template block<3,1>(0,0)=R.transpose()*diffV[0];
    S->template block<3,1>(0,1)=R.transpose()*diffV[1];
    S->template block<3,1>(0,2)=R.transpose()*diffV[2];
  }

  if(vJ) {
    vJ->setZero();
    vJ->template segment<3>(0)=R.transpose()*(diffV[0]*(*dq)[0]+diffV[1]*(*dq)[1]+diffV[2]*(*dq)[2]);
  }

  if(dvJ) {
    dvJ->setZero();
    if(dq) {
      dvJ->template segment<3>(0)=R.transpose()*
                                  (ddiffV[0]*(*dq)[0]*(*dq)[0]+ddiffV[1]*(*dq)[1]*(*dq)[0]+ddiffV[2]*(*dq)[2]*(*dq)[0]+
                                   ddiffV[3]*(*dq)[0]*(*dq)[1]+ddiffV[4]*(*dq)[1]*(*dq)[1]+ddiffV[5]*(*dq)[2]*(*dq)[1]+
                                   ddiffV[6]*(*dq)[0]*(*dq)[2]+ddiffV[7]*(*dq)[1]*(*dq)[2]+ddiffV[8]*(*dq)[2]*(*dq)[2]);
    }
    if(ddq)
      dvJ->template segment<3>(0)+=R.transpose()*(diffV[0]*(*ddq)[0]+diffV[1]*(*ddq)[1]+diffV[2]*(*ddq)[2]);
  }

  if(DvJDq) {
    DvJDq->setZero();
    DvJDq->template block<3,1>(0,0)=R.transpose()*(ddiffV[0]*(*dq)[0]+ddiffV[3]*(*dq)[1]+ddiffV[6]*(*dq)[2]);
    DvJDq->template block<3,1>(0,1)=R.transpose()*(ddiffV[1]*(*dq)[0]+ddiffV[4]*(*dq)[1]+ddiffV[7]*(*dq)[2]);
    DvJDq->template block<3,1>(0,2)=R.transpose()*(ddiffV[2]*(*dq)[0]+ddiffV[5]*(*dq)[1]+ddiffV[8]*(*dq)[2]);
  }

  if(DdvJDq) {
    DdvJDq->setZero();
    if(dq) {
      expWGradVLambda<T,Vec3T>(q,*dq,ddiffVLambda);
      DdvJDq->template block<3,1>(0,0)=R.transpose()*(ddiffVLambda[0]*(*dq)[0]+ddiffVLambda[3]*(*dq)[1]+ddiffVLambda[6]*(*dq)[2]);
      DdvJDq->template block<3,1>(0,1)=R.transpose()*(ddiffVLambda[1]*(*dq)[0]+ddiffVLambda[4]*(*dq)[1]+ddiffVLambda[7]*(*dq)[2]);
      DdvJDq->template block<3,1>(0,2)=R.transpose()*(ddiffVLambda[2]*(*dq)[0]+ddiffVLambda[5]*(*dq)[1]+ddiffVLambda[8]*(*dq)[2]);
    }
    if(ddq) {
      DdvJDq->template block<3,1>(0,0)+=R.transpose()*(ddiffV[0]*(*ddq)[0]+ddiffV[3]*(*ddq)[1]+ddiffV[6]*(*ddq)[2]);
      DdvJDq->template block<3,1>(0,1)+=R.transpose()*(ddiffV[1]*(*ddq)[0]+ddiffV[4]*(*ddq)[1]+ddiffV[7]*(*ddq)[2]);
      DdvJDq->template block<3,1>(0,2)+=R.transpose()*(ddiffV[2]*(*ddq)[0]+ddiffV[5]*(*ddq)[1]+ddiffV[8]*(*ddq)[2]);
    }
  }

  if(DdvJDdq) {
    DdvJDdq->setZero();
    DdvJDdq->template block<3,1>(0,0)=R.transpose()*(ddiffV[0]*(*dq)[0]+ddiffV[1]*(*dq)[1]+ddiffV[2]*(*dq)[2]+ddiffV[0]*(*dq)[0]+ddiffV[3]*(*dq)[1]+ddiffV[6]*(*dq)[2]);
    DdvJDdq->template block<3,1>(0,1)=R.transpose()*(ddiffV[3]*(*dq)[0]+ddiffV[4]*(*dq)[1]+ddiffV[5]*(*dq)[2]+ddiffV[1]*(*dq)[0]+ddiffV[4]*(*dq)[1]+ddiffV[7]*(*dq)[2]);
    DdvJDdq->template block<3,1>(0,2)=R.transpose()*(ddiffV[6]*(*dq)[0]+ddiffV[7]*(*dq)[1]+ddiffV[8]*(*dq)[2]+ddiffV[2]*(*dq)[0]+ddiffV[5]*(*dq)[1]+ddiffV[8]*(*dq)[2]);
  }
  return R;
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialExpWDSTDqf(const TCV& q,TTAU& DSTDqTf,const TF& f) {
  typedef typename Eigen::Matrix<T,3,1> Vec3T;
  typedef typename Eigen::Matrix<T,3,3> Mat3T;
  Vec3T diffV[3],ddiffV[9];
  Mat3T R=expWGradV<T,Vec3T>(q,diffV,ddiffV);
  Vec3T Rf=R*f.template segment<3>(0);
  DSTDqTf(0,0)+=ddiffV[0].dot(Rf)-diffV[0].cross(diffV[0]).dot(Rf);
  DSTDqTf(1,0)+=ddiffV[3].dot(Rf)-diffV[0].cross(diffV[1]).dot(Rf);
  DSTDqTf(2,0)+=ddiffV[6].dot(Rf)-diffV[0].cross(diffV[2]).dot(Rf);

  DSTDqTf(0,1)+=ddiffV[1].dot(Rf)-diffV[1].cross(diffV[0]).dot(Rf);
  DSTDqTf(1,1)+=ddiffV[4].dot(Rf)-diffV[1].cross(diffV[1]).dot(Rf);
  DSTDqTf(2,1)+=ddiffV[7].dot(Rf)-diffV[1].cross(diffV[2]).dot(Rf);

  DSTDqTf(0,2)+=ddiffV[2].dot(Rf)-diffV[2].cross(diffV[0]).dot(Rf);
  DSTDqTf(1,2)+=ddiffV[5].dot(Rf)-diffV[2].cross(diffV[1]).dot(Rf);
  DSTDqTf(2,2)+=ddiffV[8].dot(Rf)-diffV[2].cross(diffV[2]).dot(Rf);
}
template <typename T,typename TCV,typename TTAU,typename TF>
void spatialExpWDSDqf(const TCV& q,TTAU& DSTDqf,const TF& f) {
  typedef typename Eigen::Matrix<T,3,1> Vec3T;
  typedef typename Eigen::Matrix<T,3,3> Mat3T;
  Vec3T diffV[3],ddiffV[9];
  Mat3T R=expWGradV<T,Vec3T>(q,diffV,ddiffV);
  Vec3T wf=diffV[0]*f[0]+diffV[1]*f[1]+diffV[2]*f[2];
  DSTDqf.template block<3,1>(0,0)+=R.transpose()*(ddiffV[0]*f[0]+ddiffV[3]*f[1]+ddiffV[6]*f[2]-diffV[0].cross(wf));
  DSTDqf.template block<3,1>(0,1)+=R.transpose()*(ddiffV[1]*f[0]+ddiffV[4]*f[1]+ddiffV[7]*f[2]-diffV[1].cross(wf));
  DSTDqf.template block<3,1>(0,2)+=R.transpose()*(ddiffV[2]*f[0]+ddiffV[5]*f[1]+ddiffV[8]*f[2]-diffV[2].cross(wf));
}
}

#endif
