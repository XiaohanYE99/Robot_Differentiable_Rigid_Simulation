#ifndef DISTANCE_FUNCTION_H
#define DISTANCE_FUNCTION_H

#include <tuple>
#include <functional>
#include <Utils/Interp.h>
#include <Utils/SparseUtils.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
template <typename T>
bool triangleContain(const Eigen::Matrix<T,3,1>& a,
                     const Eigen::Matrix<T,3,1> b[3],
                     Eigen::Matrix<T,3,1>* bary3) {
  /*
  *Calculate the barycentric weights for point "a" using the three vertices of triangle "b"
  *If all weights sum to one and locate between [0,1] then the point is inside the triangle.
  *For details see https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Conversion_between_barycentric_and_trilinear_coordinates
  */
  //a=(b[1]-b[0])*alpha+(b[2]-b[0])*beta+b[0]
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3X2T LHS;
  LHS.col(0)=(b[1]-b[0]);
  LHS.col(1)=(b[2]-b[0]);
  Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(a-b[0]));
  if(bary3) {
    bary3->coeffRef(0)=1-bary.sum();
    bary3->coeffRef(1)=bary[0];
    bary3->coeffRef(2)=bary[1];
  }
  return bary.array().isFinite().all() && (bary.array()>=0).all() && (bary.array()<=1).all() && bary.sum()<=1;
}
template <typename T>
bool edgeEdgeIntersect(const Eigen::Matrix<T,3,1> a[2],const Eigen::Matrix<T,3,1> b[2]) {
  /*
  *Solve this: a[0]+(a[1]-a[0])*alpha=b[0]+(b[1]-b[0])*beta
  *If there's a finite solution and is between [0,1],then they are intersected
  */
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3X2T LHS;
  LHS.col(0)=-(a[1]-a[0]);
  LHS.col(1)=(b[1]-b[0]);
  Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(a[0]-b[0]));
  return bary.array().isFinite().all() && (bary.array()>=0).all() && (bary.array()<=1).all();
}
template <typename T>
bool edgeTriangleIntersect(const Eigen::Matrix<T,3,1> a[2],const Eigen::Matrix<T,3,1> b[3]) {
  /*
  *Solve this: a[0]+(a[1]-a[0])*alpha=(b[1]-b[0])*beta+(b[2]-b[0])*gamma+b[0]
  *By solving the parametric line equation and the triangle plane,if a solution is found,then there's intersection
  *If there's no solution,then the coplanar scenarios need to be resolved.
  */
  //a[0]+(a[1]-a[0])*alpha=(b[1]-b[0])*beta+(b[2]-b[0])*gamma+b[0]
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3T LHS;
  LHS.col(0)=-(a[1]-a[0]);
  LHS.col(1)=(b[1]-b[0]);
  LHS.col(2)=(b[2]-b[0]);
  Vec3T RHS=a[0]-b[0];
  Vec3T bary=LHS.inverse()*RHS;
  if(bary.template cast<double>().array().isFinite().all())
    return (bary.array()>=0).all() && (bary.array()<=1).all() && bary[1]+bary[2]<=1;
  /*else {
    //coplanar cases
    if(triangleContain<T>(a[0],b,NULL))
      return true;
    Vec3T e0[2]= {b[0],b[1]};
    if(edgeEdgeIntersect(a,e0))
      return true;
    Vec3T e1[2]= {b[0],b[2]};
    if(edgeEdgeIntersect(a,e1))
      return true;
    Vec3T e2[2]= {b[1],b[1]};
    if(edgeEdgeIntersect(a,e2))
      return true;
  }*/
  return false;
}
template <typename T>
bool triangleTriangleIntersect(const Eigen::Matrix<T,3,1> a[2],const Eigen::Matrix<T,3,1> b[3]) {
  for(int d=0; d<3; d++) {
    const Eigen::Matrix<T,3,1> ea[2]= {a[d],a[(d+1)%3]};
    if(edgeTriangleIntersect(ea,b))
      return true;
    const Eigen::Matrix<T,3,1> eb[2]= {b[d],b[(d+1)%3]};
    if(edgeTriangleIntersect(eb,a))
      return true;
  }
  return false;
}
template <typename T>
T distToSqrLineSegment(const Eigen::Matrix<T,3,1>& pt,
                       const Eigen::Matrix<T,3,1> v[2],
                       Eigen::Map<Eigen::Matrix<T,2,1>> bary,
                       Eigen::Matrix<T,3,1>& cp,
                       char* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Vec3T LHSE=v[1]-v[0],RHSE=pt-v[0];
  T alpha;
  bool systemInvertible=true;
  try {
    alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
  } catch (...) {
    systemInvertible=false;
  }
  if(!isfinite((double)alpha))
    systemInvertible=false;
  if(!systemInvertible) {
    alpha=0;
    if(feat)
      *feat=0;
  } else if(alpha>=0 && alpha<=1) {
    if(feat)
      *feat=-1;
  } else if(alpha<0) {
    alpha=0;
    if(feat)
      *feat=0;
  } else {
    ASSERT(alpha>1)
    alpha=1;
    if(feat)
      *feat=1;
  }
  bary[1]=alpha;
  bary[0]=1-alpha;
  cp=bary[0]*v[0]+bary[1]*v[1];
  return (pt-cp).squaredNorm();
}
template <typename T>
T distToSqrLineSegment(const Eigen::Matrix<T,3,1> a[2],
                       const Eigen::Matrix<T,3,1> b[2],
                       Eigen::Map<Eigen::Matrix<T,2,1>> bary,
                       Eigen::Matrix<T,3,1>& cpa,
                       Eigen::Matrix<T,3,1>& cpb,
                       Eigen::Matrix<char,2,1>* feat) {
  //a[0]*(1-alpha)+a[1]*alpha=a[0]+alpha*(a[1]-a[0])
  //b[0]*(1-beta)+b[1]*beta  =b[0]+beta*(b[1]-b[0])
  //a[0]-b[0]+alpha*(a[1]-a[0])-beta*(b[1]-b[0])=0
  //a[0]*(1-alpha)+a[1]*alpha=a[0]+alpha*(a[1]-a[0])
  //b[0]*(1-beta)+b[1]*beta  =b[0]+beta*(b[1]-b[0])
  //a[0]-b[0]+alpha*(a[1]-a[0])-beta*(b[1]-b[0])=0
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Mat3X2T LHS;
  LHS.col(0)=(a[1]-a[0]);
  LHS.col(1)=-(b[1]-b[0]);
  Vec3T RHS=b[0]-a[0];
  Mat2T LTL=LHS.transpose()*LHS;
  bool systemInvertible=true;
  try {
    bary=LTL.inverse()*(LHS.transpose()*RHS);
  } catch (...) {
    systemInvertible=false;
  }
  if(!bary.array().isFinite().all())
    systemInvertible=false;
  if(!systemInvertible || bary[0]<0 || bary[0] > 1 || bary[1]<0 || bary[1] > 1) {
    T dist,minDist=std::numeric_limits<double>::max();
    T alpha;
    Vec3T LHSE,RHSE;
    bool needTestVA[2]= {true,true};
    bool needTestVB[2]= {true,true};
    for(int d=0; d<2; d++) {
      //b-vertex to a segment
      LHSE=a[1]-a[0],RHSE=b[d]-a[0];
      alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
      if(isfinite((double)alpha) && alpha>=0 && alpha<=1) {
        needTestVB[d]=false;
        dist=(LHSE*alpha-RHSE).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<char,2,1>(-1,d);
          bary[0]=alpha;
          bary[1]=d;
          minDist=dist;
        }
      }
      //a-vertex to b segment
      LHSE=b[1]-b[0],RHSE=a[d]-b[0];
      alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
      if(isfinite((double)alpha) && alpha>=0 && alpha<=1) {
        needTestVA[d]=false;
        dist=(LHSE*alpha-RHSE).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<char,2,1>(d,-1);
          bary[0]=d;
          bary[1]=alpha;
          minDist=dist;
        }
      }
    }
    //vertex
    for(int d=0; d<2; d++)
      for(int d2=0; d2<2; d2++)
        if(needTestVA[d] && needTestVB[d2]) {
          dist=(a[d]-b[d2]).squaredNorm();
          if(dist<minDist) {
            if(feat)
              *feat=Eigen::Matrix<char,2,1>(d,d2);
            bary[0]=d;
            bary[1]=d2;
            minDist=dist;
          }
        }
  } else if(feat)
    *feat=Eigen::Matrix<char,2,1>(-1,-1);
  cpa=interp1D(a[0],a[1],bary[0]);
  cpb=interp1D(b[0],b[1],bary[1]);
  return (cpa-cpb).squaredNorm();
}
template <typename T>
T distToSqrTriangle(const Eigen::Matrix<T,3,1>& pt,
                    const Eigen::Matrix<T,3,1> v[3],
                    Eigen::Map<Eigen::Matrix<T,3,1>> bary,
                    Eigen::Matrix<T,3,1>& cp,
                    Eigen::Matrix<char,2,1>* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Mat3X2T LHS;
  LHS.col(0)=v[1]-v[0];
  LHS.col(1)=v[2]-v[0];
  Vec3T RHS=pt-v[0];
  T alpha;
  //bary
  bool systemInvertible=true;
  try {
    bary.template segment<2>(1)=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*RHS);
    bary[0]=1-bary.template segment<2>(1).sum();
  } catch (...) {
    systemInvertible=false;
  }
  if(!bary.template cast<double>().array().isFinite().all())
    systemInvertible=false;
  if(!systemInvertible || bary.minCoeff()<0) {
    T dist,minDist=std::numeric_limits<double>::max();
    //edge
    bool needTestV[3]= {true,true,true};
    for(int d=0; d<3; d++) {
      //|v[(d+1)%3+1]*alpha+v[d+1]*(1-alpha)-v[0]|^2
      Vec3T LHSE=v[(d+1)%3]-v[d],RHSE=pt-v[d];
      systemInvertible=true;
      try {
        alpha=RHSE.dot(LHSE)/LHSE.dot(LHSE);
      } catch (...) {
        systemInvertible=false;
      }
      if(!isfinite((double)alpha))
        systemInvertible=false;
      if(systemInvertible && alpha>=0 && alpha<=1) {
        needTestV[(d+1)%3]=needTestV[d]=false;
        dist=(LHSE*alpha-RHSE).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<char,2,1>((d+1)%3,d);
          bary.setZero();
          bary[(d+1)%3]=alpha;
          bary[d]=1-alpha;
          minDist=dist;
        }
      }
    }
    //vertex
    for(int d=0; d<3; d++)
      if(needTestV[d]) {
        dist=(v[d]-pt).squaredNorm();
        if(dist<minDist) {
          if(feat)
            *feat=Eigen::Matrix<char,2,1>(d,-1);
          bary.setUnit(d);
          minDist=dist;
        }
      }
  } else if(feat)
    *feat=Eigen::Matrix<char,2,1>(-1,-1);
  cp=bary[0]*v[0]+bary[1]*v[1]+bary[2]*v[2];
  return (pt-cp).squaredNorm();
}
template <typename T>
T distToSqrTetrahedron(const Eigen::Matrix<T,3,1>& pt,
                       const Eigen::Matrix<T,3,1> v[4],
                       Eigen::Map<Eigen::Matrix<T,4,1>> bary,
                       Eigen::Matrix<T,3,1>& cp,
                       Eigen::Matrix<char,3,1>* feat) {
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Mat3T LHS;
  LHS.col(0)=v[1]-v[0];
  LHS.col(1)=v[2]-v[0];
  LHS.col(2)=v[3]-v[0];
  Vec3T RHS=pt-v[0];
  bool inside=true;
  try {
    bary.template segment<3>(1)=LHS.inverse()*RHS;
    bary[0]=1-bary.template segment<3>(1).sum();
  } catch (...) {
    inside=false;
  }
  if(!bary.template cast<double>().array().isFinite().all())
    inside=false;
  if(!inside || bary.minCoeff()<0) {
    T dist,minDist=std::numeric_limits<double>::max();
    //consider 4 triangles
    Eigen::Matrix<T,3,1> baryT;
    Eigen::Matrix<T,3,1> vt[3],cpT;
    Eigen::Matrix<char,2,1> featT;
    for(int d=0; d<4; d++) {
      vt[0]=v[(d+1)%4];
      vt[1]=v[(d+2)%4];
      vt[2]=v[(d+3)%4];
      dist=distToSqrTriangle(pt,vt,Eigen::Map<Vec3T>(baryT.data()),cpT,feat?&featT:NULL);
      if(dist<minDist) {
        minDist=dist;
        bary.setConstant(0);
        bary[(d+1)%4]=baryT[0];
        bary[(d+2)%4]=baryT[1];
        bary[(d+3)%4]=baryT[2];
        cp=cpT;
        if(feat) {
          if(featT==Eigen::Matrix<char,2,1>(-1,-1))
            *feat=Eigen::Matrix<char,3,1>((d+1)%4,(d+2)%4,(d+3)%4);
          else if(featT[1]==-1)
            *feat=Eigen::Matrix<char,3,1>((d+1+featT[0])%4,-1,-1);
          else
            *feat=Eigen::Matrix<char,3,1>((d+1+featT[0])%4,(d+1+featT[1])%4,-1);
        }
      }
    }
  } else if(feat)
    *feat=Eigen::Matrix<char,3,1>(-1,-1,-1);
  cp=bary[0]*v[0]+bary[1]*v[1]+bary[2]*v[2]+bary[3]*v[3];
  return (pt-cp).squaredNorm();
}
template <typename T>
void debugDistToSqrTetrahedron();
}
#endif
