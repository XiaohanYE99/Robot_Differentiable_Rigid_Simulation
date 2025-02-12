#include "DistanceEnergy.h"
#include <Environment/DistanceFunction.h>
#include <Utils/DebugGradient.h>
#include <random>

namespace PHYSICSMOTION {
#define MOVE_TO_CENTER
//VTDistanceEnergy
template <typename T>
VTDistanceEnergy<T>::VTDistanceEnergy(const Vec3T triA,const Vec3T& triB,const Vec3T& triC,const Vec3T& point):_point(point) {
  _tri[0]=triA;
  _tri[1]=triB;
  _tri[2]=triC;
}
template <typename T>
VTDistanceEnergy<T>::VTDistanceEnergy(const Vec3T tri[3],const Vec3T& point):_point(point) {
  _tri[0]=tri[0];
  _tri[1]=tri[1];
  _tri[2]=tri[2];
}
template <typename T>
bool VTDistanceEnergy<T>::eval(T* E,Vec12T* G,Mat12T* H) {
  Vec3T v[4],bary,cp;
  v[0]=_point;
  v[1]=_tri[0];
  v[2]=_tri[1];
  v[3]=_tri[2];
#ifdef MOVE_TO_CENTER
  cp=(v[0]+v[1]+v[2]+v[3])/4;
  for(int j=0; j<4; j++)
    v[j]-=cp;
#endif
  T dSqr=distToSqrTriangle(v[0],v+1,Eigen::Map<Vec3T>(bary.data()),cp,&_feat),d=sqrt(dSqr);
  if(E)
    (*E)=d;
  if(_feat[0]==-1) {
    //closest feature is face
    Vec3T a=v[0]-v[1],b=v[2]-v[1],c=v[3]-v[1],n=(b).cross(c);
    char nSgn=(v[0]-cp).dot(n)<0?-1:1;
    d*=nSgn;
    T invNLen=1/n.norm();
    Vec3T nn=n*invNLen,bnn,cnn;
    Mat3T projN;
    Vec9T g;
    Vec12T GSum;
    Mat9T h;
    Mat12T HSum;
    if(G || H) {
      g.template segment<3>(0)=nn;
      bnn=b.cross(nn),cnn=c.cross(nn);
      g.template segment<3>(3)=(c.cross(a)-d*cnn)*invNLen;
      g.template segment<3>(6)=(a.cross(b)+d*bnn)*invNLen;
    }
    if(H) {
      h.template block<3,3>(0,0).setZero();
      projN=Mat3T::Identity()-nn*nn.transpose();
      //ga,b
      h.template block<3,3>(0,3)=-projN*cross<T>(c)*invNLen;
      h.template block<3,3>(3,0)=h.template block<3,3>(0,3).transpose();
      //ga,c
      h.template block<3,3>(0,6)=projN*cross<T>(b)*invNLen;
      h.template block<3,3>(6,0)=h.template block<3,3>(0,6).transpose();
      //gb,b
      h.template block<3,3>(3,3)=-(cnn*g.template segment<3>(3).transpose()+g.template segment<3>(3)*cnn.transpose())*invNLen;
      h.template block<3,3>(3,3)-=cross<T>(c)*h.template block<3,3>(0,3)*d*invNLen;
      //gc,c
      h.template block<3,3>(6,6)=(bnn*g.template segment<3>(6).transpose()+g.template segment<3>(6)*bnn.transpose())*invNLen;
      h.template block<3,3>(6,6)+=cross<T>(b)*h.template block<3,3>(0,6)*d*invNLen;
      //gb,c
      h.template block<3,3>(3,6)=(-cross<T>(a)-d*cross<T>(c)*h.template block<3,3>(0,6)+d*cross<T>(nn)-cnn*g.template segment<3>(6).transpose()+g.template segment<3>(3)*bnn.transpose())*invNLen;
      h.template block<3,3>(6,3)=h.template block<3,3>(3,6).transpose();
    }
    if(G)
      *G=getIF()*g*nSgn;
    if(H)
      *H=(getIF()*h*getIF().transpose())*nSgn;
  } else if(_feat[1]==-1) {
    Vec3T g;
    Mat3T h;
    if(G || H)
      g=(v[0]-cp)/d;
    if(H)
      h=(Mat3T::Identity()-g*g.transpose())/d;
    if(G) {
      G->setZero();
      G->template segment<3>(0)=g;
      G->template segment<3>((1+_feat[0])*3)=-g;
    }
    if(H) {
      H->setZero();
      H->template block<3,3>(0,0)=h;
      H->template block<3,3>(0,(1+_feat[0])*3)=-h;
      H->template block<3,3>((1+_feat[0])*3,0)=-h;
      H->template block<3,3>((1+_feat[0])*3,(1+_feat[0])*3)=h;
    }
  } else {
    std::array<int,3> id= {0,_feat[0]+1,_feat[1]+1};
    //closest _feature is edge
    Vec3T vE[3]= {v[0],v[_feat[0]+1],v[_feat[1]+1]};
    Vec3T a=vE[0]-vE[1],b=vE[0]-vE[2],c=vE[1]-vE[2];
    Vec3T n=a.cross(b);
    T nLen=n.norm(),invNLen=1/nLen,invCLen=1/c.norm();
    Vec3T nn=n*invNLen,cn=c*invCLen;
    Mat3T projN;
    Vec9T g,GSum;
    Mat9T h,HSum;
    if(G || H) {
      g.template segment<3>(0)=b.cross(nn)*invCLen;
      g.template segment<3>(3)=-a.cross(nn)*invCLen;
      g.template segment<3>(6)=-cn*nLen*invCLen*invCLen;
    }
    if(H) {
      projN=Mat3T::Identity()-nn*nn.transpose();
      //ga,a
      h.template block<3,3>(0,0)=-cross<T>(b)*projN*cross<T>(b)*invNLen*invCLen;
      //ga,b
      h.template block<3,3>(0,3)=cross<T>(b)*projN*cross<T>(a)*invNLen*invCLen-cross<T>(nn)*invCLen;
      h.template block<3,3>(3,0)=h.template block<3,3>(0,3).transpose();
      //ga,c
      h.template block<3,3>(0,6)=-g.template segment<3>(0)*cn.transpose()*invCLen;
      h.template block<3,3>(6,0)=h.template block<3,3>(0,6).transpose();
      //gb,b
      h.template block<3,3>(3,3)=-cross<T>(a)*projN*cross<T>(a)*invNLen*invCLen;
      //gb,c
      h.template block<3,3>(3,6)=-g.template segment<3>(3)*cn.transpose()*invCLen;
      h.template block<3,3>(6,3)=h.template block<3,3>(3,6).transpose();
      //gc,c
      h.template block<3,3>(6,6)=(3*cn*cn.transpose()-Mat3T::Identity())*nLen*invCLen*invCLen*invCLen;
    }
    if(G) {
      G->setZero();
      GSum=getIE()*g;
      for(int i=0; i<(int)id.size(); ++i)
        G->template segment<3>(id[i]*3)=GSum.template segment<3>(i*3);
    }
    if(H) {
      H->setZero();
      HSum=getIE()*h*getIE().transpose();
      for(int i=0; i<(int) id.size(); ++i)
        for(int j=0; j<(int) id.size(); ++j)
          H->template block<3,3>(id[i]*3,id[j]*3)=HSum.template block<3,3>(i*3,j*3);
    }
  }
  return true;
}
template <typename T>
const typename VTDistanceEnergy<T>::SMatT& VTDistanceEnergy<T>::getIF() {
  static SMatT I;
  OMP_CRITICAL_{
    if(I.size()==0) {
      STrips trips;
      //a
      addBlockId<T>(trips,0,0,3,1);
      addBlockId<T>(trips,3,0,3,-1);
      //b
      addBlockId<T>(trips,6,3,3,1);
      addBlockId<T>(trips,3,3,3,-1);
      //c
      addBlockId<T>(trips,9,6,3,1);
      addBlockId<T>(trips,3,6,3,-1);
      I.resize(12,9);
      I.setFromTriplets(trips.begin(),trips.end());
    }
  }
  return I;
}
template <typename T>
const typename VTDistanceEnergy<T>::SMatT& VTDistanceEnergy<T>::getIE() {
  static SMatT I;
  OMP_CRITICAL_{
    if(I.size()==0) {
      STrips trips;
      //a
      addBlockId<T>(trips,3,0,3,-1);
      addBlockId<T>(trips,0,0,3,1);
      //b
      addBlockId<T>(trips,6,3,3,-1);
      addBlockId<T>(trips,0,3,3,1);
      //c
      addBlockId<T>(trips,3,6,3,1);
      addBlockId<T>(trips,6,6,3,-1);
      I.resize(9,9);
      I.setFromTriplets(trips.begin(),trips.end());
    }
  }
  return I;
}
template <typename T>
bool VTDistanceEnergy<T>::debug(const Vec3T tri[3],const Vec3T& point,T perturbRange,const Eigen::Matrix<char,2,1>& feat) {
  _tri[0]=tri[0];
  _tri[1]=tri[1];
  _tri[2]=tri[2];
  _point=point;
  T E,E2;
  Mat12T H;
  Vec12T G,G2,dx;
  DEFINE_NUMERIC_DELTA_T(T)
  dx=Vec12T::Random();
  _point=_point+Vec::Random(3)*perturbRange;
  for(int i=0; i<3; ++i)
    _tri[i]+=Vec::Random(3)*perturbRange;
  if(!eval(&E,&G,&H) || E==0 || !isfinite(E))
    return false;
  _point+=dx.template segment<3>(0)*DELTA;
  for(int i=0; i<3; ++i)
    _tri[i]+=dx.template segment<3>((i+1)*3)*DELTA;
  if(!eval(&E2,&G2,NULL))
    return false;
  if(_feat!=feat)
    return false;
  std::cout << "Feat=(" << (int)_feat[0] << "," << (int)_feat[1] << ")" << std::endl;
  DEBUG_GRADIENT("dE",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("dG",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
  return true;
}
template <typename T>
void VTDistanceEnergy<T>::debug(T perturbRange,const Eigen::Matrix<char,2,1>& feat) {
  while(true) {
    Vec3T point=Vec3T::Random();
    Vec3T tri[3]= {Vec3T::Random(),Vec3T::Random(),Vec3T::Random()};
    if(debug(tri,point,perturbRange,feat))
      break;
  }
}
//EEDistanceEnergy
template <typename T>
EEDistanceEnergy<T>::EEDistanceEnergy(const Vec3T edge1A,const Vec3T edge1B,const Vec3T edge2A,const Vec3T edge2B,bool JTJApprox):_JTJApprox(JTJApprox) {
  _edge1[0]=edge1A;
  _edge1[1]=edge1B;
  _edge2[0]=edge2A;
  _edge2[1]=edge2B;
}
template <typename T>
EEDistanceEnergy<T>::EEDistanceEnergy(const Vec3T edge1[2],const Vec3T edge2[2],bool JTJApprox):_JTJApprox(JTJApprox) {
  _edge1[0]=edge1[0];
  _edge1[1]=edge1[1];
  _edge2[0]=edge2[0];
  _edge2[1]=edge2[1];
}
template <typename T>
bool EEDistanceEnergy<T>::eval(T* E,Vec12T* G,Mat12T* H) {
  Vec2T bary;
  Vec3T v[4],cpa,cpb;
  v[0]=_edge1[0];
  v[1]=_edge1[1];
  v[2]=_edge2[0];
  v[3]=_edge2[1];
#ifdef MOVE_TO_CENTER
  cpa=(v[0]+v[1]+v[2]+v[3])/4;
  for(int j=0; j<4; j++)
    v[j]-=cpa;
#endif
  T dSqr=distToSqrLineSegment(v,v+2,Eigen::Map<Vec2T>(bary.data()),cpa,cpb,&_feat),d=sqrt(dSqr);
  if(E)
    (*E)=d;
  if(_feat[0]==-1&& _feat[1]==-1) {
    //closest feature is face
    Vec3T a=v[0]-v[2],b=v[0]-v[1],c=v[2]-v[3],n=(b).cross(c);
    char nSgn=a.dot(n)<0?-1:1;
    d*= nSgn;
    T invNLen=1/n.norm();
    Vec3T nn=n*invNLen,bnn,cnn;
    Mat3T projN;
    Vec9T g;
    Vec12T GSum;
    Mat9T h;
    Mat12T HSum;
    if(G || H) {
      g.resize(9);
      g.template segment<3>(0)=nn;
      bnn=b.cross(nn),cnn=c.cross(nn);
      g.template segment<3>(3)=(c.cross(a)-d*cnn)*invNLen;
      g.template segment<3>(6)=(a.cross(b)+d*bnn)*invNLen;
    }
    if(H) {
      h.resize(9,9);
      h.template block<3,3>(0,0).setZero();
      projN=Mat3T::Identity()-nn*nn.transpose();
      //ga,b
      h.template block<3,3>(0,3)=-projN*cross<T>(c)*invNLen;
      h.template block<3,3>(3,0)=h.template block<3,3>(0,3).transpose();
      //ga,c
      h.template block<3,3>(0,6)=projN*cross<T>(b)*invNLen;
      h.template block<3,3>(6,0)=h.template block<3,3>(0,6).transpose();
      //gb,b
      h.template block<3,3>(3,3)=-(cnn*g.template segment<3>(3).transpose()+g.template segment<3>(3)*cnn.transpose())*invNLen;
      h.template block<3,3>(3,3)-=cross<T>(c)*h.template block<3,3>(0,3)*d*invNLen;
      //gc,c
      h.template block<3,3>(6,6)=(bnn*g.template segment<3>(6).transpose()+g.template segment<3>(6)*bnn.transpose())*invNLen;
      h.template block<3,3>(6,6)+=cross<T>(b)*h.template block<3,3>(0,6)*d*invNLen;
      //gb,c
      h.template block<3,3>(3,6)=(-cross<T>(a)-d*cross<T>(c)*h.template block<3,3>(0,6)+d*cross<T>(nn)-cnn*g.template segment<3>(6).transpose()+g.template segment<3>(3)*bnn.transpose())*invNLen;
      h.template block<3,3>(6,3)=h.template block<3,3>(3,6).transpose();
    }
    if(G)
      *G=getIEE()*g*nSgn;
    if(H)
      *H=(getIEE()*h*getIEE().transpose())*nSgn;
  } else if(_feat[0]>=0&& _feat[1]>=0) {
    Vec3T vv[2]= {v[_feat[0]+0],v[_feat[1]+2]};
    Vec3T g;
    Mat3T h;
    if(G || H)
      g=(vv[0]-vv[1])/d;
    if(H)
      h=(Mat3T::Identity()-g*g.transpose())/d;
    if(G) {
      G->setZero();
      G->template segment<3>((_feat[0]+0)*3)=g;
      G->template segment<3>((_feat[1]+2)*3)=-g;
    }
    if(H) {
      H->setZero();
      H->template block<3,3>((_feat[0]+0)*3,(_feat[0]+0)*3)=h;
      H->template block<3,3>((_feat[0]+0)*3,(_feat[1]+2)*3)=-h;
      H->template block<3,3>((_feat[1]+2)*3,(_feat[0]+0)*3)=-h;
      H->template block<3,3>((_feat[1]+2)*3,(_feat[1]+2)*3)=h;
    }
  } else {
    //closest _feature is edge
    Eigen::Matrix<int,3,1>
    _featE=_feat[0]==-1?Eigen::Matrix<int,3,1>(2+_feat[1],0,1):Eigen::Matrix<int,3,1>(_feat[0],2,3);
    //std::shared_ptr<FEMVertex> VE[3]={_EE._v[_featE[0]],_EE._v[_featE[1]],_EE._v[_featE[2]]};
    std::array<int,3> id= {_featE[0],_featE[1],_featE[2]};
    Vec3T vE[3]= {v[_featE[0]],v[_featE[1]],v[_featE[2]]};
    Vec3T a=vE[0]-vE[1],b=vE[0]-vE[2],c=vE[1]-vE[2];
    Vec3T n=a.cross(b);
    T nLen=n.norm(),invNLen=1/nLen,invCLen=1/c.norm();
    Vec3T nn=n*invNLen,cn=c*invCLen;
    Mat3T projN;
    Vec9T g,GSum;
    Mat9T h,HSum;
    if(G || H) {
      g.template segment<3>(0)=b.cross(nn)*invCLen;
      g.template segment<3>(3)=-a.cross(nn)*invCLen;
      g.template segment<3>(6)=-cn*nLen*invCLen*invCLen;
    }
    if(H) {
      projN=Mat3T::Identity()-nn*nn.transpose();
      //ga,a
      h.template block<3,3>(0,0)=-cross<T>(b)*projN*cross<T>(b)*invNLen*invCLen;
      //ga,b
      h.template block<3,3>(0,3)=cross<T>(b)*projN*cross<T>(a)*invNLen*invCLen-cross<T>(nn)*invCLen;
      h.template block<3,3>(3,0)=h.template block<3,3>(0,3).transpose();
      //ga,c
      h.template block<3,3>(0,6)=-g.template segment<3>(0)*cn.transpose()*invCLen;
      h.template block<3,3>(6,0)=h.template block<3,3>(0,6).transpose();
      //gb,b
      h.template block<3,3>(3,3)=-cross<T>(a)*projN*cross<T>(a)*invNLen*invCLen;
      //gb,c
      h.template block<3,3>(3,6)=-g.template segment<3>(3)*cn.transpose()*invCLen;
      h.template block<3,3>(6,3)=h.template block<3,3>(3,6).transpose();
      //gc,c
      h.template block<3,3>(6,6)=(3*cn*cn.transpose()-Mat3T::Identity())*nLen*invCLen*invCLen*invCLen;
    }
    if(G) {
      G->setZero();
      GSum=VTDistanceEnergy<T>::getIE()*g;
      for(int i=0; i<(int) id.size(); ++i)
        G->template segment<3>(id[i]*3)=GSum.template segment<3>(i*3);
    }
    if(H) {
      H->setZero();
      HSum=VTDistanceEnergy<T>::getIE()*h*VTDistanceEnergy<T>::getIE().transpose();
      for(int i=0; i<(int) id.size(); ++i)
        for(int j=0; j<(int) id.size(); ++j)
          H->template block<3,3>(id[i]*3,id[j]*3)=HSum.template block<3,3>(i*3,j*3);
    }
  }
  return true;
}
template <typename T>
const typename EEDistanceEnergy<T>::SMatT& EEDistanceEnergy<T>::getIEE() {
  static SMatT I;
  OMP_CRITICAL_{
    if(I.size()==0) {
      STrips trips;
      //a
      addBlockId<T>(trips,0,0,3,1);
      addBlockId<T>(trips,6,0,3,-1);
      //b
      addBlockId<T>(trips,0,3,3,1);
      addBlockId<T>(trips,3,3,3,-1);
      //c
      addBlockId<T>(trips,6,6,3,1);
      addBlockId<T>(trips,9,6,3,-1);
      I.resize(12,9);
      I.setFromTriplets(trips.begin(),trips.end());
    }
  }
  return I;
}
template <typename T>
bool EEDistanceEnergy<T>::debug(const Vec3T edge1[2],const Vec3T edge2[2],T perturbRange,const Eigen::Matrix<char,2,1>& feat) {
  _edge1[0]=edge1[0];
  _edge1[1]=edge1[1];
  _edge2[0]=edge2[0];
  _edge2[1]=edge2[1];
  T E,E2;
  Mat12T H;
  Vec12T G,G2,dx;
  DEFINE_NUMERIC_DELTA_T(T)
  dx=Vec12T::Random();
  _edge1[0]+=Vec::Random(3)*perturbRange;
  _edge1[1]+=Vec::Random(3)*perturbRange;
  _edge2[0]+=Vec::Random(3)*perturbRange;
  _edge2[1]+=Vec::Random(3)*perturbRange;
  if(!eval(&E,&G,&H) || E==0 || !isfinite(E))
    return false;
  _edge1[0]+=dx.template segment<3>(0*3)*DELTA;
  _edge1[1]+=dx.template segment<3>(1*3)*DELTA;
  _edge2[0]+=dx.template segment<3>(2*3)*DELTA;
  _edge2[1]+=dx.template segment<3>(3*3)*DELTA;
  if(!eval(&E2,&G2,NULL))
    return false;
  if(_feat!=feat)
    return false;
  std::cout << "Feat=(" << (int)_feat[0] << "," << (int)_feat[1] << ")" << std::endl;
  DEBUG_GRADIENT("dE",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("dG",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
  return true;
}
template <typename T>
void EEDistanceEnergy<T>::debug(T perturbRange,const Eigen::Matrix<char,2,1>& feat) {
  while(true) {
    Vec3T edge1[2]= {Vec3T::Random(),Vec3T::Random()};
    Vec3T edge2[2]= {Vec3T::Random(),Vec3T::Random()};
    if(debug(edge1,edge2,perturbRange,feat))
      break;
  }
}
//VTBarrierEnergy
template <typename T,typename PFunc>
VTBarrierEnergy<T,PFunc>::VTBarrierEnergy(const Vec3T triA,const Vec3T& triB,const Vec3T& triC,const Vec3T& point,const PFunc& p,T d0,T coef)
  :VTDistanceEnergy<T>(triA,triB,triC,point),_p(p),_d0(d0),_coef(coef) {}
template <typename T,typename PFunc>
VTBarrierEnergy<T,PFunc>::VTBarrierEnergy(const Vec3T tri[3],const Vec3T& point,const PFunc& p,T d0,T coef)
  :VTDistanceEnergy<T>(tri,point),_p(p),_d0(d0),_coef(coef) {}
template <typename T,typename PFunc>
bool VTBarrierEnergy<T,PFunc>::eval(T* E,Vec12T* G,Mat12T* H) {
  T e,distance=0,D=0,DD=0;
  if(!VTDistanceEnergy<T>::eval(&distance,G,H))
    return false;
  e=_p.eval(distance,G?&D:NULL,H?&DD:NULL,_d0,_coef);
  if(E)
    *E=e;
  if(H)
    *H=*H*D+DD*(*G)*G->transpose();
  if(G)
    *G*=D;
  return true;
}
//EEBarrierEnergy
template <typename T,typename PFunc>
EEBarrierEnergy<T,PFunc>::EEBarrierEnergy(const Vec3T edge1A,const Vec3T edge1B,const Vec3T edge2A,const Vec3T edge2B,bool JTJApprox,const PFunc& p,T d0,T coef,T mollifierCoef)
  :EEDistanceEnergy<T>(edge1A,edge1B,edge2A,edge2B,JTJApprox),_p(p),_d0(d0),_coef(coef),_mollifierCoef(mollifierCoef) {}
template <typename T,typename PFunc>
EEBarrierEnergy<T,PFunc>::EEBarrierEnergy(const Vec3T edge1[2],const Vec3T edge2[2],bool JTJApprox,const PFunc& p,T d0,T coef,T mollifierCoef)
  :EEDistanceEnergy<T>(edge1[0],edge1[1],edge2[0],edge2[1],JTJApprox),_p(p),_d0(d0),_coef(coef),_mollifierCoef(mollifierCoef) {}
template <typename T,typename PFunc>
bool EEBarrierEnergy<T,PFunc>::eval(T* E,Vec12T* G,Mat12T* H) {
  Vec12T GM;
  Mat12T HM;
  T e,eM=0,distance=0,D=0,DD=0;
  if(!EEDistanceEnergy<T>::eval(&distance,G,H))
    return false;
  e=_p.eval(distance,G?&D:NULL,H?&DD:NULL,_d0,_coef);
  //mollifier
  if(mollifier(&eM,G?&GM:NULL,H?&HM:NULL)) {
    if(E)
      *E=e*eM;
    if(H)
      *H=(*H*D+DD*(*G)*G->transpose())*eM+HM*e+(*G*GM.transpose()+GM*G->transpose())*D;
    if(G)
      *G=*G*D*eM+GM*e;
  } else {
    if(E)
      *E=e;
    if(H)
      *H=*H*D+DD*(*G)*G->transpose();
    if(G)
      *G=*G*D;
  }
  return true;
}
template <typename T,typename PFunc>
bool EEBarrierEnergy<T,PFunc>::mollifier(T* E,Vec12T* G,Mat12T* H) const {
  if(_mollifierCoef<=0)
    return false;
  Vec3T v01=_edge1[0]-_edge1[1];
  Vec3T v23=_edge2[0]-_edge2[1];
  //compute energy
  Vec3T crossVar=v01.cross(v23);
  T c=crossVar.squaredNorm();
  if(c>=_mollifierCoef)
    return false;
  if(E)
    (*E)=c*2/_mollifierCoef-c*c/(_mollifierCoef*_mollifierCoef);
  if(G||H) {
    Mat3T c01=cross<T>(v01);
    Mat3T c23=cross<T>(v23);
    Vec3T D01=c01*crossVar;
    Vec3T D23=c23*crossVar;
    Vec3T D[4]= {D23,-D23,-D01,D01};
    T coefG=2/_mollifierCoef*(1-c/_mollifierCoef)*2;
    if(G)
      for(int d=0; d<4; d++)
        G->template segment<3>(d*3)=D[d]*coefG;
    if(H) {
      Mat3T HBlk[2][2];
      HBlk[0][0]=c23*c23.transpose();
      HBlk[1][1]=c01*c01.transpose();
      HBlk[0][1]=c23*c01-cross<T>(crossVar);
      HBlk[1][0]=HBlk[0][1].transpose();
      T coefH=-2/(_mollifierCoef*_mollifierCoef)*2*2;
      for(int d=0; d<4; d++)
        for(int d2=0; d2<4; d2++) {
          Mat3T blk=D[d]*D[d2].transpose()*coefH;
          if(!_JTJApprox)
            blk+=HBlk[d/2][d2/2]*coefG*(d%2?-1:1)*(d2%2?-1:1);
          H->template block<3,3>(d*3,d2*3)=blk;
        }
    }
  }
  return true;
}
//instance
template class VTDistanceEnergy<FLOAT>;
template class EEDistanceEnergy<FLOAT>;
template class VTBarrierEnergy<FLOAT,Px>;
template class EEBarrierEnergy<FLOAT,Px>;
template class VTBarrierEnergy<FLOAT,Logx>;
template class EEBarrierEnergy<FLOAT,Logx>;
template class VTBarrierEnergy<FLOAT,CLogx>;
template class EEBarrierEnergy<FLOAT,CLogx>;
template class VTBarrierEnergy<FLOAT,InvQuadraticx>;
template class EEBarrierEnergy<FLOAT,InvQuadraticx>;
template class VTBarrierEnergy<FLOAT,Cubicx>;
template class EEBarrierEnergy<FLOAT,Cubicx>;
}

