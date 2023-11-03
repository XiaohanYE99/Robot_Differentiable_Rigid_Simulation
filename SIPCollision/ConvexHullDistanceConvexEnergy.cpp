#include "ConvexHullDistanceConvexEnergy.h"
#include "TrajectorySIPEnergy.h"
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
//CCBarrierEnergy
template <typename T,typename PFunc,typename TH>
CCBarrierConvexEnergy<T,PFunc,TH>::CCBarrierConvexEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit)
  :CCBarrierEnergy<T,PFunc,TH>(p1,p2,p,d0,grad,coef,implicit) {}
template <typename T,typename PFunc,typename TH>
bool CCBarrierConvexEnergy<T,PFunc,TH>::initialize(Vec4T* res,const ArticulatedBody* body) {
  //initial guess
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
  //std::cout<<"dist= "<<dist<<std::endl;
  if(dist<=_d0) {
    OMP_CRITICAL_{
      std::cout<<"Penetrated in CCBarrierEnergy! d="<<dist<<std::endl;
    }
    return false;
  }
  //normal,project
  Vec3T n=Vec3T(p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]).normalized();
  Vec2T r1=GJKPolytope<T>::project(_p1,n);
  Vec2T r2=GJKPolytope<T>::project(_p2,n);
  if(r1[1]<r2[0]-_d0) {
    _x=Vec4TH((TH)-n[0],(TH)-n[1],(TH)-n[2],(TH)-(r1[1]+r2[0])/2);
  } else {
    if(r2[1]>=r1[0]-_d0) {
      OMP_CRITICAL_{
        std::cout<<"Strange error, GJK predicts separation by "<<dist
                 <<" but the projected intervals are not separate, "
                 <<"range=["<<r1[0]<<","<<r1[1]<<"] and ["<<r2[0]<<","<<r2[1]<<"]"<<std::endl;
        GJKPolytope<T>::writeConfigVTK("GJKConfig.vtk",_p1,_p2,p1,p2,body);
        GJKPolytope<T>::writeConfig("GJKConfig.dat",_p1,_p2);
        //ASSERT(false)
      }
      return false;
    }
    _x=Vec4TH((TH)n[0],(TH)n[1],(TH)n[2],(TH)(r2[1]+r1[0])/2);
  }
  //SVM margin
  TH alpha=(_d0Half+1)/(dist/2+1);
  _x*=alpha;
  if(res)
    *res=_x.template cast<T>();
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierConvexEnergy<T,PFunc,TH>::debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0,bool output) {
  T E,E2,coef;
  MatT HTheta;
  PFunc barrier;
  Vec GTheta,GTheta2;
  GJKPolytope<T> p2;
  CollisionGradInfo<T> info,info2;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  for(int pass=0; pass<2; pass++)
    while(true) {
      coef=rand()/(T)RAND_MAX;
      Vec x=Vec::Random(body.nrDOF());
      Vec dx=Vec::Random(body.nrDOF());
      //evaluate energy/gradient
      info.reset(body,x);
      p2=GJKPolytope<T>(JID,body,info);
      //p.writeVTK("poly1.vtk",&body);
      //p2.writeVTK("poly2.vtk",&body);
      CCBarrierConvexEnergy<T,PFunc,TH> e(pass? p: p2,pass? p2: p,barrier,d0,&info,coef,implicit);
      e.setOutput(output);
      if(!implicit)
        e.initialize(NULL,&body);
      if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
        continue;
      //evaluate again
      info2.reset(body,x+dx*DELTA);
      p2=GJKPolytope<T>(JID,body,info2);
      CCBarrierConvexEnergy<T,PFunc,TH> e2(pass? p: p2,pass? p2: p,barrier,d0,&info2,coef,implicit);
      e2.setOutput(output);
      if(!implicit)
        e2.initialize(e._x.template cast<T>());
      if(!e2.eval(&E2,&body,&info2,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dE"+std::string(implicit? "Implicit": "Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG"+std::string(implicit? "Implicit": "Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
      break;
    }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierConvexEnergy<T,PFunc,TH>::debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
  T E,E2,coef;
  MatT HTheta;
  PFunc barrier;
  Vec GTheta,GTheta2;
  GJKPolytope<T> p,p2;
  CollisionGradInfo<T> info,info2;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  while(true) {
    coef=rand()/(T)RAND_MAX;
    Vec x=Vec::Random(body.nrDOF());
    Vec dx=Vec::Random(body.nrDOF());
    //evaluate energy/gradient
    info.reset(body,x);
    p=GJKPolytope<T>(JID,body,info);
    p2=GJKPolytope<T>(JID2,body,info);
    //p.writeVTK("poly1.vtk",&body);
    //p2.writeVTK("poly2.vtk",&body);
    CCBarrierConvexEnergy<T,PFunc,TH> e(p,p2,barrier,d0,&info,coef,implicit);
    e.setOutput(output);
    T dist=GJKPolytope<T>::distance(p,p2,NULL,NULL);
    std::cout<<dist<<std::endl;
    if(!implicit)
      e.initialize(NULL,&body);
    if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
      continue;
    //evaluate again
    info2.reset(body,x+dx*DELTA);
    p=GJKPolytope<T>(JID,body,info2);
    p2=GJKPolytope<T>(JID2,body,info2);
    CCBarrierConvexEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info2,coef,implicit);
    e2.setOutput(output);
    if(!implicit)
      e2.initialize(e._x.template cast<T>());
    if(!e2.eval(&E2,&body,&info2,&GTheta2,NULL))
      continue;
    DEBUG_GRADIENT("dE"+std::string(implicit? "Implicit": "Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("dG"+std::string(implicit? "Implicit": "Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    break;
  }
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierConvexEnergy<T,PFunc,TH>::energy(const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const {
  if(!CCBarrierEnergy<T,PFunc,TH>::energy(x,E,G,H))
    return false;
  //normal length<1
  TH e,D=0,DD=0,len=x.template segment<3>(0).norm();
  Vec3TH n=x.template segment<3>(0)/len;
  e=_p.template eval<TH>(1-len,G? &D: NULL,H? &DD: NULL,0,1);
  if(!isfinite(e))
    return false;
  if(e==0)
    return true;
  E+=e;
  if(G)
    G->template segment<3>(0)+=-D*n;
  if(H)
    H->template block<3,3>(0,0)+=DD*n*n.transpose()+D*(n*n.transpose()-Mat3TH::Identity())/len;
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierConvexEnergy<T,PFunc,TH>::optimize(Vec4TH& x,TH& E,Vec4TH& G,Mat4TH& H) const {
  Vec4TH D,G2,x2;
  Eigen::FullPivLU<Mat4TH> invH;
  Mat4TH id=Mat4TH::Identity(),H2;
  TH E2,alpha=1,alphaDec=0.5,alphaInc=3.0,alphaMax=1e10;
  //assemble
  //debugEnergy(x);
  if(!energy(x,E,&G,&H))
    return false;
  //main loop
  alpha*=std::max<TH>(1,H.cwiseAbs().maxCoeff());
  alphaMax*=std::max<TH>(1,H.cwiseAbs().maxCoeff());
//  int i;
//  for(i=0;i<1e5;++i){
  while(alpha<alphaMax) {
    if(E==0)
      break;
    //search direction
    invH.compute(H+id*alpha);
    D=invH.solve(G);
    //termination
    if(G.cwiseAbs().maxCoeff()<=Epsilon<TH>::defaultEps())
      break;
    //test
    x2=x-D;
    if(energy(x2,E2,&G2,&H2) && E2<E) {
      alpha=std::max<TH>(alpha*alphaDec,Epsilon<TH>::defaultEps());
      x=x2;
      E=E2;
      G=G2;
      H=H2;
    } else {
      alpha*=alphaInc;
    }
  }
  return true;
}
//gradient
template <typename T,typename PFunc,typename TH>
void CCBarrierConvexEnergy<T,PFunc,TH>::sensitivity(const Vec4T&,const Vec4T&,const Mat4T& H,Mat4T& S) const {
  S=H.inverse();
}
//instance
template class CCBarrierConvexEnergy<FLOAT,Px>;
template class CCBarrierConvexEnergy<FLOAT,Logx>;
template class CCBarrierConvexEnergy<FLOAT,CLogx>;
template class CCBarrierConvexEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierConvexEnergy<FLOAT,Cubicx>;
}
