#include "ConvexHullDistanceEnergy.h"
#include "TrajectorySIPEnergy.h"
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
//CCDistanceEnergy
template <typename T>
CCDistanceEnergy<T>::CCDistanceEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2):_p1(p1),_p2(p2) {}
template <typename T>
bool CCDistanceEnergy<T>::eval(T* E,Vec12T* G,Mat12T* H) {
  ASSERT_MSG(!G && !H,"Convex hull distance function is not differentiable!")
  T dist=GJKPolytope<T>::distance(_p1,_p2,NULL,NULL);
  if(E)
    *E=dist;
  return true;
}
//CCBarrierEnergy
template <typename T,typename PFunc,typename TH>
CCBarrierEnergy<T,PFunc,TH>::CCBarrierEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit)
  :CCDistanceEnergy<T>(p1,p2),_p(p),_d0(d0),_d0Half(d0/2),_coef(coef),_implicit(implicit),_output(false) {
  //globalVss1
  if(grad && p1.jid()>=0)
    _globalVss1=&(grad->_globalVss[p1.jid()]);
  else {
    ASSERT_MSG(grad || p1.jid()==-1,"CollisionGradInfo must be provided to extract globalVss!")
    _globalVss1=NULL;
  }
  //globalVss2
  if(grad && p2.jid()>=0)
    _globalVss2=&(grad->_globalVss[p2.jid()]);
  else {
    ASSERT_MSG(grad || p2.jid()==-1,"CollisionGradInfo must be provided to extract globalVss!")
    _globalVss2=NULL;
  }
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::update(Vec4T* res) {
  TH E;
  Vec4TH G;
  Mat4TH H;
  bool succ=optimize(_x,E,G,H);
  if(res)
    *res=_x.template cast<T>();
  return succ;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::initialize(Vec4T* res,const ArticulatedBody* body) {
  //initial guess
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
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
  if(res)
    *res=_x.template cast<T>();
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::initialize(const Vec4T& x0) {
  _x=x0.template cast<TH>();
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta) {
  TH E2;
  Vec4TH G2;
  Mat4TH H2;
  if(_implicit) {
    if(!initialize(NULL,body)) {
      std::cout<<"Cannot initialize"<<std::endl;
      return false;
    }
    //optimize
    if(!optimize(_x,E2,G2,H2)) {
      std::cout<<"Cannot optimize"<<std::endl;
      return false;
    }
    if(E)
      *E=_coef*(T)E2;
    if(E2==0)
      return true;
  } else {
    //initialize();
    //compute energy
    if(!energy(_x,E2,NULL,NULL))
      return false;
    if(E)
      *E=_coef*(T)E2;
    if(E2==0)
      return true;
  }
  //compute gradient and hessian
  if(body && grad) {
    Mat3XT DTG;
    Vec4T G=G2.template cast<T>();
    Mat4T H=H2.template cast<T>();
    bool hessian=grad && grad->_HTheta.size()>0;
    computeDTGH(*body,*grad,_x.template cast<T>(),
                hessian? &G: NULL,hessian? &H: NULL);
    /*if(GTheta) {
      GTheta->setZero(body->nrDOF());
      grad->_info.DTG(*body,mapM(DTG=grad->_DTG),mapV(*GTheta));
    }
    if(HTheta) {
      *HTheta=grad->_HTheta;
      grad->_info.toolB(*body,mapM(DTG=grad->_DTG),[&](int r,int c,T val) {
        (*HTheta)(r,c)+=val;
      });
    }*/
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0,bool output) {
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
      CCBarrierEnergy<T,PFunc,TH> e(pass? p: p2,pass? p2: p,barrier,d0,&info,coef,implicit);
      e.setOutput(output);
      if(!implicit)
        e.initialize(NULL,&body);
      if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
        continue;
      //evaluate again
      info2.reset(body,x+dx*DELTA);
      p2=GJKPolytope<T>(JID,body,info2);
      CCBarrierEnergy<T,PFunc,TH> e2(pass? p: p2,pass? p2: p,barrier,d0,&info2,coef,implicit);
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
void CCBarrierEnergy<T,PFunc,TH>::debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
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
    CCBarrierEnergy<T,PFunc,TH> e(p,p2,barrier,d0,&info,coef,implicit);
    e.setOutput(output);
    if(!implicit)
      e.initialize(NULL,&body);
    if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
      continue;
    //evaluate again
    info2.reset(body,x+dx*DELTA);
    p=GJKPolytope<T>(JID,body,info2);
    p2=GJKPolytope<T>(JID2,body,info2);
    CCBarrierEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info2,coef,implicit);
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
void CCBarrierEnergy<T,PFunc,TH>::setOutput(bool output) {
  _output=output;
}
//helper
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::energyP(const Vec3TH& v,const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const {
  TH e,D=0,DD=0,d=v.dot(x.template segment<3>(0))-x[3];
  Vec4TH x2=Vec4TH(v[0],v[1],v[2],-1);
  e=_p.template eval<TH>(d,G? &D: NULL,H? &DD: NULL,(TH)_d0Half,1);
  if(!isfinite(e))
    return false;
  if(e==0)
    return true;
  E+=e;
  if(G)
    *G+=x2*D;
  if(H)
    *H+=x2*x2.transpose()*DD;
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::energyN(const Vec3TH& v,const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const {
  TH e,D=0,DD=0,d=x[3]-v.dot(x.template segment<3>(0));
  Vec4TH x2=Vec4TH(-v[0],-v[1],-v[2],1);
  e=_p.template eval<TH>(d,G? &D: NULL,H? &DD: NULL,(TH)_d0Half,1);
  if(!isfinite(e))
    return false;
  if(e==0)
    return true;
  E+=e;
  if(G)
    *G+=x2*D;
  if(H)
    *H+=x2*x2.transpose()*DD;
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::debugEnergyP(const Vec3TH& v,const Vec4TH& x,T perturbRange) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Vec4TH G,G2,dx,x2;
  dx.setRandom();
  Mat4TH H,H2;
  // Debug energyP
  E=E2=0;
  G.setZero();
  G2.setZero();
  H.setZero();
  H2.setZero();
  energyP(v,x,E,&G,&H);
  x2=x+dx*DELTA;
  energyP(v,x2,E2,&G2,&H2);
  if(E==0 && E2==0)return;
  DEBUG_GRADIENT("energyP E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("energyP G",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::debugEnergyN(const CCBarrierEnergy::Vec3TH& v,const CCBarrierEnergy::Vec4TH& x,T perturbRange) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Vec4TH G,G2,dx,x2;
  dx.setRandom();
  Mat4TH H,H2;

  // Debug energyN
  E=E2=0;
  G.setZero();
  G2.setZero();
  H.setZero();
  H2.setZero();
  energyN(v,x,E,&G,&H);
  x2=x+dx*DELTA;
  energyN(v,x2,E2,&G2,&H2);
  if(E==0 && E2==0)return;
  DEBUG_GRADIENT("energyN E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("energyN G",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::energy(const Vec4TH& x,TH& E,Vec4TH* G,Mat4TH* H) const {
  E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  //positive
  if(_globalVss1) {
    for(int c=0; c<_globalVss1->cols(); c++) {
//      debugEnergyP(_globalVss1->col(c).template cast<TH>(),x,1e-9);
      if(!energyP(_globalVss1->col(c).template cast<TH>(),x,E,G,H))
        return false;
    }
  } else {
    for(int c=0; c<(int)_p1.mesh()->vss().size(); c++) {
//      debugEnergyP(_p1.mesh()->vss()[c].template cast<TH>(),x,1e-9);
      if(!energyP(_p1.mesh()->vss()[c].template cast<TH>(),x,E,G,H))
        return false;
    }
  }
  //negative
  if(_globalVss2) {
    for(int c=0; c<_globalVss2->cols(); c++) {
//      debugEnergyN(_globalVss2->col(c).template cast<TH>(),x,1e-9);
      if(!energyN(_globalVss2->col(c).template cast<TH>(),x,E,G,H))
        return false;
    }
  } else {
    for(int c=0; c<(int)_p2.mesh()->vss().size(); c++) {
//      debugEnergyN(_p2.mesh()->vss()[c].template cast<TH>(),x,1e-9);
      if(!energyN(_p2.mesh()->vss()[c].template cast<TH>(),x,E,G,H))
        return false;
    }
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::optimize(Vec4TH& x,TH& E,Vec4TH& G,Mat4TH& H) const {
  Vec4TH invHN,N,D,G2,GPrj,x2;
  Eigen::FullPivLU<Mat4TH> invH;
  Mat4TH id=Mat4TH::Identity(),H2;
  TH E2,alpha=1,alphaDec=0.5,alphaInc=3.0,alphaMax=1e10;
  //assemble
  x.template segment<3>(0).normalize();//project
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
    N=Vec4TH(x[0],x[1],x[2],0);
    invH.compute(H+id*alpha);
    invHN=invH.solve(N);
    D=invH.solve(G)-invHN*invHN.dot(G)/invHN.dot(N);
    //termination
    GPrj=G-G.dot(N)*N;
    if(GPrj.cwiseAbs().maxCoeff()<=Epsilon<TH>::defaultEps())
      break;
    //test
    x2=x-D;
    x2.template segment<3>(0).normalize();
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
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::debugEnergy(const Vec4TH& x) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Mat4TH H;
  Vec4TH G,G2,dx;
  dx.setRandom();
  energy(x,E,&G,&H);
  energy(x+dx*DELTA,E2,&G2,NULL);
  DEBUG_GRADIENT("E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("G",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
//gradient
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::sensitivity(const Vec4T& x,const Vec4T& G,const Mat4T& H,Mat4T& S) const {
  S=H;
  S.template block<3,3>(0,0).diagonal().array()-=G.template segment<3>(0).dot(x.template segment<3>(0));
  S=S.inverse().eval();
  Vec4T SN=S.template block<4,3>(0,0)*x.template segment<3>(0);
  S-=SN*SN.transpose()/SN.template segment<3>(0).dot(x.template segment<3>(0));
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::energyPDTGH(const Vec3T& v,const Vec3T& vl,const Vec4T& x,Mat3X4T* DTG,
    const Vec3T& Rvl,Mat3X4T* LRH,Mat3X4T* wLRH,
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt) const {
  Mat3T nntDD,C=cross<T>(Rvl);
  Mat3X4T DDEDXDTheta;
  T D=0,DD=0,d=v.dot(x.template segment<3>(0))-x[3];
  T e=_p.template eval<T>(d,&D,&DD,_d0Half,1);
  if(e==0)
    return;
  if(DTG)
    SIPEnergy<T>::template parallelAdd<3,4>(*DTG,0,0,SIPEnergy<T>::computeDTG(x.template segment<3>(0)*D*_coef,vl));
  if(LRH && wLRH) {
    ROT(DDEDXDTheta)=(x.template segment<3>(0)*v.transpose())*DD+Mat3T::Identity()*D;
    CTR(DDEDXDTheta)=-x.template segment<3>(0)*DD;
    SIPEnergy<T>::template parallelAdd<3,4>(*LRH,0,0,DDEDXDTheta);
    SIPEnergy<T>::template parallelAdd<3,4>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(Mww && Mtw && Mwt && Mtt) {
    nntDD=x.template segment<3>(0)*x.template segment<3>(0).transpose()*DD*_coef;
    SIPEnergy<T>::template parallelAdd<3,3>(*Mww,0,0,C*nntDD*C.transpose());
    SIPEnergy<T>::template parallelAdd<3,3>(*Mtw,0,0,nntDD*C.transpose());
    SIPEnergy<T>::template parallelAdd<3,3>(*Mwt,0,0,C*nntDD);
    SIPEnergy<T>::template parallelAdd<3,3>(*Mtt,0,0,nntDD);
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::energyNDTGH(const Vec3T& v,const Vec3T& vl,const Vec4T& x,Mat3X4T* DTG,
    const Vec3T& Rvl,Mat3X4T* LRH,Mat3X4T* wLRH,
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt) const {
  Mat3T nntDD,C=cross<T>(Rvl);
  Mat3X4T DDEDXDTheta;
  T D=0,DD=0,d=x[3]-v.dot(x.template segment<3>(0));
  T e=_p.template eval<T>(d,&D,&DD,_d0Half,1);
  if(e==0)
    return;
  if(DTG)
    SIPEnergy<T>::template parallelAdd<3,4>(*DTG,0,0,SIPEnergy<T>::computeDTG(-x.template segment<3>(0)*D*_coef,vl));
  if(LRH && wLRH) {
    ROT(DDEDXDTheta)=(x.template segment<3>(0)*v.transpose())*DD-Mat3T::Identity()*D;
    CTR(DDEDXDTheta)=-x.template segment<3>(0)*DD;
    SIPEnergy<T>::template parallelAdd<3,4>(*LRH,0,0,DDEDXDTheta);
    SIPEnergy<T>::template parallelAdd<3,4>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(Mww && Mtw && Mwt && Mtt) {
    nntDD=x.template segment<3>(0)*x.template segment<3>(0).transpose()*DD*_coef;
    SIPEnergy<T>::template parallelAdd<3,3>(*Mww,0,0,C*nntDD*C.transpose());
    SIPEnergy<T>::template parallelAdd<3,3>(*Mtw,0,0,nntDD*C.transpose());
    SIPEnergy<T>::template parallelAdd<3,3>(*Mwt,0,0,C*nntDD);
    SIPEnergy<T>::template parallelAdd<3,3>(*Mtt,0,0,nntDD);
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::computeDTGH(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec4T& x,const Vec4T* G,const Mat4T* H) const {
  Mat3X4T wDDEDXDTheta1,tDDEDXDTheta1,DTG;
  Mat3X4T wDDEDXDTheta2,tDDEDXDTheta2;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat4T S,Hxx;
  if(G && H && _implicit) {
    sensitivity(x,*G,*H,S);
    Hxx=-S*_coef;
  }
  //positive
  if(_globalVss1) {
    if(G && H) {
      wDDEDXDTheta1.setZero();
      tDDEDXDTheta1.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p1.jid());
    }
    //stage 1
    DTG.setZero();
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    body.joint(_p1.jid())._mesh->getMesh(vss,iss);
    std::shared_ptr<MeshExact> local;
    local.reset(new MeshExact(vss,iss,false));
    //std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p1.jid())._mesh);
    for(int c=0; c<_globalVss1->cols(); c++)
      energyPDTGH(_globalVss1->col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta1: NULL,&wDDEDXDTheta1,
                  (G && H)? &Mww: NULL,&Mtw,&Mwt,&Mtt);
    //TRANSI(info._DTG,_p1.jid())+=DTG;
    for(int r=0; r<3; r++)

      for(int c=0; c<4; c++)

        SIPEnergy<T>::parallelAdd(info._DTG(r,c+_p1.jid()*4),DTG(r,c));
    //stage 2
    if(G && H) {
      if(_implicit) {
        Mww+=wDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        Mtw+=tDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        Mwt+=wDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        Mtt+=tDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
      }
      contractHTheta(_p1.jid(),_p1.jid(),body,info,Mww,Mtw,Mwt,Mtt);
    }
  }
  //negative
  if(_globalVss2) {
    if(G && H) {
      wDDEDXDTheta2.setZero();
      tDDEDXDTheta2.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p2.jid());
    }
    //stage 1
    DTG.setZero();
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    body.joint(_p2.jid())._mesh->getMesh(vss,iss);
    std::shared_ptr<MeshExact> local;
    local.reset(new MeshExact(vss,iss,false));
    
    //std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_globalVss2->cols(); c++)
      energyNDTGH(_globalVss2->col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta2: NULL,&wDDEDXDTheta2,
                  (G && H)? &Mww: NULL,&Mtw,&Mwt,&Mtt);
    //TRANSI(info._DTG,_p2.jid())+=DTG;
    for(int r=0; r<3; r++)

      for(int c=0; c<4; c++)

        SIPEnergy<T>::parallelAdd(info._DTG(r,c+_p2.jid()*4),DTG(r,c));
    //stage 2
    if(G && H) {
      if(_implicit) {
        Mww+=wDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        Mtw+=tDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        Mwt+=wDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        Mtt+=tDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
      }
      contractHTheta(_p2.jid(),_p2.jid(),body,info,Mww,Mtw,Mwt,Mtt);
    }
  }
  //positive/negative
  if(_globalVss1 && _globalVss2 && _implicit) {
    //stage 2
    if(G && H) {
      Mww=wDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose();
      Mtw=tDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose();
      Mwt=wDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose();
      Mtt=tDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose();
      contractHTheta(_p1.jid(),_p2.jid(),body,info,Mww,Mtw,Mwt,Mtt);
      contractHTheta(_p2.jid(),_p1.jid(),body,info,Mww.transpose(),Mwt.transpose(),Mtw.transpose(),Mtt.transpose());
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::contractHTheta
(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,
 const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const {
  info._info.JRCSparse(body,kL,[&](int row,const Vec3T& JRL) {
    info._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      SIPEnergy<T>::parallelAdd(info._HTheta(row,col),JRL.dot(Mww*JRR));
    },[&](int col,const Vec3T& JCR) {
      SIPEnergy<T>::parallelAdd(info._HTheta(row,col),JRL.dot(Mwt*JCR));
    });
  },[&](int row,const Vec3T& JCL) {
    info._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      SIPEnergy<T>::parallelAdd(info._HTheta(row,col),JCL.dot(Mtw*JRR));
    },[&](int col,const Vec3T& JCR) {
      SIPEnergy<T>::parallelAdd(info._HTheta(row,col),JCL.dot(Mtt*JCR));
    });
  });
}
//instance
template class CCDistanceEnergy<FLOAT>;
template class CCBarrierEnergy<FLOAT,Px>;
template class CCBarrierEnergy<FLOAT,Logx>;
template class CCBarrierEnergy<FLOAT,CLogx>;
template class CCBarrierEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierEnergy<FLOAT,Cubicx>;
}
