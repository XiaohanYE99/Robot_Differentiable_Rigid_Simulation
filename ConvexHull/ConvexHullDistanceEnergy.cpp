#include "ConvexHullDistanceEnergy.h"
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
  :CCDistanceEnergy<T>(p1,p2),_p(p),_d0(d0),_d0Half(d0/2),_coef(coef),_implicit(implicit),_output(false),_grad(grad) {}
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
Eigen::Matrix<T,4,1> CCBarrierEnergy<T,PFunc,TH>::getX() const {
  return _x.template cast<T>();
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::initialize(Vec4T* res,const ArticulatedBody* body) {
  //initial guess
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
  if(dist<=_d0) {
    /*OMP_CRITICAL_{
      std::cout<<"Penetrated in CCBarrierEnergy! d="<<dist<<std::endl;
    }*/
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
      /*OMP_CRITICAL_{
        std::cout<<"Strange error, GJK predicts separation by "<<dist
                 <<" but the projected intervals are not separate, "
                 <<"range=["<<r1[0]<<","<<r1[1]<<"] and ["<<r2[0]<<","<<r2[1]<<"]"<<std::endl;
        GJKPolytope<T>::writeConfigVTK("GJKConfig.vtk",_p1,_p2,p1,p2);
        GJKPolytope<T>::writeConfig("GJKConfig.dat",_p1,_p2);
        //ASSERT(false)
      }*/
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
bool CCBarrierEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,Vec4T* x) {
  TH E2;
  Vec4TH G2;
  Mat4TH H2;
  if(_implicit) {
    if(!initialize(NULL,body)) {
      //std::cout<<"Cannot initialize"<<std::endl;
      return false;
    }
    //optimize
    if(!optimize(_x,E2,G2,H2)) {
      //std::cout<<"Cannot optimize"<<std::endl;
      return false;
    }
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  } else {
    //initialize(NULL,body);
    //compute energy
    if(!energy(_x,E2,NULL,NULL))
      return false;
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  }
  if(x)
    (*x)=_x.template cast<T>();
  //compute gradient and hessian
  if(body && grad) {
    Mat3XT DTG;
    Vec4T G=G2.template cast<T>();
    Mat4T H=H2.template cast<T>();
    bool hessian=grad && grad->_HTheta.size()>0;
    computeDTGH(*body,*grad,_x.template cast<T>(),
                hessian? &G: NULL,hessian? &H: NULL,DNDX);
    if(GTheta) {
      GTheta->setZero(body->nrDOF());
      grad->_info.DTG(*body,mapM(DTG=grad->_DTG),mapV(*GTheta));
    }
    if(HTheta) {
      *HTheta=grad->_HTheta;
      grad->_info.toolB(*body,mapM(DTG=grad->_DTG),[&](int r,int c,T val) {
        (*HTheta)(r,c)+=val;
      });
    }
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::evalLRI(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,
     Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g,Vec4T* x) {
  TH E2;
  Vec4TH G2;
  Mat4TH H2;
  if(_implicit) {
    if(!initialize(NULL,body)) {
      //std::cout<<"Cannot initialize"<<std::endl;
      return false;
    }
    //optimize
    if(!optimize(_x,E2,G2,H2)) {
      std::cout<<"Cannot optimize"<<std::endl;
      return false;
    }
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  } else {
    //initialize(NULL,body);
    //compute energy
    if(!energy(_x,E2,NULL,NULL))
      return false;
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  }
  if(x)
    (*x)=_x.template cast<T>();
  //compute gradient and hessian
  if(body && grad) {
    Vec4T G=G2.template cast<T>();
    Mat4T H=H2.template cast<T>();
    bool hessian=grad && grad->_HTheta.size()>0;
    computeDTGHLRI(*body,*grad,_x.template cast<T>(),
                hessian? &G: NULL,hessian? &H: NULL,DNDX,
                DTG1,DTG2,m,g);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierEnergy<T,PFunc,TH>::evalBackward(const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) {
  TH E2;
  Vec4TH G2;
  Mat4TH H2;
  if(_implicit) {
    if(!initialize(NULL,body)) {
      //std::cout<<"Cannot initialize"<<std::endl;
      return false;
    }
    //optimize
    if(!optimize(_x,E2,G2,H2)) {
      std::cout<<"Cannot optimize"<<std::endl;
      return false;
    }
    if(E2==0) return true;
  } else {
    //initialize();
    //compute energy
    if(!energy(_x,E2,NULL,NULL))
      return false;
    if(E2==0) return true;
  }
  {
    Vec4T G=G2.template cast<T>();
    Mat4T H=H2.template cast<T>();
    //std::vector<MatX3T> HThetaD1,HThetaD2;
    computeHBackward(*body,*grad,_x.template cast<T>(),G,H,HThetaD1,HThetaD2);
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
  std::vector<Mat3X4T> DNDX;
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
      CCBarrierEnergy<T,PFunc,TH> e(pass?p:p2,pass?p2:p,barrier,d0,&info,coef,implicit);
      e.setOutput(output);
      if(!implicit)
        e.initialize(NULL,&body);
      if(!e.eval(&E,&body,&info,&DNDX,&GTheta,&HTheta))
        continue;
      //evaluate again
      info2.reset(body,x+dx*DELTA);
      p2=GJKPolytope<T>(JID,body,info2);
      CCBarrierEnergy<T,PFunc,TH> e2(pass?p:p2,pass?p2:p,barrier,d0,&info2,coef,implicit);
      e2.setOutput(output);
      if(!implicit)
        e2.initialize(e._x.template cast<T>());
      if(!e2.eval(&E2,&body,&info2,&DNDX,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dE"+std::string(implicit?"Implicit":"Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG"+std::string(implicit?"Implicit":"Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
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
  std::vector<Mat3X4T> DNDX;
  std::vector<MatX3T> HThetaD1,HThetaD2;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  GEOMETRY_SCALAR DELTA_GEOMETRY=Epsilon<GEOMETRY_SCALAR>::finiteDifferenceEps();
  std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(JID)._mesh);
  std::shared_ptr<MeshExact> local2=std::dynamic_pointer_cast<MeshExact>(body.joint(JID2)._mesh);
  while(true) {
    coef=rand()/(T)RAND_MAX;
    Vec x=Vec::Random(body.nrDOF());
    Vec dx=Vec::Random(body.nrDOF());
    //evaluate energy/gradient
    Vec4T u,u2;
    u.setZero();
    u2.setZero();
    info.reset(body,x);
    p=GJKPolytope<T>(JID,body,info);
    p2=GJKPolytope<T>(JID2,body,info);
    //p.writeVTK("poly1.vtk",&body);
    //p2.writeVTK("poly2.vtk",&body);
    CCBarrierEnergy<T,PFunc,TH> e(p,p2,barrier,d0,&info,coef,implicit);
    e.setOutput(output);
    if(!implicit)
      e.initialize(NULL,&body);
    if(!e.eval(&E,&body,&info,&DNDX,&GTheta,&HTheta))
      continue;
    if(E==0)
      continue;
    //evaluate forward/backward
    {
      e.evalBackward(&body,&info,&HThetaD1,&HThetaD2);
      //evaluate again
      info2.reset(body,x+dx*DELTA);
      p=GJKPolytope<T>(JID,body,info2);
      p2=GJKPolytope<T>(JID2,body,info2);
      CCBarrierEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info2,coef,implicit);
      e2.setOutput(output);
      if(!implicit)
        e2.initialize(e._x.template cast<T>());
      if(!e2.eval(&E2,&body,&info2,&DNDX,&GTheta2,NULL,&u2))
        continue;
      DEBUG_GRADIENT("dE"+std::string(implicit? "Implicit": "Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG"+std::string(implicit? "Implicit": "Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    }
    //evaluate X
    p2=GJKPolytope<T>(JID2,body,info);
    for(int i=0; i<p.globalVss().cols(); i++) {
      //save vertex position and perturb
      dx.setRandom(3);
      MeshExact::Vec3T tmp=local->vss()[i];
      local->vssNonConst()[i]+=dx.template cast<GEOMETRY_SCALAR>()*DELTA_GEOMETRY;
      //debug backward derivative
      info.reset(body,x);
      p=GJKPolytope<T>(JID,body,info);
      CCBarrierEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info,coef,implicit);
      e2.setOutput(output);
      if(!implicit)
        e2.initialize(e._x.template cast<T>());
      e2.eval(&E2,&body,&info,&DNDX,&GTheta2,NULL);
      DEBUG_GRADIENT("dGX-P"+std::string(implicit? "Implicit": "Explicit"),(HThetaD1.at(i)*dx).norm(),(HThetaD1.at(i)*dx-(GTheta2-GTheta)/DELTA_GEOMETRY).norm())
      //load vertex position
      local->vssNonConst()[i]=tmp;
    }
    p=GJKPolytope<T>(JID,body,info);
    for(int i=0; i<p2.globalVss().cols(); i++) {
      //save vertex position and perturb
      dx.setRandom(3);
      MeshExact::Vec3T tmp=local2->vss()[i];
      local2->vssNonConst()[i]+=dx.template cast<GEOMETRY_SCALAR>()*DELTA_GEOMETRY;
      //debug backward derivative
      info.reset(body,x);
      p2=GJKPolytope<T>(JID2,body,info);
      CCBarrierEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info,coef,implicit);
      e2.setOutput(output);
      if(!implicit)
        e2.initialize(e._x.template cast<T>());
      e2.eval(&E2,&body,&info,&DNDX,&GTheta2,NULL);
      DEBUG_GRADIENT("dGX-N"+std::string(implicit? "Implicit": "Explicit"),(HThetaD2.at(i)*dx).norm(),(HThetaD2.at(i)*dx-(GTheta2-GTheta)/DELTA_GEOMETRY).norm())
      //load vertex position
      local2->vssNonConst()[i]=tmp;
    }
    p2=GJKPolytope<T>(JID2,body,info);
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
  for(int c=0; c<_p1.globalVss().cols(); c++) {
    //debugEnergyP(_p1.globalVss().col(c).template cast<TH>(),x,1e-9);
    if(!energyP(_p1.globalVss().col(c).template cast<TH>(),x,E,G,H))
      return false;
  }
  //negative
  for(int c=0; c<_p2.globalVss().cols(); c++) {
    //debugEnergyN(_p2.globalVss().col(c).template cast<TH>(),x,1e-9);
    if(!energyN(_p2.globalVss().col(c).template cast<TH>(),x,E,G,H))
      return false;
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
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,GAll* g) const {
  Mat3T nntDD,C=cross<T>(Rvl);
  Mat3X4T DDEDXDTheta;
  T D=0,DD=0,d=v.dot(x.template segment<3>(0))-x[3];
  T e=_p.template eval<T>(d,&D,&DD,_d0Half,1);
  if(e==0)
    return;
  if(DTG)
    parallelAdd<T,3,4>(*DTG,0,0,computeDTG<T>(x.template segment<3>(0)*D*_coef,vl));
  if(g){
    g->_g1._w+=C*(x.template segment<3>(0)*D*_coef);
    g->_g1._t+=(x.template segment<3>(0)*D*_coef);
  }
  if(LRH && wLRH) {
    ROT(DDEDXDTheta)=(x.template segment<3>(0)*v.transpose())*DD+Mat3T::Identity()*D;
    CTR(DDEDXDTheta)=-x.template segment<3>(0)*DD;
    parallelAdd<T,3,4>(*LRH,0,0,DDEDXDTheta);
    parallelAdd<T,3,4>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(Mww && Mtw && Mwt && Mtt) {
    nntDD=x.template segment<3>(0)*x.template segment<3>(0).transpose()*DD*_coef;
    parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*C.transpose());
    parallelAdd<T,3,3>(*Mtw,0,0,nntDD*C.transpose());
    parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
    parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::energyNDTGH(const Vec3T& v,const Vec3T& vl,const Vec4T& x,Mat3X4T* DTG,
    const Vec3T& Rvl,Mat3X4T* LRH,Mat3X4T* wLRH,
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,GAll* g) const {
  Mat3T nntDD,C=cross<T>(Rvl);
  Mat3X4T DDEDXDTheta;
  T D=0,DD=0,d=x[3]-v.dot(x.template segment<3>(0));
  T e=_p.template eval<T>(d,&D,&DD,_d0Half,1);
  if(e==0)
    return;
  if(DTG)
    parallelAdd<T,3,4>(*DTG,0,0,computeDTG<T>(-x.template segment<3>(0)*D*_coef,vl));
  if(g){
    g->_g2._w+=C*(-x.template segment<3>(0)*D*_coef);
    g->_g2._t+=(-x.template segment<3>(0)*D*_coef);
  }
  if(LRH && wLRH) {
    ROT(DDEDXDTheta)=(x.template segment<3>(0)*v.transpose())*DD-Mat3T::Identity()*D;
    CTR(DDEDXDTheta)=-x.template segment<3>(0)*DD;
    parallelAdd<T,3,4>(*LRH,0,0,DDEDXDTheta);
    parallelAdd<T,3,4>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(Mww && Mtw && Mwt && Mtt) {
    nntDD=x.template segment<3>(0)*x.template segment<3>(0).transpose()*DD*_coef;
    parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*C.transpose());
    parallelAdd<T,3,3>(*Mtw,0,0,nntDD*C.transpose());
    parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
    parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::computeDTGH(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec4T& x,const Vec4T* G,const Mat4T* H,std::vector<Mat3X4T>* DNDX) const {
  Mat3X4T wDDEDXDTheta1,tDDEDXDTheta1,DTG;
  Mat3X4T wDDEDXDTheta2,tDDEDXDTheta2;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat4T S,Hxx;
  if(DNDX) DNDX->resize(4);
  if(G && H && _implicit) {
    sensitivity(x,*G,*H,S);
    Hxx=-S*_coef;
  }
  //positive
  if(_p1.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta1.setZero();
      tDDEDXDTheta1.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p1.jid());
    }
    //stage 1
    DTG.setZero();
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p1.jid())._mesh);
    for(int c=0; c<_p1.globalVss().cols(); c++)
      energyPDTGH(_p1.globalVss().col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta1: NULL,&wDDEDXDTheta1,
                  (G && H)? &Mww: NULL,&Mtw,&Mwt,&Mtt);
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        parallelAdd(info._DTG(r,c+_p1.jid()*4),DTG(r,c));
    //stage 2
    if(G && H) {
      if(_implicit) {
        Mww+=wDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        Mtw+=tDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        Mwt+=wDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        Mtt+=tDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        if(DNDX){
          DNDX->at(0)=wDDEDXDTheta1*Hxx;
          DNDX->at(1)=tDDEDXDTheta1*Hxx;
        } 
      }
      contractHTheta(_p1.jid(),_p1.jid(),body,info,Mww,Mtw,Mwt,Mtt);
    }
  }
  //negative
  if(_p2.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta2.setZero();
      tDDEDXDTheta2.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p2.jid());
    }
    //stage 1
    DTG.setZero();
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_p2.globalVss().cols(); c++)
      energyNDTGH(_p2.globalVss().col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta2: NULL,&wDDEDXDTheta2,
                  (G && H)? &Mww: NULL,&Mtw,&Mwt,&Mtt);
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        parallelAdd(info._DTG(r,c+_p2.jid()*4),DTG(r,c));
    //stage 2
    if(G && H) {
      if(_implicit) {
        Mww+=wDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        Mtw+=tDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        Mwt+=wDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        Mtt+=tDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        if(DNDX){
          DNDX->at(2)=wDDEDXDTheta2*Hxx;
          DNDX->at(3)=tDDEDXDTheta2*Hxx;
        } 
      }
      contractHTheta(_p2.jid(),_p2.jid(),body,info,Mww,Mtw,Mwt,Mtt);
    }
  }
  //positive/negative
  if(_p1.jid()>=0 && _p2.jid()>=0 && _implicit) {
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
void CCBarrierEnergy<T,PFunc,TH>::computeHBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec4T& x,const Vec4T& G,const Mat4T& H,
    std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) const {
  MatX4T HThetaX=MatX4T::Zero(info._info._xM.size(),4);
  Mat3X4T wDDEDXDTheta1,tDDEDXDTheta1,Mwx,DTG;
  Mat3X4T wDDEDXDTheta2,tDDEDXDTheta2,Mtx;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat4T S,Hxx;
  if(_implicit) {
    sensitivity(x,G,H,S);
    Hxx=-S*_coef;
  }
  //positive
  if(_p1.jid()>=0 && _implicit) {
    wDDEDXDTheta1.setZero();
    tDDEDXDTheta1.setZero();
    R=ROTI(info._info._TM,_p1.jid());
    //stage 1
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p1.jid())._mesh);
    for(int c=0; c<_p1.globalVss().cols(); c++)
      energyPDTGH(_p1.globalVss().col(c),local->vss()[c].template cast<T>(),x,NULL,
                  R*local->vss()[c].template cast<T>(),&tDDEDXDTheta1,&wDDEDXDTheta1,
                  NULL,NULL,NULL,NULL);
    //stage 2
    Mwx=wDDEDXDTheta1*Hxx;
    Mtx=tDDEDXDTheta1*Hxx;
    contractHBackward(_p1.jid(),body,info,HThetaX,Mwx,Mtx);
  }
  //negative
  if(_p2.jid()>=0 && _implicit) {
    wDDEDXDTheta2.setZero();
    tDDEDXDTheta2.setZero();
    R=ROTI(info._info._TM,_p2.jid());
    //stage 1
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_p2.globalVss().cols(); c++)
      energyNDTGH(_p2.globalVss().col(c),local->vss()[c].template cast<T>(),x,NULL,
                  R*local->vss()[c].template cast<T>(),&tDDEDXDTheta2,&wDDEDXDTheta2,
                  NULL,NULL,NULL,NULL);
    //stage 2
    Mwx=wDDEDXDTheta2*Hxx;
    Mtx=tDDEDXDTheta2*Hxx;
    contractHBackward(_p2.jid(),body,info,HThetaX,Mwx,Mtx);
  }
  
  //positive
  if(_p1.jid()>=0 && HThetaD1) {
    R=ROTI(info._info._TM,_p1.jid());
    HThetaD1->resize(_p1.globalVss().cols());
    //stage 3
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p1.jid())._mesh);
    for(int c=0; c<_p1.globalVss().cols(); c++) {
      DTG.setZero();
      wDDEDXDTheta1.setZero();
      tDDEDXDTheta1.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      energyPDTGH(_p1.globalVss().col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),&tDDEDXDTheta1,&wDDEDXDTheta1,
                  &Mww,&Mtw,&Mwt,&Mtt);
      Mwt=Mwt*R-cross<T>(DTG.col(3))*R;
      Mtt=Mtt*R;
      HThetaD1->at(c)=HThetaX*tDDEDXDTheta1.transpose()*R;
      contractHBackward(_p1.jid(),body,info,HThetaD1->at(c),Mwt,Mtt);
      for(int i=0; i<body.nrDOF(); i++)
        for(int j=0; j<3; j++)
          parallelAdd(info._HThetaD(i,j+(_p1.getVertexId()[0]+c)*3),HThetaD1->at(c)(i,j));
    }
  }
  //negative
  if(_p2.jid()>=0 && HThetaD2) {
    R=ROTI(info._info._TM,_p2.jid());
    HThetaD2->resize(_p2.globalVss().cols());
    //stage 3
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_p2.globalVss().cols(); c++) {
      DTG.setZero();
      wDDEDXDTheta2.setZero();
      tDDEDXDTheta2.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      energyNDTGH(_p2.globalVss().col(c),local->vss()[c].template cast<T>(),x,&DTG,
                  R*local->vss()[c].template cast<T>(),&tDDEDXDTheta2,&wDDEDXDTheta2,
                  &Mww,&Mtw,&Mwt,&Mtt);
      Mwt=Mwt*R-cross<T>(DTG.col(3))*R;
      Mtt=Mtt*R;
      HThetaD2->at(c)=HThetaX*tDDEDXDTheta2.transpose()*R;
      contractHBackward(_p2.jid(),body,info,HThetaD2->at(c),Mwt,Mtt);
      for(int i=0; i<body.nrDOF(); i++)
        for(int j=0; j<3; j++)
          parallelAdd(info._HThetaD(i,j+(_p2.getVertexId()[0]+c)*3),HThetaD2->at(c)(i,j));
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::computeDTGHLRI(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec4T& x,const Vec4T* G,const Mat4T* H,std::vector<Mat3X4T>* DNDX,
    Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g) const {
  Mat3X4T wDDEDXDTheta1,tDDEDXDTheta1;
  Mat3X4T wDDEDXDTheta2,tDDEDXDTheta2;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat4T S,Hxx;
  if(DNDX) DNDX->resize(4);
  if(G && H && _implicit) {
    sensitivity(x,*G,*H,S);
    Hxx=-S*_coef;
  }
  //positive
  if(_p1.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta1.setZero();
      tDDEDXDTheta1.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p1.jid());
      Mat3X4T trans1=TRANSI(info._info._TM,_p1.jid());
    }
    //stage 1
    DTG1->setZero();
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(_p1.mesh());
    for(int c=0; c<_p1.globalVss().cols(); c++)
      energyPDTGH(_p1.globalVss().col(c),local->vss()[c].template cast<T>(),x,DTG1,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta1: NULL,&wDDEDXDTheta1,
                  (G && H)? &(m._m11._Mww): NULL,&(m._m11._Mtw),&(m._m11._Mwt),&(m._m11._Mtt),&g);
    //stage 2
    if(G && H) {
      if(_implicit) {
        m._m11._Mww+=wDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        m._m11._Mtw+=tDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        m._m11._Mwt+=wDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        m._m11._Mtt+=tDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        if(DNDX){
          DNDX->at(0)=wDDEDXDTheta1*Hxx;
          DNDX->at(1)=tDDEDXDTheta1*Hxx;
        } 
      }
    }
  }
  //negative
  if(_p2.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta2.setZero();
      tDDEDXDTheta2.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p2.jid());
    }
    //stage 1
    DTG2->setZero();
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(_p2.mesh());
    for(int c=0; c<_p2.globalVss().cols(); c++)
      energyNDTGH(_p2.globalVss().col(c),local->vss()[c].template cast<T>(),x,DTG2,
                  R*local->vss()[c].template cast<T>(),(G && H && _implicit)? &tDDEDXDTheta2: NULL,&wDDEDXDTheta2,
                  (G && H)? &(m._m22._Mww): NULL,&(m._m22._Mtw),&(m._m22._Mwt),&(m._m22._Mtt),&g);
    //stage 2
    if(G && H) {
      if(_implicit) {
        m._m22._Mww+=wDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        m._m22._Mtw+=tDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        m._m22._Mwt+=wDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        m._m22._Mtt+=tDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        if(DNDX){
          DNDX->at(2)=wDDEDXDTheta2*Hxx;
          DNDX->at(3)=tDDEDXDTheta2*Hxx;
        } 
      }
    }
  }
  //positive/negative
  if(_p1.jid()>=0 && _p2.jid()>=0 && _implicit) {
    //stage 2
    if(G && H) {
      m._m12._Mww=wDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose();
      m._m12._Mtw=tDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose();
      m._m12._Mwt=wDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose();
      m._m12._Mtt=tDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose();
      m._m21._Mww=(wDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose()).transpose();
      m._m21._Mwt=(tDDEDXDTheta1*Hxx*wDDEDXDTheta2.transpose()).transpose();
      m._m21._Mtw=(wDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose()).transpose();
      m._m21._Mtt=(tDDEDXDTheta1*Hxx*tDDEDXDTheta2.transpose()).transpose();
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::contractHTheta
(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,
 const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const {
  info._info.JRCSparse(body,kL,[&](int row,const Vec3T& JRL) {
    info._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      parallelAdd(info._HTheta(row,col),JRL.dot(Mww*JRR));
    },[&](int col,const Vec3T& JCR) {
      parallelAdd(info._HTheta(row,col),JRL.dot(Mwt*JCR));
    });
  },[&](int row,const Vec3T& JCL) {
    info._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      parallelAdd(info._HTheta(row,col),JCL.dot(Mtw*JRR));
    },[&](int col,const Vec3T& JCR) {
      parallelAdd(info._HTheta(row,col),JCL.dot(Mtt*JCR));
    });
  });
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::contractHThetaL
(int kL,int kR,const ArticulatedBody& body,CollisionGradInfo<T>& info,CollisionGradInfo<T>& infoL,
 const Mat3T& Mww,const Mat3T& Mtw,const Mat3T& Mwt,const Mat3T& Mtt) const {
  info._info.JRCSparse(body,kL,[&](int row,const Vec3T& JRL) {
    infoL._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      parallelAdd(info._HThetaL(row,col),JRL.dot(Mww*JRR));
    },[&](int col,const Vec3T& JCR) {
      parallelAdd(info._HThetaL(row,col),JRL.dot(Mwt*JCR));
    });
  },[&](int row,const Vec3T& JCL) {
    infoL._info.JRCSparse(body,kR,[&](int col,const Vec3T& JRR) {
      parallelAdd(info._HThetaL(row,col),JCL.dot(Mtw*JRR));
    },[&](int col,const Vec3T& JCR) {
      parallelAdd(info._HThetaL(row,col),JCL.dot(Mtt*JCR));
    });
  });
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::contractHBackward
(int k,const ArticulatedBody& body,const CollisionGradInfo<T>& info,
 MatX4T& HThetaX,const Mat3X4T& Mwx,const Mat3X4T& Mtx) const {
  info._info.JRCSparse(body,k,[&](int row,const Vec3T& JR) {
    parallelAdd<T,4>(HThetaX,row,0,JR.transpose()*Mwx);
  },[&](int row,const Vec3T& JC) {
    parallelAdd<T,4>(HThetaX,row,0,JC.transpose()*Mtx);
  });
}
template <typename T,typename PFunc,typename TH>
void CCBarrierEnergy<T,PFunc,TH>::contractHBackward
(int k,const ArticulatedBody& body,const CollisionGradInfo<T>& info,
 MatX3T& HThetaX,const Mat3T& Mwt,const Mat3T& Mtt) const {
  info._info.JRCSparse(body,k,[&](int row,const Vec3T& JR) {
    parallelAdd<T,3>(HThetaX,row,0,JR.transpose()*Mwt);
  },[&](int row,const Vec3T& JC) {
    parallelAdd<T,3>(HThetaX,row,0,JC.transpose()*Mtt);
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
