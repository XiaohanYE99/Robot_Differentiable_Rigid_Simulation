#include "ConvexHullDistanceFrictionEnergy.h"
#include <Utils/CrossSpatialUtils.h>
#include <Utils/SparseUtils.h>

namespace PHYSICSMOTION {
//CCBarrierEnergy
template <typename T,typename PFunc,typename TH>
CCBarrierFrictionEnergy<T,PFunc,TH>::CCBarrierFrictionEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const GJKPolytope<T>& pl1,const GJKPolytope<T>& pl2,const Vec4T& x,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,T dt,bool implicit)
  :CCBarrierEnergy<T,PFunc,TH>(p1,p2,p,d0,grad,coef,implicit),_pl1(pl1),_pl2(pl2),_dt(dt),_eps(1e-6),_fri(0.8) {
  _x=Vec4TH((TH)x[0],(TH)x[1],(TH)x[2],(TH)x[3]);
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::initialize(Vec4T* res,const ArticulatedBody* body) {
  FUNCTION_NOT_IMPLEMENTED
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::debugGradient(bool implicit,const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0,bool output) {
  T E,E2,ETmp,coef,dt=0.01f;
  PFunc barrier;
  Vec GTheta,GTheta2;
  MatT HTheta,HThetaL;
  GJKPolytope<T> p2,pl2;
  const GJKPolytope<T>& pl=p;
  CollisionGradInfo<T> info,info2,infoL;
  std::vector<Mat3X4T> DNDX;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  for(int pass=0; pass<2; pass++)
    while(true) {
      coef=rand()/(T)RAND_MAX;
      Vec x=Vec::Random(body.nrDOF());
      Vec xL=x+dt*Vec::Random(body.nrDOF());
      Vec dx=Vec::Random(body.nrDOF());
      //evaluate energy/gradient
      info.reset(body,x);
      infoL.reset(body,xL);
      p2=GJKPolytope<T>(JID,body,info);
      pl2=GJKPolytope<T>(JID,body,infoL);
      //calculate separating plane
      CCBarrierConvexEnergy<T,PFunc,TH> el(pass?pl:pl2,pass?pl2:pl,barrier,d0,&info,coef);
      if(!el.eval(&ETmp,&body,&infoL,&DNDX,NULL,NULL))
        continue;
      if(ETmp==0)
        continue;
      CCBarrierFrictionEnergy<T,PFunc,TH> e(pass?p:p2,pass?p2:p,pass?pl:pl2,pass?pl2:pl,el.getX(),barrier,d0,&info,coef,dt,implicit);
      e.setOutput(output);
      if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
        continue;
      //evaluate backward
      {
        std::vector<MatX3T> HThetaD1,HThetaD2;
        e.evalBackward(&body,&info,&infoL,&DNDX,&HThetaD1,&HThetaD2);
        HThetaL=info._HThetaL;
        //evaluate again
        info.reset(body,x);
        info2.reset(body,xL+dx*DELTA);
        pl2=GJKPolytope<T>(JID,body,info2);
        //calculate separating plane
        CCBarrierConvexEnergy<T,PFunc,TH> el2(pass?pl:pl2,pass?pl2:pl,barrier,d0,&info2,coef);
        el2.eval(&ETmp,&body,&info2,&DNDX,NULL,NULL);
        CCBarrierFrictionEnergy<T,PFunc,TH> e2(pass?p:p2,pass?p2:p,pass?pl:pl2,pass?pl2:pl,el2.getX(),barrier,d0,&info,coef,dt,implicit);
        e2.setOutput(output);
        if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
          continue;
        DEBUG_GRADIENT("dGL"+std::string(implicit?"Implicit":"Explicit"),(HThetaL*dx).norm(),(HThetaL*dx-(GTheta2-GTheta)/DELTA).norm())
      }
      //eval forward
      {
        info.reset(body,x+dx*DELTA);
        infoL.reset(body,xL);
        p2=GJKPolytope<T>(JID,body,info);
        pl2=GJKPolytope<T>(JID,body,infoL);
        //calculate separating plane
        CCBarrierConvexEnergy<T,PFunc,TH> el2(pass?pl:pl2,pass?pl2:pl,barrier,d0,&info2,coef);
        el2.eval(&ETmp,&body,&info2,&DNDX,NULL,NULL);
        CCBarrierFrictionEnergy<T,PFunc,TH> e2(pass?p:p2,pass?p2:p,pass?pl:pl2,pass?pl2:pl,el.getX(),barrier,d0,&info,coef,dt,implicit);
        e2.setOutput(output);
        if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
          continue;
        DEBUG_GRADIENT("dE"+std::string(implicit? "Implicit": "Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
        DEBUG_GRADIENT("dG"+std::string(implicit? "Implicit": "Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
      }
      break;
    }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::debugGradient(bool implicit,const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
  T E,E2,ETmp,coef,dt=0.01f;
  PFunc barrier;
  Vec GTheta,GTheta2;
  MatT HTheta,HThetaL;
  GJKPolytope<T> p,pl,p2,pl2;
  CollisionGradInfo<T> info,info2,infoL;
  std::vector<Mat3X4T> DNDX;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  while(true) {
    coef=rand()/(T)RAND_MAX;
    Vec x=Vec::Random(body.nrDOF());
    Vec xL=x+dt*Vec::Random(body.nrDOF());
    Vec dx=Vec::Random(body.nrDOF());
    Vec3T u,u2;
    u.setZero();
    u2.setZero();
    //evaluate energy/gradient
    info.reset(body,x);
    infoL.reset(body,xL);
    p=GJKPolytope<T>(JID,body,info);
    pl=GJKPolytope<T>(JID,body,infoL);
    p2=GJKPolytope<T>(JID2,body,info);
    pl2=GJKPolytope<T>(JID2,body,infoL);
    //calculate separating plane
    CCBarrierConvexEnergy<T,PFunc,TH> el(pl,pl2,barrier,d0,&infoL,coef);
    if(!el.eval(&ETmp,&body,&infoL,&DNDX,NULL,NULL))
      continue;
    if(ETmp==0)
      continue;
    CCBarrierFrictionEnergy<T,PFunc,TH> e(p,p2,pl,pl2,el.getX(),barrier,d0,&info,coef,dt,implicit);
    e.setOutput(output);
    if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
      continue;
    //evaluate backward
    {
      std::vector<MatX3T> HThetaD1,HThetaD2;
      e.evalBackward(&body,&info,&infoL,&DNDX,&HThetaD1,&HThetaD2);
      HThetaL=info._HThetaL;
      //evaluate again
      info.reset(body,x);
      info2.reset(body,xL+dx*DELTA);
      pl=GJKPolytope<T>(JID,body,info2);
      pl2=GJKPolytope<T>(JID2,body,info2);
      //calculate separating plane
      CCBarrierConvexEnergy<T,PFunc,TH> el2(pl,pl2,barrier,d0,&info2,coef);
      el2.eval(&ETmp,&body,&info2,&DNDX,NULL,NULL);
      CCBarrierFrictionEnergy<T,PFunc,TH> e2(p,p2,pl,pl2,el2.getX(),barrier,d0,&info,coef,dt,implicit);
      e2.setOutput(output);
      if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dGL"+std::string(implicit? "Implicit": "Explicit"),(HThetaL*dx).norm(),(HThetaL*dx-(GTheta2-GTheta)/DELTA).norm())
    }
    //eval forward
    {
      info.reset(body,x+dx*DELTA);
      infoL.reset(body,xL);
      p=GJKPolytope<T>(JID,body,info);
      pl=GJKPolytope<T>(JID,body,infoL);
      p2=GJKPolytope<T>(JID2,body,info);
      pl2=GJKPolytope<T>(JID2,body,infoL);
      //calculate separating plane
      CCBarrierConvexEnergy<T,PFunc,TH> el2(pl,pl2,barrier,d0,&info2,coef);
      el2.eval(&ETmp,&body,&info2,&DNDX,NULL,NULL);
      CCBarrierFrictionEnergy<T,PFunc,TH> e2(p,p2,pl,pl2,el.getX(),barrier,d0,&info,coef,dt,implicit);
      e2.setOutput(output);
      if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dE"+std::string(implicit? "Implicit": "Explicit"),GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG"+std::string(implicit? "Implicit": "Explicit"),(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    }
    break;
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::debugEnergyP(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,T perturbRange) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Vec3TH G,G2,du,u2;
  du.setRandom();
  Mat3TH H,H2;
  // Debug energyP
  E=E2=0;
  G.setZero();
  G2.setZero();
  H.setZero();
  H2.setZero();
  energyP(v,vLast,u,x,E,&G,&H);
  u2=u+du*DELTA;
  energyP(v,vLast,u2,x,E2,&G2,&H2);
  if(E==0 && E2==0)return;
  DEBUG_GRADIENT("energyP E",G.dot(du),G.dot(du)-(E2-E)/DELTA)
  DEBUG_GRADIENT("energyP G",(H*du).norm(),(H*du-(G2-G)/DELTA).norm())
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::debugEnergyN(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,T perturbRange) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Vec3TH G,G2,du,u2;
  du.setRandom();
  Mat3TH H,H2;
  // Debug energyP
  E=E2=0;
  G.setZero();
  G2.setZero();
  H.setZero();
  H2.setZero();
  energyN(v,vLast,u,x,E,&G,&H);
  u2=u+du*DELTA;
  energyN(v,vLast,u2,x,E2,&G2,&H2);
  if(E==0 && E2==0)return;
  DEBUG_GRADIENT("energyP E",G.dot(du),G.dot(du)-(E2-E)/DELTA)
  DEBUG_GRADIENT("energyP G",(H*du).norm(),(H*du-(G2-G)/DELTA).norm())
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::debugEnergy(const Vec3TH& x) const {
  DEFINE_NUMERIC_DELTA_T(TH)
  TH E,E2;
  Mat3TH H;
  Vec3TH G,G2,dx;
  dx.setRandom();
  energy(x,E,&G,&H);
  energy(x+dx*DELTA,E2,&G2,NULL);
  DEBUG_GRADIENT("E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("G",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
template <typename T,typename PFunc,typename TH>
T CCBarrierFrictionEnergy<T,PFunc,TH>::fri() const {
  return _fri;
}
template <typename T,typename PFunc,typename TH>
T& CCBarrierFrictionEnergy<T,PFunc,TH>::fri() {
  return _fri;
}
//helper
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::energyP(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,TH& E,Vec3TH* G,Mat3TH* H) const {
  TH e,D_=0,D=0,DD=0,d=vLast.dot(x.template segment<3>(0))-x[3];
  e=_p.template eval<TH>(d,&D_,NULL,(TH)_d0Half,1);
  if(!isfinite(e))
    return false;
  if(e==0)
    return true;
  Vec3TH f_n=D_*(x.template segment<3>(0));
  Vec3TH x_v=(Mat3TH::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm())*(v-vLast)*TH(1.0/_dt);
  E+=TH(_fri*_dt)*f_n.norm()*sqrt((x_v-u).squaredNorm()+TH(_eps));
  Vec3TH DX=-2*(x_v-u);
  Mat3TH DDX=2*Mat3TH::Identity();
  if(G) {
    D=0.5/sqrt((x_v-u).squaredNorm()+TH(_eps));
    *G+=TH(_fri*_dt)*f_n.norm()*D*DX;
  }
  if(H) {
    DD=-0.25/pow((x_v-u).squaredNorm()+TH(_eps),1.5);
    *H+=TH(_fri*_dt)*f_n.norm()*(DX*DX.transpose()*DD+D*DDX);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::energyN(const Vec3TH& v,const Vec3TH& vLast,const Vec3TH& u,const Vec4TH& x,TH& E,Vec3TH* G,Mat3TH* H) const {
  TH e,D_=0,D=0,DD=0,d=x[3]-vLast.dot(x.template segment<3>(0));
  e=_p.template eval<TH>(d,&D_,NULL,(TH)_d0Half,1);
  if(!isfinite(e))
    return false;
  if(e==0)
    return true;
  Vec3TH f_n=D_*(x.template segment<3>(0));
  Vec3TH x_v=(Mat3TH::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm())*(v-vLast)*TH(1.0/_dt);
  E+=TH(_fri*_dt)*f_n.norm()*sqrt((x_v-u).squaredNorm()+TH(_eps));
  Vec3TH DX=-2*(x_v-u);
  Mat3TH DDX=2*Mat3TH::Identity();
  if(G) {
    D=0.5/sqrt((x_v-u).squaredNorm()+TH(_eps));
    *G+=TH(_fri*_dt)*f_n.norm()*D*DX;
  }
  if(H) {
    DD=-0.25/pow((x_v-u).squaredNorm()+TH(_eps),1.5);
    *H+=TH(_fri*_dt)*f_n.norm()*(DX*DX.transpose()*DD+D*DDX);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::energy(const Vec3TH& u,TH& E,Vec3TH* G,Mat3TH* H) const {
  E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  //positive
  for(int c=0; c<_p1.globalVss().cols(); c++) {
    //debugEnergyP(_p1.globalVss().col(c).template cast<TH>(),_pl1.globalVss().col(c).template cast<TH>(),u,_x,1e-9);
    if(!energyP(_p1.globalVss().col(c).template cast<TH>(),_pl1.globalVss().col(c).template cast<TH>(),u,_x,E,G,H))
      return false;
  }
  //negative
  for(int c=0; c<_p2.globalVss().cols(); c++) {
    //debugEnergyN(_p2.globalVss().col(c).template cast<TH>(),_pl2.globalVss().col(c).template cast<TH>(),u,_x,1e-9);
    if(!energyN(_p2.globalVss().col(c).template cast<TH>(),_pl2.globalVss().col(c).template cast<TH>(),u,_x,E,G,H))
      return false;
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta) {
  TH E2;
  Vec3TH G2;
  Mat3TH H2;
  if(_implicit) {
    //optimize
    _u.setZero();
    if(!optimize(_u,E2,G2,H2)) {
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
    //compute energy
    _u.setZero();
    if(!energy(_u,E2,NULL,NULL))
      return false;
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  }
  if(body && grad) {
    Mat3XT DTG;
    Vec3T G=G2.template cast<T>();
    Mat3T H=H2.template cast<T>();
    bool hessian=grad && grad->_HTheta.size()>0;
    computeDTGH(*body,*grad,_u.template cast<T>(),_x.template cast<T>(),
                hessian?&G:NULL,hessian?&H:NULL);
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
bool CCBarrierFrictionEnergy<T,PFunc,TH>::evalLRI(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta,
               Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g) {
  TH E2;
  Vec3TH G2;
  Mat3TH H2;
  if(_implicit) {
    //optimize
    _u.setZero();
    if(!optimize(_u,E2,G2,H2)) {
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
    //compute energy
    _u.setZero();
    if(!energy(_u,E2,NULL,NULL))
      return false;
    if(E)
      *E=_coef*(T)E2;
    if(E2==0) {
      if(GTheta)
        GTheta->setZero(body->nrDOF());
      return true;
    }
  }
  if(body && grad) {
    Mat3XT DTG;
    Vec3T G=G2.template cast<T>();
    Mat3T H=H2.template cast<T>();
    bool hessian=grad && grad->_HTheta.size()>0;
    computeDTGHLRI(*body,*grad,_u.template cast<T>(),_x.template cast<T>(),
                hessian?&G:NULL,hessian?&H:NULL,DTG1,DTG2,m,g);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::evalBackward(const ArticulatedBody* body,CollisionGradInfo<T>* grad,CollisionGradInfo<T>* Pos,
    std::vector<Mat3X4T>* DNDX,std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) {
  TH E2;
  Vec3TH G2;
  Mat3TH H2;
  if(_implicit) {
    _u.setZero();
    //optimize
    if(!optimize(_u,E2,G2,H2)) {
      std::cout<<"Cannot optimize"<<std::endl;
      return false;
    }
    if(E2==0) return true;
  } else {
    //compute energy
    _u.setZero();
    if(!energy(_u,E2,NULL,NULL))
      return false;
    if(E2==0) return true;
  }
  {
    Vec3T G=G2.template cast<T>();
    Mat3T H=H2.template cast<T>();
    Mat3XT DTG;

    computeHBackward(*body,*grad,_u.template cast<T>(),_x.template cast<T>(),G,H,HThetaD1,HThetaD2);
    computeHLBackward(*body,*grad,*Pos,_u.template cast<T>(),_x.template cast<T>(),&G,&H,DNDX);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierFrictionEnergy<T,PFunc,TH>::optimize(Vec3TH& u,TH& E,Vec3TH& G,Mat3TH& H) const {
  Vec3TH D,G2,u2;
  Eigen::FullPivLU<Mat3TH> invH;
  Mat3TH id=Mat3TH::Identity(),H2;
  TH E2,alpha=1,alphaDec=0.5,alphaInc=3.0,alphaMax=1e10;
  //debugEnergy(u);
  //assemble
  if(!energy(u,E,&G,&H))
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
    u2=u-D;
    if(energy(u2,E2,&G2,&H2) && E2<E) {
      alpha=std::max<TH>(alpha*alphaDec,Epsilon<TH>::defaultEps());
      u=u2;
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
void CCBarrierFrictionEnergy<T,PFunc,TH>::computeDTGHLRI(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec3T& u,const Vec4T& x,const Vec3T* G,const Mat3T* H,
    Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g)const {
  Mat3T wDDEDXDTheta1,tDDEDXDTheta1;
  Mat3X4T DTG;
  Mat3T wDDEDXDTheta2,tDDEDXDTheta2;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat3T S,Hxx;
  if(G && H && _implicit) {
    sensitivity(u,*G,*H,S);
    Hxx=-_coef*S;
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
    DTG1->setZero();
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(_p1.mesh());
    for(int c=0; c<_p1.globalVss().cols(); c++)
      energyPDTGH(_p1.globalVss().col(c),_pl1.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,DTG1,
                  R*local->vss()[c].template cast<T>(),NULL,(G && H && _implicit)?&tDDEDXDTheta1:NULL,&wDDEDXDTheta1,NULL,NULL,NULL,NULL,
                  (G && H)? &(m._m11._Mww): NULL,&(m._m11._Mtw),&(m._m11._Mwt),&(m._m11._Mtt),NULL,&g);
    //stage 2
    if(G && H) {
      if(_implicit) {
        m._m11._Mww+=wDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        m._m11._Mtw+=tDDEDXDTheta1*(Hxx*wDDEDXDTheta1.transpose());
        m._m11._Mwt+=wDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
        m._m11._Mtt+=tDDEDXDTheta1*(Hxx*tDDEDXDTheta1.transpose());
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
      energyNDTGH(_p2.globalVss().col(c),_pl2.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,DTG2,
                  R*local->vss()[c].template cast<T>(),NULL,(G && H && _implicit)?&tDDEDXDTheta2:NULL,&wDDEDXDTheta2,NULL,NULL,NULL,NULL,
                  (G && H)? &(m._m22._Mww): NULL,&(m._m22._Mtw),&(m._m22._Mwt),&(m._m22._Mtt),NULL,&g);
    //stage 2
    if(G && H) {
      if(_implicit) {
        m._m22._Mww+=wDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        m._m22._Mtw+=tDDEDXDTheta2*(Hxx*wDDEDXDTheta2.transpose());
        m._m22._Mwt+=wDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
        m._m22._Mtt+=tDDEDXDTheta2*(Hxx*tDDEDXDTheta2.transpose());
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
void CCBarrierFrictionEnergy<T,PFunc,TH>::computeDTGH(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec3T& u,const Vec4T& x,const Vec3T* G,const Mat3T* H) const {
  Mat3T wDDEDXDTheta1,tDDEDXDTheta1;
  Mat3X4T DTG;
  Mat3T wDDEDXDTheta2,tDDEDXDTheta2;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat3T S,Hxx;
  if(G && H && _implicit) {
    sensitivity(u,*G,*H,S);
    Hxx=-_coef*S;
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
      energyPDTGH(_p1.globalVss().col(c),_pl1.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,&DTG,
                  R*local->vss()[c].template cast<T>(),NULL,(G && H && _implicit)?&tDDEDXDTheta1:NULL,&wDDEDXDTheta1,NULL,NULL,NULL,NULL,
                  (G && H)?&Mww:NULL,&Mtw,&Mwt,&Mtt);
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
      }
      CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p1.jid(),_p1.jid(),body,info,Mww,Mtw,Mwt,Mtt);
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
      energyNDTGH(_p2.globalVss().col(c),_pl2.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,&DTG,
                  R*local->vss()[c].template cast<T>(),NULL,(G && H && _implicit)?&tDDEDXDTheta2:NULL,&wDDEDXDTheta2,NULL,NULL,NULL,NULL,
                  (G && H)?&Mww:NULL,&Mtw,&Mwt,&Mtt);
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
      }
      CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p2.jid(),_p2.jid(),body,info,Mww,Mtw,Mwt,Mtt);
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
      CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p1.jid(),_p2.jid(),body,info,Mww,Mtw,Mwt,Mtt);
      CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p2.jid(),_p1.jid(),body,info,Mww.transpose(),Mwt.transpose(),Mtw.transpose(),Mtt.transpose());
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::computeHLBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,CollisionGradInfo<T>& infoL,
    const Vec3T& u,const Vec4T& x,const Vec3T* G,const Mat3T* H,const std::vector<Mat3X4T>* DNDX) const {
  Mat3T wDDEDXDTheta1,tDDEDXDTheta1,wDDEDXDThetaL1,tDDEDXDThetaL1;
  Mat3T wDDEDXDTheta2,tDDEDXDTheta2,wDDEDXDThetaL2,tDDEDXDThetaL2;
  Mat3X4T wDDEDXDThetaN1,tDDEDXDThetaN1,wDDEDXDThetaN2,tDDEDXDThetaN2,DDEDXDThetaN;
  Mat3T Mwu,Mtu;
  MatX3T HThetaU2=MatX3T::Zero(info._info._xM.size(),3);
  MatX3T HThetaUL2=MatX3T::Zero(info._info._xM.size(),3);
  MatX4T HThetaN2=MatX4T::Zero(info._info._xM.size(),4);
  MatX4T HThetaNL2=MatX4T::Zero(info._info._xM.size(),4);
  MatX3T HThetaU1=MatX3T::Zero(info._info._xM.size(),3);
  MatX3T HThetaUL1=MatX3T::Zero(info._info._xM.size(),3);
  MatX4T HThetaN1=MatX4T::Zero(info._info._xM.size(),4);
  MatX4T HThetaNL1=MatX4T::Zero(info._info._xM.size(),4);
  Mat3T Mww,Mtw,Mwt,Mtt,R,RL;
  Mat3T S,Hxx;
  MatT tmp,tmpi;
  DDEDXDThetaN.setZero();
  if(G && H && _implicit) {
    sensitivity(u,*G,*H,S);
    Hxx=-S*_coef;
  }
  //cal DUDN
  for(int c=0; c<_p1.globalVss().cols(); c++) energyPDUDN(_p1.globalVss().col(c),_pl1.globalVss().col(c),u,x,&DDEDXDThetaN);
  for(int c=0; c<_p2.globalVss().cols(); c++) energyNDUDN(_p2.globalVss().col(c),_pl2.globalVss().col(c),u,x,&DDEDXDThetaN);
  //positive
  if(_p1.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta1.setZero();
      tDDEDXDTheta1.setZero();
      wDDEDXDThetaL1.setZero();
      tDDEDXDThetaL1.setZero();
      wDDEDXDThetaN1.setZero();
      tDDEDXDThetaN1.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p1.jid());
      RL=ROTI(infoL._info._TM,_pl1.jid());
    }
    //stage 1
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p1.jid())._mesh);
    for(int c=0; c<_p1.globalVss().cols(); c++) {
      Vec3T Rvll=RL*local->vss()[c].template cast<T>();
      energyPDTGH(_p1.globalVss().col(c),_pl1.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,NULL,
                  R*local->vss()[c].template cast<T>(),&Rvll,&tDDEDXDTheta1,&wDDEDXDTheta1,&tDDEDXDThetaL1,&wDDEDXDThetaL1,&tDDEDXDThetaN1,&wDDEDXDThetaN1,
                  &Mww,&Mtw,&Mwt,&Mtt,&Hxx);
    }
    //stage 2
    if(G && H) {
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,infoL,HThetaNL1,DNDX->at(0),DNDX->at(1));
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,info,HThetaN1,wDDEDXDThetaN1,tDDEDXDThetaN1);
      CCBarrierEnergy<T,PFunc,TH>::contractHThetaL(_p1.jid(),_p1.jid(),body,info,infoL,Mww,Mtw,Mwt,Mtt);
      tmp=HThetaN1*HThetaNL1.transpose();
      parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmp);
      //std::cout<<DNDX->at(0).norm()<<" "<<DNDX->at(1).norm()<<" "<<_pl1.jid()<<std::endl;
      if(_implicit) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,info,HThetaU1,wDDEDXDTheta1,tDDEDXDTheta1);
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,infoL,HThetaUL1,wDDEDXDThetaL1.transpose()*Hxx,tDDEDXDThetaL1.transpose()*Hxx);
        tmpi=HThetaU1*(HThetaUL1.transpose()+Hxx*DDEDXDThetaN*HThetaNL1.transpose()/_coef);
        parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmpi);
      }
    }
  }
  //negative
  if(_p2.jid()>=0) {
    if(G && H) {
      wDDEDXDTheta2.setZero();
      tDDEDXDTheta2.setZero();
      wDDEDXDThetaL2.setZero();
      tDDEDXDThetaL2.setZero();
      wDDEDXDThetaN2.setZero();
      tDDEDXDThetaN2.setZero();
      Mww=Mtw=Mwt=Mtt=Mat3T::Zero();
      R=ROTI(info._info._TM,_p2.jid());
      RL=ROTI(infoL._info._TM,_pl2.jid());
    }
    //stage 1
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_p2.globalVss().cols(); c++) {
      Vec3T Rvll=RL*local->vss()[c].template cast<T>();
      energyNDTGH(_p2.globalVss().col(c),_pl2.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,NULL,
                  R*local->vss()[c].template cast<T>(),&Rvll,&tDDEDXDTheta2,&wDDEDXDTheta2,&tDDEDXDThetaL2,&wDDEDXDThetaL2,&tDDEDXDThetaN2,&wDDEDXDThetaN2,
                  &Mww,&Mtw,&Mwt,&Mtt,&Hxx);
    }
    //stage 2
    if(G && H) {
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,infoL,HThetaNL2,DNDX->at(2),DNDX->at(3));
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,info,HThetaN2,wDDEDXDThetaN2,tDDEDXDThetaN2);
      CCBarrierEnergy<T,PFunc,TH>::contractHThetaL(_p2.jid(),_p2.jid(),body,info,infoL,Mww,Mtw,Mwt,Mtt);
      tmp=HThetaN2*HThetaNL2.transpose();
      parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmp);
      //std::cout<<_pl2.wDNDX().norm()<<" "<<_pl2.tDNDX().norm()<<" "<<_pl2.jid()<<std::endl;
      if(_implicit) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,info,HThetaU2,wDDEDXDTheta2,tDDEDXDTheta2);
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,infoL,HThetaUL2,wDDEDXDThetaL2.transpose()*Hxx,tDDEDXDThetaL2.transpose()*Hxx);
        tmpi=HThetaU2*(HThetaUL2.transpose()+Hxx*DDEDXDThetaN*HThetaNL2.transpose()/_coef);
        parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmpi);
      }
    }
  }
  //positive/negative
  if(_p1.jid()>=0 && _p2.jid()>=0) {
    //stage 2
    if(G && H) {
      if(_implicit) {
        tmpi=HThetaU1*(HThetaUL2.transpose()+Hxx*DDEDXDThetaN*HThetaNL2.transpose()/_coef)+HThetaU2*(HThetaUL1.transpose()+Hxx*DDEDXDThetaN*HThetaNL1.transpose()/_coef);
        parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmpi);
      }
      tmp=HThetaN1*HThetaNL2.transpose()+HThetaN2*HThetaNL1.transpose();
      parallelAdd<T,-1,-1>(info._HThetaL,0,0,tmp);
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::computeHBackward(const ArticulatedBody& body,CollisionGradInfo<T>& info,
    const Vec3T& u,const Vec4T& x,const Vec3T& G,const Mat3T& H,
    std::vector<MatX3T>* HThetaD1,std::vector<MatX3T>* HThetaD2) const {
  MatX3T HThetaX=MatX3T::Zero(info._info._xM.size(),3);
  Mat3T wDDEDXDTheta1,tDDEDXDTheta1,Mwx;
  Mat3X4T DTG;
  Mat3T wDDEDXDTheta2,tDDEDXDTheta2,Mtx;
  Mat3T Mww,Mtw,Mwt,Mtt,R;
  Mat3T S,Hxx;
  if(_implicit) {
    sensitivity(u,G,H,S);
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
      energyPDTGH(_p1.globalVss().col(c),_pl1.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,NULL,
                  R*local->vss()[c].template cast<T>(),NULL,&tDDEDXDTheta1,&wDDEDXDTheta1,NULL,NULL,NULL,NULL,
                  NULL,NULL,NULL,NULL);
    //stage 2
    Mwx=wDDEDXDTheta1*Hxx;
    Mtx=tDDEDXDTheta1*Hxx;
    CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,info,HThetaX,Mwx,Mtx);
  }
  //negative
  if(_p2.jid()>=0 && _implicit) {
    wDDEDXDTheta2.setZero();
    tDDEDXDTheta2.setZero();
    R=ROTI(info._info._TM,_p2.jid());
    //stage 1
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(_p2.jid())._mesh);
    for(int c=0; c<_p2.globalVss().cols(); c++)
      energyNDTGH(_p2.globalVss().col(c),_pl2.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,NULL,
                  R*local->vss()[c].template cast<T>(),NULL,&tDDEDXDTheta2,&wDDEDXDTheta2,NULL,NULL,NULL,NULL,
                  NULL,NULL,NULL,NULL);
    //stage 2
    Mwx=wDDEDXDTheta2*Hxx;
    Mtx=tDDEDXDTheta2*Hxx;
    CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,info,HThetaX,Mwx,Mtx);
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
      energyPDTGH(_p1.globalVss().col(c),_pl1.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,&DTG,
                  R*local->vss()[c].template cast<T>(),NULL,&tDDEDXDTheta1,&wDDEDXDTheta1,NULL,NULL,NULL,NULL,
                  &Mww,&Mtw,&Mwt,&Mtt);
      Mwt=Mwt*R-cross<T>(DTG.col(3))*R;
      Mtt=Mtt*R;
      HThetaD1->at(c)=HThetaX*tDDEDXDTheta1.transpose()*R;
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,info,HThetaD1->at(c),Mwt,Mtt);
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
      energyNDTGH(_p2.globalVss().col(c),_pl2.globalVss().col(c),local->vss()[c].template cast<T>(),u,x,&DTG,
                  R*local->vss()[c].template cast<T>(),NULL,&tDDEDXDTheta2,&wDDEDXDTheta2,NULL,NULL,NULL,NULL,
                  &Mww,&Mtw,&Mwt,&Mtt);
      Mwt=Mwt*R-cross<T>(DTG.col(3))*R;
      Mtt=Mtt*R;
      HThetaD2->at(c)=HThetaX*tDDEDXDTheta2.transpose()*R;
      CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,info,HThetaD2->at(c),Mwt,Mtt);
      for(int i=0; i<body.nrDOF(); i++)
        for(int j=0; j<3; j++)
          parallelAdd(info._HThetaD(i,j+(_p2.getVertexId()[0]+c)*3),HThetaD2->at(c)(i,j));
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::energyPDTGH(const Vec3T& v,const Vec3T& vLast,const Vec3T& vl,const Vec3T& u,const Vec4T& x,Mat3X4T* DTG,
    const Vec3T& Rvl,const Vec3T* Rvll,Mat3T* LRH,Mat3T* wLRH,Mat3T* LRH1,Mat3T* wLRH1,Mat3X4T* LRH2,Mat3X4T* wLRH2,
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,Mat3T* Hxx,GAll* g) const {
  Mat3T nntDD,C=cross<T>(Rvl),CL;
  if(Rvll) CL=cross<T>(*Rvll);
  Mat3T DDEDXDTheta;
  Mat3X4T DDEDXDTheta2;
  T e,D_=0,DD_=0,D=0,DD=0,d=vLast.dot(x.template segment<3>(0))-x[3];
  e=_p.template eval<T>(d,&D_,&DD_,_d0Half,1);
  if(e==0)
    return;
  Vec3T f_n=D_*(x.template segment<3>(0));
  Mat3T B=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm());
  Vec3T x_v=B*(v-vLast)*(1.0/_dt);
  Vec3T DU=-2*(x_v-u);
  D=0.5/sqrt((x_v-u).squaredNorm()+_eps);
  DD=-0.25/pow((x_v-u).squaredNorm()+_eps,1.5);
  Vec3T DX=2*B*(x_v-u)/_dt;
  Mat3T DDX=2*B*B.transpose()/(_dt*_dt);
  Mat3T DXDU=-2*B/_dt;
  if(DTG)
    parallelAdd<T,3,4>(*DTG,0,0,computeDTG<T>(_coef*_fri*_dt*f_n.norm()*D*DX,vl));
  if(LRH && wLRH) {
    DDEDXDTheta=_fri*_dt*f_n.norm()*(DD*DX*DU.transpose()+D*DXDU);
    parallelAdd<T,3,3>(*LRH,0,0,DDEDXDTheta);
    parallelAdd<T,3,3>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(!LRH1) {
    if(Mww && Mtw && Mwt && Mtt) {
      nntDD=_coef*_fri*_dt*f_n.norm()*(DX*DX.transpose()*DD+D*DDX);
      parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*C.transpose());
      parallelAdd<T,3,3>(*Mtw,0,0,nntDD*C.transpose());
      parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
      parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
    }
  } else {
    Vec3T DFDX=DD_*x.template segment<3>(0)*x.template segment<3>(0).transpose()*f_n/f_n.norm();
    Vec3T DXL=-2*(B/_dt)*(x_v-u);
    Mat3T DXDXL=-2*B*B/(_dt*_dt);
    Mat3T DUDXL=2*B/_dt;
    Mat3X4T CDXDN;
    ROT(CDXDN)=(x.template segment<3>(0)*vLast.transpose())*DD_+Mat3T::Identity()*D_;
    CTR(CDXDN)=-x.template segment<3>(0)*DD_;
    Vec4T DFDN=CDXDN.transpose()*f_n/f_n.norm();
    Mat3T NV=x.template segment<3>(0)*(v-vLast).transpose()/x.template segment<3>(0).norm();
    NV.diagonal().array()+=x.template segment<3>(0).dot(v-vLast)/x.template segment<3>(0).norm();
    NV*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
    Mat3T NV1=x.template segment<3>(0)*(x_v-u).transpose()/x.template segment<3>(0).norm();
    NV1.diagonal().array()+=x.template segment<3>(0).dot(x_v-u)/x.template segment<3>(0).norm();
    NV1*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
    Vec3T DNtmp=-2*NV.transpose()*(x_v-u)/_dt;
    Vec4T DN(DNtmp[0],DNtmp[1],DNtmp[2],(T)0);//right
    Mat3X4T DXDN;
    DXDN.setZero();
    ROT(DXDN)=-2*B*NV/(_dt*_dt)-2/_dt*(NV1);
    Mat3X4T DUDN;
    DUDN.setZero();
    ROT(DUDN)=2*NV/_dt;
    DDEDXDTheta=_fri*_dt*(DU*D*DFDX.transpose()+f_n.norm()*(DD*DU*DXL.transpose()+D*DUDXL));
    parallelAdd<T,3,3>(*LRH1,0,0,DDEDXDTheta);
    parallelAdd<T,3,3>(*wLRH1,0,0,DDEDXDTheta*CL.transpose());
    DDEDXDTheta2=_fri*_dt*(DX*D*DFDN.transpose()+f_n.norm()*(DD*DX*DN.transpose()+D*DXDN));
    parallelAdd<T,3,4>(*LRH2,0,0,DDEDXDTheta2);
    parallelAdd<T,3,4>(*wLRH2,0,0,C*DDEDXDTheta2);
    if(Mww && Mtw && Mwt && Mtt) {
      nntDD=_coef*_fri*_dt*(DX*D*DFDX.transpose()+f_n.norm()*(DD*DX*DXL.transpose()+D*DXDXL));
      parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*CL.transpose());
      parallelAdd<T,3,3>(*Mtw,0,0,nntDD*CL.transpose());
      parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
      parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::energyNDTGH(const Vec3T& v,const Vec3T& vLast,const Vec3T& vl,const Vec3T& u,const Vec4T& x,Mat3X4T* DTG,
    const Vec3T& Rvl,const Vec3T* Rvll,Mat3T* LRH,Mat3T* wLRH,Mat3T* LRH1,Mat3T* wLRH1,Mat3X4T* LRH2,Mat3X4T* wLRH2,
    Mat3T* Mww,Mat3T* Mtw,Mat3T* Mwt,Mat3T* Mtt,Mat3T *Hxx,GAll* g) const {
  Mat3T nntDD,C=cross<T>(Rvl),CL;
  if(Rvll) CL=cross<T>(*Rvll);
  Mat3T DDEDXDTheta;
  Mat3X4T DDEDXDTheta2;
  T e,D_=0,DD_=0,D=0,DD=0,d=x[3]-vLast.dot(x.template segment<3>(0));
  e=_p.template eval<T>(d,&D_,&DD_,_d0Half,1);
  if(e==0)
    return;
  Vec3T f_n=D_*(-x.template segment<3>(0));
  Mat3T B=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm());
  Vec3T x_v=B*(v-vLast)*(1.0/_dt);
  Vec3T DU=-2*(x_v-u);
  D=0.5/sqrt((x_v-u).squaredNorm()+_eps);
  DD=-0.25/pow((x_v-u).squaredNorm()+_eps,1.5);
  Vec3T DX=2*B*(x_v-u)/_dt;
  Mat3T DDX=2*B*B.transpose()/(_dt*_dt);
  Mat3T DXDU=-2*B/_dt;
  if(DTG)
    parallelAdd<T,3,4>(*DTG,0,0,computeDTG<T>(_coef*_fri*_dt*f_n.norm()*D*DX,vl));
  if(LRH && wLRH) {
    DDEDXDTheta=_fri*_dt*f_n.norm()*(DD*DX*DU.transpose()+D*DXDU);
    parallelAdd<T,3,3>(*LRH,0,0,DDEDXDTheta);
    parallelAdd<T,3,3>(*wLRH,0,0,C*DDEDXDTheta);
  }
  if(!LRH1) {
    if(Mww && Mtw && Mwt && Mtt) {
      nntDD=_fri*_dt*f_n.norm()*(DX*DX.transpose()*DD+D*DDX)*_coef;
      parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*C.transpose());
      parallelAdd<T,3,3>(*Mtw,0,0,nntDD*C.transpose());
      parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
      parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
    }
  } else {
    Vec3T DFDX=DD_*x.template segment<3>(0)*x.template segment<3>(0).transpose()*f_n/f_n.norm();
    Vec3T DXL=-2*(B/_dt)*(x_v-u);
    Mat3T DXDXL=-2*B*B/(_dt*_dt);
    Mat3T DUDXL=2*B/_dt;
    Mat3X4T CDXDN;
    ROT(CDXDN)=(x.template segment<3>(0)*vLast.transpose())*DD_-Mat3T::Identity()*D_;
    CTR(CDXDN)=-x.template segment<3>(0)*DD_;
    Vec4T DFDN=CDXDN.transpose()*f_n/f_n.norm();
    Mat3T NV=x.template segment<3>(0)*(v-vLast).transpose()/x.template segment<3>(0).norm();
    NV.diagonal().array()+=x.template segment<3>(0).dot(v-vLast)/x.template segment<3>(0).norm();
    NV*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
    Mat3T NV1=x.template segment<3>(0)*(x_v-u).transpose()/x.template segment<3>(0).norm();
    NV1.diagonal().array()+=x.template segment<3>(0).dot(x_v-u)/x.template segment<3>(0).norm();
    NV1*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
    Vec3T DNtmp=-2*NV.transpose()*(x_v-u)/_dt;
    Vec4T DN(DNtmp[0],DNtmp[1],DNtmp[2],(T)0);//right
    Mat3X4T DXDN;
    DXDN.setZero();
    ROT(DXDN)=-2*B*NV/(_dt*_dt)-2/_dt*(NV1);
    Mat3X4T DUDN;
    DUDN.setZero();
    ROT(DUDN)=2*NV/_dt;
    DDEDXDTheta=_fri*_dt*(DU*D*DFDX.transpose()+f_n.norm()*(DD*DU*DXL.transpose()+D*DUDXL));
    parallelAdd<T,3,3>(*LRH1,0,0,DDEDXDTheta);
    parallelAdd<T,3,3>(*wLRH1,0,0,DDEDXDTheta*CL.transpose());
    DDEDXDTheta2=_fri*_dt*(DX*D*DFDN.transpose()+f_n.norm()*(DD*DX*DN.transpose()+D*DXDN));
    parallelAdd<T,3,4>(*LRH2,0,0,DDEDXDTheta2);
    parallelAdd<T,3,4>(*wLRH2,0,0,C*DDEDXDTheta2);
    if(Mww && Mtw && Mwt && Mtt) {
      nntDD=_coef*_fri*_dt*(DX*D*DFDX.transpose()+f_n.norm()*(DD*DX*DXL.transpose()+D*DXDXL));
      parallelAdd<T,3,3>(*Mww,0,0,C*nntDD*CL.transpose());
      parallelAdd<T,3,3>(*Mtw,0,0,nntDD*CL.transpose());
      parallelAdd<T,3,3>(*Mwt,0,0,C*nntDD);
      parallelAdd<T,3,3>(*Mtt,0,0,nntDD);
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::energyPDUDN(const Vec3T& v,const Vec3T& vLast,const Vec3T& u,const Vec4T& x,Mat3X4T* LRH) const {
  T e,D_=0,DD_=0,D=0,DD=0,d=vLast.dot(x.template segment<3>(0))-x[3];
  e=_p.template eval<T>(d,&D_,&DD_,_d0Half,1);
  if(e==0)
    return;
  Vec3T f_n=D_*(x.template segment<3>(0));
  Mat3T B=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm());
  Vec3T x_v=B*(v-vLast)*(1.0/_dt);
  Vec3T DU=-2*(x_v-u);
  D=0.5/sqrt((x_v-u).squaredNorm()+_eps);
  DD=-0.25/pow((x_v-u).squaredNorm()+_eps,1.5);
  Mat3X4T CDXDN;
  ROT(CDXDN)=(x.template segment<3>(0)*vLast.transpose())*DD_+Mat3T::Identity()*D_;
  CTR(CDXDN)=-x.template segment<3>(0)*DD_;
  Vec4T DFDN=CDXDN.transpose()*f_n/f_n.norm();
  Mat3T NV=x.template segment<3>(0)*(v-vLast).transpose()/x.template segment<3>(0).norm();
  NV.diagonal().array()+=x.template segment<3>(0).dot(v-vLast)/x.template segment<3>(0).norm();
  NV*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
  Vec3T DNtmp=-2*NV.transpose()*(x_v-u)/_dt;
  Vec4T DN(DNtmp[0],DNtmp[1],DNtmp[2],(T)0);//right
  Mat3X4T DUDN;
  DUDN.setZero();
  ROT(DUDN)=2*NV/_dt;
  Mat3X4T DDEDXDTheta1=_fri*_dt*(DU*D*DFDN.transpose()+f_n.norm()*(DD*DU*DN.transpose()+D*DUDN));
  parallelAdd<T,3,4>(*LRH,0,0,DDEDXDTheta1);
}
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::energyNDUDN(const Vec3T& v,const Vec3T& vLast,const Vec3T& u,const Vec4T& x,Mat3X4T* LRH) const {
  T e,D_=0,DD_=0,D=0,DD=0,d=x[3]-vLast.dot(x.template segment<3>(0));
  e=_p.template eval<T>(d,&D_,&DD_,_d0Half,1);
  if(e==0)
    return;
  Vec3T f_n=D_*(-x.template segment<3>(0));
  Mat3T B=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/(x.template segment<3>(0)).squaredNorm());
  Vec3T x_v=B*(v-vLast)*(1.0/_dt);
  Vec3T DU=-2*(x_v-u);
  D=0.5/sqrt((x_v-u).squaredNorm()+_eps);
  DD=-0.25/pow((x_v-u).squaredNorm()+_eps,1.5);
  Mat3X4T CDXDN;
  ROT(CDXDN)=(x.template segment<3>(0)*vLast.transpose())*DD_-Mat3T::Identity()*D_;
  CTR(CDXDN)=-x.template segment<3>(0)*DD_;
  Vec4T DFDN=CDXDN.transpose()*f_n/f_n.norm();
  Mat3T NV=x.template segment<3>(0)*(v-vLast).transpose()/x.template segment<3>(0).norm();
  NV.diagonal().array()+=x.template segment<3>(0).dot(v-vLast)/x.template segment<3>(0).norm();
  NV*=(Mat3T::Identity()-x.template segment<3>(0)*x.template segment<3>(0).transpose()/x.template segment<3>(0).squaredNorm())/x.template segment<3>(0).norm();
  Vec3T DNtmp=-2*NV.transpose()*(x_v-u)/_dt;
  Vec4T DN(DNtmp[0],DNtmp[1],DNtmp[2],(T)0);//right
  Mat3X4T DUDN;
  DUDN.setZero();
  ROT(DUDN)=2*NV/_dt;
  Mat3X4T DDEDXDTheta=_fri*_dt*(DU*D*DFDN.transpose()+f_n.norm()*(DD*DU*DN.transpose()+D*DUDN));
  parallelAdd<T,3,4>(*LRH,0,0,DDEDXDTheta);
}
//gradient
template <typename T,typename PFunc,typename TH>
void CCBarrierFrictionEnergy<T,PFunc,TH>::sensitivity(const Vec3T&,const Vec3T&,const Mat3T& H,Mat3T& S) const {
  S=H.inverse();
}
//instance
template class CCBarrierFrictionEnergy<FLOAT,Px>;
template class CCBarrierFrictionEnergy<FLOAT,Logx>;
template class CCBarrierFrictionEnergy<FLOAT,CLogx>;
template class CCBarrierFrictionEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierFrictionEnergy<FLOAT,Cubicx>;
}
