#include "ConvexHullMeshDistanceFrictionEnergy.h"
#include "ConvexHullMeshDistanceEnergy.h"
#include <Environment/DistanceFunction.h>
#include <Environment/ConvexHullExact.h>
#include <Utils/CrossSpatialUtils.h>
#include <Utils/SparseUtils.h>

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH>
CCBarrierMeshFrictionEnergy<T,PFunc,TH>::FrictionTerm::FrictionTerm() {}
template <typename T,typename PFunc,typename TH>
CCBarrierMeshFrictionEnergy<T,PFunc,TH>::FrictionTerm::FrictionTerm(const Eigen::Matrix<int,4,1>& pss,const Eigen::Matrix<int,4,1>& vid,const Vec3T V[4],const Vec4T& bary,const Vec3T& n,T D)
  :_pss(pss),_vid(vid),_bary(bary),_n(n),_D(D) {
  for(int d=0; d<4; d++)
    _V[d]=V[d];
}
//CCBarrierMeshFrictionEnergy
template <typename T,typename PFunc,typename TH>
CCBarrierMeshFrictionEnergy<T,PFunc,TH>::CCBarrierMeshFrictionEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const GJKPolytope<T>& pl1,const GJKPolytope<T>& pl2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,T dt,const std::vector<FrictionTerm>* terms)
  :CCBarrierMeshEnergy<T,PFunc,TH>(p1,p2,p,d0,grad,coef,false),_pl1(pl1),_pl2(pl2),_dt(dt),_eps(1e-4f),_fri(1.f) {
  if(terms) {
    _terms=*terms;
    return;
  }
  std::shared_ptr<MeshExact> c1=_p1.mesh();
  std::shared_ptr<MeshExact> c2=_p2.mesh();
  //triangle to triangle
  for(int i=0; i<(int)c1->iss().size(); i++) {
    for(int j=0; j<(int)c2->iss().size(); j++) {
      //edge to edge
      for(int ei=0,cnti=3; ei<3; ei++,cnti++) {
        for(int ej=0,cntj=3; ej<3; ej++,cntj++) {
          int v0=c1->iss()[i][ei];
          int v1=c1->iss()[i][(ei+1)%3];
          //edge does not belong here
          if((c1->bss()[i]&(1<<cnti))==0)
            continue;
          int v2=c2->iss()[j][ej];
          int v3=c2->iss()[j][(ej+1)%3];
          //edge does not belong here
          if((c2->bss()[j]&(1<<cntj))==0)
            continue;
          Vec2T bary;
          Vec3T cpa,cpb,n,V[4]= {
            _pl1.globalVss().col(v0),
            _pl1.globalVss().col(v1),
            _pl2.globalVss().col(v2),
            _pl2.globalVss().col(v3),
          };
          T dSqr=distToSqrLineSegment<T>(V,V+2,Eigen::Map<Vec2T>(bary.data()),cpa,cpb,NULL),d=sqrt(dSqr),D=0;
          if(_p.template eval<T>(d,&D,NULL,_d0,_coef)==0)
            continue;
          //construct term
          Vec4T b(1-bary[0],bary[0],-(1-bary[1]),-bary[1]);
          n=V[0]*b[0]+V[1]*b[1]+V[2]*b[2]+V[3]*b[3];
          _terms.push_back(FrictionTerm(Eigen::Matrix<int,4,1>(1,1,2,2),Eigen::Matrix<int,4,1>(v0,v1,v2,v3),V,b,n.normalized(),abs(D)));
        }
      }
      //vertex to triangle
      for(int vi=0,cnti=0; vi<3; vi++,cnti++) {
        //vertex does not belong here
        if((c1->bss()[i]&(1<<cnti))==0)
          continue;
        int v0=c1->iss()[i][vi];
        int v1=c2->iss()[j][0];
        int v2=c2->iss()[j][1];
        int v3=c2->iss()[j][2];
        Vec3T bary,cp,n,V[4]= {
          _pl1.globalVss().col(v0),
          _pl2.globalVss().col(v1),
          _pl2.globalVss().col(v2),
          _pl2.globalVss().col(v3),
        };
        T dSqr=distToSqrTriangle<T>(V[0],V+1,Eigen::Map<Vec3T>(bary.data()),cp,NULL),d=sqrt(dSqr),D=0;
        if(_p.template eval<T>(d,&D,NULL,_d0,_coef)==0)
          continue;
        //construct term
        Vec4T b(-1,bary[0],bary[1],bary[2]);
        n=V[0]*b[0]+V[1]*b[1]+V[2]*b[2]+V[3]*b[3];
        _terms.push_back(FrictionTerm(Eigen::Matrix<int,4,1>(1,2,2,2),Eigen::Matrix<int,4,1>(v0,v1,v2,v3),V,b,n.normalized(),abs(D)));
      }
      //triangle to vertex
      for(int vj=0,cntj=0; vj<3; vj++,cntj++) {
        //vertex does not belong here
        if((c2->bss()[j]&(1<<cntj))==0)
          continue;
        int v0=c2->iss()[j][vj];
        int v1=c1->iss()[i][0];
        int v2=c1->iss()[i][1];
        int v3=c1->iss()[i][2];
        Vec3T bary,cp,n,V[4]= {
          _pl2.globalVss().col(v0),
          _pl1.globalVss().col(v1),
          _pl1.globalVss().col(v2),
          _pl1.globalVss().col(v3),
        };
        T dSqr=distToSqrTriangle<T>(V[0],V+1,Eigen::Map<Vec3T>(bary.data()),cp,NULL),d=sqrt(dSqr),D=0;
        if(_p.template eval<T>(d,&D,NULL,_d0,_coef)==0)
          continue;
        //construct term
        Vec4T b(-1,bary[0],bary[1],bary[2]);
        n=V[0]*b[0]+V[1]*b[1]+V[2]*b[2]+V[3]*b[3];
        _terms.push_back(FrictionTerm(Eigen::Matrix<int,4,1>(2,1,1,1),Eigen::Matrix<int,4,1>(v0,v1,v2,v3),V,b,n.normalized(),abs(D)));
      }
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshFrictionEnergy<T,PFunc,TH>::debugGradient(const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0,bool output) {
  T E,E2,ETmp,coef,dt=0.01f;
  PFunc barrier;
  Vec GTheta,GTheta2;
  MatT HTheta,HThetaL;
  GJKPolytope<T> p2,pl2;
  const GJKPolytope<T>& pl=p;
  CollisionGradInfo<T> info,infoL;
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
      CCBarrierMeshEnergy<T,PFunc,TH> el(pass?pl:pl2,pass?pl2:pl,barrier,d0,&infoL,coef);
      if(!el.eval(&ETmp,&body,&infoL,NULL,NULL,NULL))
        continue;
      if(ETmp==0)
        continue;
      CCBarrierMeshFrictionEnergy<T,PFunc,TH> e(pass?p:p2,pass?p2:p,pass?pl:pl2,pass?pl2:pl,barrier,d0,&info,coef,dt);
      e.setOutput(output);
      if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
        continue;
      //eval forward
      info.reset(body,x+dx*DELTA);
      p2=GJKPolytope<T>(JID,body,info);
      //calculate separating plane
      CCBarrierMeshFrictionEnergy<T,PFunc,TH> e2(pass?p:p2,pass?p2:p,pass?pl:pl2,pass?pl2:pl,barrier,d0,&info,coef,dt);
      e2.setOutput(output);
      if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
      break;
    }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshFrictionEnergy<T,PFunc,TH>::debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
  T E,E2,ETmp,coef,dt=0.01f;
  PFunc barrier;
  Vec GTheta,GTheta2;
  MatT HTheta,HThetaL;
  GJKPolytope<T> p,pl,p2,pl2;
  CollisionGradInfo<T> info,infoL;
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
    CCBarrierMeshEnergy<T,PFunc,TH> el(pl,pl2,barrier,d0,&infoL,coef);
    if(!el.eval(&ETmp,&body,&infoL,NULL,NULL,NULL))
      continue;
    if(ETmp==0)
      continue;
    CCBarrierMeshFrictionEnergy<T,PFunc,TH> e(p,p2,pl,pl2,barrier,d0,&info,coef,dt);
    e.setOutput(output);
    if(!e.eval(&E,&body,&info,&GTheta,&HTheta))
      continue;
    //eval forward
    info.reset(body,x+dx*DELTA);
    infoL.reset(body,xL);
    p=GJKPolytope<T>(JID,body,info);
    pl=GJKPolytope<T>(JID,body,infoL);
    p2=GJKPolytope<T>(JID2,body,info);
    pl2=GJKPolytope<T>(JID2,body,infoL);
    //calculate separating plane
    CCBarrierMeshFrictionEnergy<T,PFunc,TH> e2(p,p2,pl,pl2,barrier,d0,&info,coef,dt);
    e2.setOutput(output);
    if(!e2.eval(&E2,&body,&info,&GTheta2,NULL))
      continue;
    DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    break;
  }
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshFrictionEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,Vec* GTheta,MatT* HTheta) {
  if(E)
    *E=0;
  MAll m;
  for(const FrictionTerm& term:_terms) {
    T coef=_dt*term._D*_fri;
    Vec3T dv=Vec3T::Zero();
    Mat3T P=Mat3T::Identity()-term._n*term._n.transpose();
    for(int d=0; d<4; d++) {
      const GJKPolytope<T>& p=term._pss[d]==1?_p1:_p2;
      dv+=term._bary[d]*P*(p.globalVss().col(term._vid[d])-term._V[d]);
    }
    T vNorm=sqrt((dv/_dt).squaredNorm()+_eps);
    if(E)
      *E+=coef*vNorm;
    if(body && grad) {
      T tmp=(vNorm*_dt*_dt);
      Vec3T G=P*dv*coef/tmp;
      Mat3T H=P*(Mat3T::Identity()*coef/tmp-coef*dv*dv.transpose()/(tmp*tmp*vNorm))*P;
      if(body && grad)
        computeDTGH(term,*body,*grad,G,H,m);
    }
  }
  if(body && grad)
    CCBarrierMeshEnergy<T,PFunc,TH>::contractHAll(*body,*grad,m);
  //compute gradient and hessian
  if(body && grad) {
    Mat3XT DTG;
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
bool CCBarrierMeshFrictionEnergy<T,PFunc,TH>::evalbackward(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad) {
  if(E)
    *E=0;
  MAll m;
  for(const FrictionTerm& term:_terms) {
    T coef=_dt*term._D*_fri;
    Vec3T dv=Vec3T::Zero();
    Mat3T P=Mat3T::Identity()-term._n*term._n.transpose();
    for(int d=0; d<4; d++) {
      const GJKPolytope<T>& p=term._pss[d]==1?_p1:_p2;
      dv+=term._bary[d]*P*(p.globalVss().col(term._vid[d])-term._V[d]);
    }
    T vNorm=sqrt((dv/_dt).squaredNorm()+_eps);
    if(E)
      *E+=coef*vNorm;
    if(body && grad) {
      T tmp=(vNorm*_dt*_dt);
      Vec3T G=P*dv*coef/tmp;
      Mat3T H=P*(Mat3T::Identity()*coef/tmp-coef*dv*dv.transpose()/(tmp*tmp*vNorm))*P;
      if(body && grad)
        computeHBackward(term,*body,*grad,G,H,m);
    }
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
const std::vector<typename CCBarrierMeshFrictionEnergy<T,PFunc,TH>::FrictionTerm>& CCBarrierMeshFrictionEnergy<T,PFunc,TH>::terms() const {
  return _terms;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshFrictionEnergy<T,PFunc,TH>::computeDTGH(const FrictionTerm& term,const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec3T& G,const Mat3T& H,MAll& m) const {
  for(int i=0; i<4; i++) {
    GJKPolytopePtr pssi=term._pss[i]==1?&_p1:&_p2;
    if(pssi->jid()<0)
      continue;
    Vec3T xi=pssi->mesh()->vss()[term._vid[i]].template cast<T>();
    Mat3T Rxi=cross<T>(pssi->globalVss().col(term._vid[i])-CTRI(grad._info._TM,pssi->jid()));
    //G
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        parallelAdd(grad._DTG(r,c+pssi->jid()*4),G[r]*(c<3?xi[c]:1)*term._bary[i]);
    //H
    for(int j=0; j<4; j++) {
      GJKPolytopePtr pssj=term._pss[j]==1?&_p1:&_p2;
      if(pssj->jid()<0)
        continue;
      Mat3T Rxj=cross<T>(pssj->globalVss().col(term._vid[j])-CTRI(grad._info._TM,pssj->jid()));
      MPair* mp=NULL;
      if(pssi==&_p1 && pssj==&_p1)
        mp=&(m._m11);
      else if(pssi==&_p1 && pssj==&_p2)
        mp=&(m._m12);
      else if(pssi==&_p2 && pssj==&_p1)
        mp=&(m._m21);
      else if(pssi==&_p2 && pssj==&_p2)
        mp=&(m._m22);
      else {
        ASSERT(false)
      }
      mp->_Mww+=Rxi*H*Rxj.transpose()*term._bary[i]*term._bary[j];
      mp->_Mwt+=Rxi*H*term._bary[i]*term._bary[j];
      mp->_Mtw+=H*Rxj.transpose()*term._bary[i]*term._bary[j];
      mp->_Mtt+=H*term._bary[i]*term._bary[j];
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshFrictionEnergy<T,PFunc,TH>::computeHBackward(const FrictionTerm& term,const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec3T& G,const Mat3T& H,MAll& m) const {
  for(int i=0; i<4; i++) {
    GJKPolytopePtr pssi=term._pss[i]==1?&_p1:&_p2;
    if(pssi->jid()<0)
      continue;
    Vec3T xi=pssi->mesh()->vss()[term._vid[i]].template cast<T>();
    Mat3T Rxi=cross<T>(pssi->globalVss().col(term._vid[i])-CTRI(grad._info._TM,pssi->jid()));
    //Mat3T Ri=ROTI(grad._info._TM,pssi->jid());
    //G
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        parallelAdd(grad._DTG(r,c+pssi->jid()*4),G[r]*(c<3?xi[c]:1)*term._bary[i]);
    //H
    for(int j=0; j<4; j++) {
      GJKPolytopePtr pssj=term._pss[j]==1?&_p1:&_p2;
      if(pssj->jid()<0)
        continue;
      MatX3T HThetaD;
      Mat3T Mwt,Mtt;
      HThetaD.setZero(body.nrDOF(),3);
      //Mat3T Rxj=cross<T>(pssj->globalVss().col(term._vid[j])-CTRI(grad._info._TM,pssj->jid()));
      Mat3T Rj=ROTI(grad._info._TM,pssj->jid());
      Mwt=Rxi*H*term._bary[i]*term._bary[j]*Rj;//-cross<T>(G.template block<3,1>(3*i))*Ri;
      Mtt=H*term._bary[i]*term._bary[j]*Rj;
      if(i==j) {
        Vec3T Q=G*term._bary[i];
        Mwt=Mwt-cross<T>(Q)*Rj;
      }
      int c=term._vid[j];
      if(pssi==&_p1) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,grad,HThetaD,Mwt,Mtt);
      } else if(pssi==&_p2) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,grad,HThetaD,Mwt,Mtt);
      }
      for(int ii=0; ii<body.nrDOF(); ii++)
        for(int jj=0; jj<3; jj++) {
          if(pssj==&_p1) parallelAdd(grad._HThetaD(ii,jj+(_p1.getVertexId()[0]+c)*3),HThetaD(ii,jj));
          else if(pssj==&_p2) parallelAdd(grad._HThetaD(ii,jj+(_p2.getVertexId()[0]+c)*3),HThetaD(ii,jj));
        }
    }
  }
}
//instance
template class CCBarrierMeshFrictionEnergy<FLOAT,Px>;
template class CCBarrierMeshFrictionEnergy<FLOAT,Logx>;
template class CCBarrierMeshFrictionEnergy<FLOAT,CLogx>;
template class CCBarrierMeshFrictionEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierMeshFrictionEnergy<FLOAT,Cubicx>;
}
