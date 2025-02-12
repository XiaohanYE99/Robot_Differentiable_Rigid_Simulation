#include "ConvexHullMeshDistanceEnergy.h"
#include <Environment/ConvexHullExact.h>
#include <Environment/DistanceFunction.h>
#include "DistanceEnergy.h"
#include <stack>

namespace PHYSICSMOTION {
template <typename T,typename PFunc,typename TH>
CCBarrierMeshEnergy<T,PFunc,TH>::CCBarrierMeshEnergy(const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,const PFunc& p,T d0,const CollisionGradInfo<T>* grad,T coef,bool implicit)
  :CCBarrierEnergy<T,PFunc,TH>(p1,p2,p,d0,grad,coef,implicit) {}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,Vec4T* x) {
  if(E)
    *E=0;
  std::shared_ptr<MeshExact> c1=_p1.mesh();
  std::shared_ptr<MeshExact> c2=_p2.mesh();
  //intersection check
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
  if(dist<=_d0)
    return false;
  //energy computation
  if(_useLRI) {
    if(!evalBsh(c1,c2,E,body,grad))
      return false;
  }
  else{
    if(_useBVH) {
      if(!evalBvh(c1,c2,E,body,grad))
        return false;
    } else {
      if(!evalBF(c1,c2,E,body,grad))
        return false;
    }
  }
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
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalbackward(T *E,const ArticulatedBody* body,CollisionGradInfo<T>* grad) {
  std::shared_ptr<MeshExact> c1=_p1.mesh();
  std::shared_ptr<MeshExact> c2=_p2.mesh();
  //intersection check
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
  if(dist<=_d0)
    return false;
  //energy computation
  if(_useBVH) {
    if(!evalBvh(c1,c2,E,body,grad,true))
      return false;
  } else {
    if(!evalBF(c1,c2,E,body,grad,true))
      return false;
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::debugGradient(const GJKPolytope<T>& p,const ArticulatedBody& body,int JID,T x0,T d0,bool output) {
  T E,E2,E3,coef;
  MatT HTheta;
  PFunc barrier;
  Vec GTheta,GTheta2;
  GJKPolytope<T> p2;
  CollisionGradInfo<T> info,info2;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
  std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(JID)._mesh);
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
      CCBarrierMeshEnergy<T,PFunc,TH> e(pass?p:p2,pass?p2:p,barrier,d0,&info,coef);
      e.setOutput(output);
      e._useBVH=true;
      if(!e.eval(&E,&body,&info,NULL,&GTheta,&HTheta))
        continue;
      if(E==0)
        continue;
      //E2
      e._useBVH=false;
      e.eval(&E2,&body,&info,NULL,NULL,NULL);
      DEBUG_GRADIENT("ERef",E,E-E2)
      //E2
      e._useBVH=true;
      e.eval(&E3,&body,&info,NULL,NULL,NULL);
      DEBUG_GRADIENT("ERef2",E,E-E3)
      //evaluate again
      info2.reset(body,x+dx*DELTA);
      p2=GJKPolytope<T>(JID,body,info2);
      CCBarrierMeshEnergy<T,PFunc,TH> e2(pass?p:p2,pass?p2:p,barrier,d0,&info2,coef);
      e2.setOutput(output);
      if(!e2.eval(&E2,&body,&info2,NULL,&GTheta2,NULL))
        continue;
      DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
      DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
      break;
    }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
  T E,E2,E3,coef;
  MatT HTheta;
  PFunc barrier;
  Vec GTheta,GTheta2;
  GJKPolytope<T> p,p2;
  CollisionGradInfo<T> info,info2;
  barrier._x0=(double)x0;
  DEFINE_NUMERIC_DELTA_T(T)
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
    CCBarrierMeshEnergy<T,PFunc,TH> e(p,p2,barrier,d0,&info,coef);
    e.setOutput(output);
    e._useBVH=true;
    if(!e.eval(&E,&body,&info,NULL,&GTheta,&HTheta))
      continue;
    if(E==0)
      continue;
    //E2
    e._useBVH=false;
    e.eval(&E2,&body,&info,NULL,NULL,NULL);
    DEBUG_GRADIENT("ERef",E,E-E2)
    //E2
    e._useBVH=true;
    e.eval(&E3,&body,&info,NULL,NULL,NULL);
    DEBUG_GRADIENT("ERef2",E,E-E3)
    //evaluate again
    info2.reset(body,x+dx*DELTA);
    p=GJKPolytope<T>(JID,body,info2);
    p2=GJKPolytope<T>(JID2,body,info2);
    CCBarrierMeshEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info2,coef);
    e2.setOutput(output);
    if(!e2.eval(&E2,&body,&info2,NULL,&GTheta2,NULL,&u2))
      continue;
    DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    break;
  }
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalBF(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward) const {
  MAll m;
  //triangle to triangle
  for(int i=0; i<(int)c1->iss().size(); i++) {
    for(int j=0; j<(int)c2->iss().size(); j++) {
      //edge to edge
      for(int ei=0,cnti=3; ei<3; ei++,cnti++) {
        for(int ej=0,cntj=3; ej<3; ej++,cntj++) {
          Eigen::Matrix<int,2,1> edgei(c1->iss()[i][ei],c1->iss()[i][(ei+1)%3]);
          //edge does not belong here
          if((c1->bss()[i]&(1<<cnti))==0)
            continue;
          Eigen::Matrix<int,2,1> edgej(c2->iss()[j][ej],c2->iss()[j][(ej+1)%3]);
          //edge does not belong here
          if((c2->bss()[j]&(1<<cntj))==0)
            continue;
          GJKPolytopePtr pss[4]= {&_p1,&_p1,&_p2,&_p2};
          int vid[4]= {edgei[0],edgei[1],edgej[0],edgej[1]};
          if(!evalEE(pss,vid,E,body,grad,m,backward))
            return false;
        }
      }
      //vertex to triangle
      for(int vi=0,cnti=0; vi<3; vi++,cnti++) {
        //vertex does not belong here
        if((c1->bss()[i]&(1<<cnti))==0)
          continue;
        GJKPolytopePtr pss[4]= {&_p1,&_p2,&_p2,&_p2};
        int vid[4]= {c1->iss()[i][vi],c2->iss()[j][0],c2->iss()[j][1],c2->iss()[j][2]};
        if(!evalVT(pss,vid,E,body,grad,m,backward))
          return false;
      }
      //triangle to vertex
      for(int vj=0,cntj=0; vj<3; vj++,cntj++) {
        //vertex does not belong here
        if((c2->bss()[j]&(1<<cntj))==0)
          continue;
        GJKPolytopePtr pss[4]= {&_p2,&_p1,&_p1,&_p1};
        int vid[4]= {c2->iss()[j][vj],c1->iss()[i][0],c1->iss()[i][1],c1->iss()[i][2]};
        if(!evalVT(pss,vid,E,body,grad,m,backward))
          return false;
      }
    }
  }
  if(body && grad && !backward)
    contractHAll(*body,*grad,m);
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalBvh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward) const {
  MAll m;
  const auto& bvh1=_p1.getBVH();
  const auto& bvh2=_p2.getBVH();
  std::stack<std::pair<int,int>> ss;
  ss.push(std::make_pair((int)bvh1.size()-1,(int)bvh2.size()-1));
  while(!ss.empty()) {
    int id1=ss.top().first;
    int id2=ss.top().second;
    ss.pop();
    if(!bvh1[id1]._bb.enlarged(Vec3T::Constant(_d0+_p._x0)).intersect(bvh2[id2]._bb))
      continue;
    else if(bvh1[id1]._cell>=0 && bvh2[id2]._cell>=0) {
      //triangle to triangle
      int i=bvh1[id1]._cell;
      int j=bvh2[id2]._cell;
      //edge to edge
      for(int ei=0,cnti=3; ei<3; ei++,cnti++) {
        for(int ej=0,cntj=3; ej<3; ej++,cntj++) {
          Eigen::Matrix<int,2,1> edgei(c1->iss()[i][ei],c1->iss()[i][(ei+1)%3]);
          //edge does not belong here
          if((c1->bss()[i]&(1<<cnti))==0)
            continue;
          Eigen::Matrix<int,2,1> edgej(c2->iss()[j][ej],c2->iss()[j][(ej+1)%3]);
          //edge does not belong here
          if((c2->bss()[j]&(1<<cntj))==0)
            continue;
          GJKPolytopePtr pss[4]= {&_p1,&_p1,&_p2,&_p2};
          int vid[4]= {edgei[0],edgei[1],edgej[0],edgej[1]};
          if(!evalEE(pss,vid,E,body,grad,m,backward))
            return false;
        }
      }
      //vertex to triangle
      for(int vi=0,cnti=0; vi<3; vi++,cnti++) {
        //vertex does not belong here
        if((c1->bss()[i]&(1<<cnti))==0)
          continue;
        GJKPolytopePtr pss[4]= {&_p1,&_p2,&_p2,&_p2};
        int vid[4]= {c1->iss()[i][vi],c2->iss()[j][0],c2->iss()[j][1],c2->iss()[j][2]};
        if(!evalVT(pss,vid,E,body,grad,m,backward))
          return false;
      }
      //triangle to vertex
      for(int vj=0,cntj=0; vj<3; vj++,cntj++) {
        //vertex does not belong here
        if((c2->bss()[j]&(1<<cntj))==0)
          continue;
        GJKPolytopePtr pss[4]= {&_p2,&_p1,&_p1,&_p1};
        int vid[4]= {c2->iss()[j][vj],c1->iss()[i][0],c1->iss()[i][1],c1->iss()[i][2]};
        if(!evalVT(pss,vid,E,body,grad,m,backward))
          return false;
      }
    } else if(bvh1[id1]._cell>=0) {
      ss.push(std::make_pair(id1,bvh2[id2]._l));
      ss.push(std::make_pair(id1,bvh2[id2]._r));
    } else if(bvh2[id2]._cell>=0) {
      ss.push(std::make_pair(bvh1[id1]._l,id2));
      ss.push(std::make_pair(bvh1[id1]._r,id2));
    } else {
      ss.push(std::make_pair(bvh1[id1]._l,bvh2[id2]._l));
      ss.push(std::make_pair(bvh1[id1]._l,bvh2[id2]._r));
      ss.push(std::make_pair(bvh1[id1]._r,bvh2[id2]._l));
      ss.push(std::make_pair(bvh1[id1]._r,bvh2[id2]._r));
    }
  }
  if(body && grad && !backward)
    contractHAll(*body,*grad,m);
  return true;
}
template <typename T,typename PFunc,typename TH>
T CCBarrierMeshEnergy<T,PFunc,TH>::evalBsh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,bool backward) const {
  MAll m;
  const auto& bvh1=_p1.getBVH();
  const auto& bvh2=_p2.getBVH();
  T P1=0,P2=0,P=0;
  int id1=bvh1.size()-1;
  int id2=bvh2.size()-1;
  Vec3T x1=bvh1[id1]._bb.center();
  Vec3T x2=bvh2[id2]._bb.center();
  if(bvh1[id1]._cell>=0 && bvh2[id2]._cell>=0) {
    P=12*(1+sqrt((x1-x2).norm())).squaredNorm();
    return P;
  }
  
  
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalEE(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward) const {
  T Eee;
  Vec12T G;
  Mat12T H;
  EEBarrierEnergy<T,PFunc> ee(
    pss[0]->globalVss().col(vid[0]),
    pss[1]->globalVss().col(vid[1]),
    pss[2]->globalVss().col(vid[2]),
    pss[3]->globalVss().col(vid[3]),
    false,_p,_d0,_coef);
  //ee energy
  if(!ee.eval(&Eee,grad?&G:NULL,grad?&H:NULL))
    return false;
  if(Eee==0)
    return true;
  if(E)
    *E+=Eee;
  if(body && grad) {
    if(!backward) computeDTGH(pss,vid,*body,*grad,G,H,m);
    else computeHBackward(pss,vid,*body,*grad,G,H,m);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalVT(GJKPolytopePtr pss[4],int vid[4],T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,MAll& m,bool backward) const {
  T Evt;
  Vec12T G;
  Mat12T H;
  VTBarrierEnergy<T,PFunc> vt(
    pss[1]->globalVss().col(vid[1]),
    pss[2]->globalVss().col(vid[2]),
    pss[3]->globalVss().col(vid[3]),
    pss[0]->globalVss().col(vid[0]),
    _p,_d0,_coef);
  //vt energy
  if(!vt.eval(&Evt,grad?&G:NULL,grad?&H:NULL))
    return false;
  if(Evt==0)
    return true;
  if(E)
    *E+=Evt;
  if(body && grad) {
    if(!backward) computeDTGH(pss,vid,*body,*grad,G,H,m);
    else computeHBackward(pss,vid,*body,*grad,G,H,m);
  }
  return true;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::computeDTGH(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const {
  for(int i=0; i<4; i++) {
    if(pss[i]->jid()<0)
      continue;
    Vec3T xi=pss[i]->mesh()->vss()[vid[i]].template cast<T>();
    Mat3T Rxi=cross<T>(pss[i]->globalVss().col(vid[i])-CTRI(grad._info._TM,pss[i]->jid()));
    //G
    for(int r=0; r<3; r++)
      for(int c=0; c<4; c++)
        parallelAdd(grad._DTG(r,c+pss[i]->jid()*4),G[i*3+r]*(c<3?xi[c]:1));
    //H
    for(int j=0; j<4; j++) {
      if(pss[j]->jid()<0)
        continue;
      Mat3T Rxj=cross<T>(pss[j]->globalVss().col(vid[j])-CTRI(grad._info._TM,pss[j]->jid()));
      MPair* mp=NULL;
      if(pss[i]==&_p1 && pss[j]==&_p1)
        mp=&(m._m11);
      else if(pss[i]==&_p1 && pss[j]==&_p2)
        mp=&(m._m12);
      else if(pss[i]==&_p2 && pss[j]==&_p1)
        mp=&(m._m21);
      else if(pss[i]==&_p2 && pss[j]==&_p2)
        mp=&(m._m22);
      else {
        ASSERT(false)
      }
      mp->_Mww+=Rxi*H.template block<3,3>(i*3,j*3)*Rxj.transpose();
      mp->_Mwt+=Rxi*H.template block<3,3>(i*3,j*3);
      mp->_Mtw+=H.template block<3,3>(i*3,j*3)*Rxj.transpose();
      mp->_Mtt+=H.template block<3,3>(i*3,j*3);
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::computeHBackward(GJKPolytopePtr pss[4],const int vid[4],const ArticulatedBody& body,CollisionGradInfo<T>& grad,const Vec12T& G,const Mat12T& H,MAll& m) const {
  for(int i=0; i<4; i++) {
    if(pss[i]->jid()<0)
      continue;
    //Vec3T xi=pss[i]->mesh()->vss()[vid[i]].template cast<T>();
    Mat3T Rxi=cross<T>(pss[i]->globalVss().col(vid[i])-CTRI(grad._info._TM,pss[i]->jid()));
    //Mat3T Ri=ROTI(grad._info._TM,pss[i]->jid());
    //H
    for(int j=0; j<4; j++) {
      if(pss[j]->jid()<0)
        continue;
      MatX3T HThetaD;
      Mat3T Mwt,Mtt;
      HThetaD.setZero(body.nrDOF(),3);
      //Mat3T Rxj=cross<T>(pss[j]->globalVss().col(vid[j])-CTRI(grad._info._TM,pss[j]->jid()));
      Mat3T Rj=ROTI(grad._info._TM,pss[j]->jid());
      Mwt=Rxi*H.template block<3,3>(i*3,j*3)*Rj;//-cross<T>(G.template block<3,1>(3*i))*Ri;
      Mtt=H.template block<3,3>(i*3,j*3)*Rj;
      if(i==j) {
        Vec3T Q(G[3*i],G[3*i+1],G[3*i+2]);
        Mwt=Mwt-cross<T>(Q)*Rj;
      }
      int c=vid[j];
      if(pss[i]==&_p1) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p1.jid(),body,grad,HThetaD,Mwt,Mtt);
      } else if(pss[i]==&_p2) {
        CCBarrierEnergy<T,PFunc,TH>::contractHBackward(_p2.jid(),body,grad,HThetaD,Mwt,Mtt);
      }
      for(int ii=0; ii<body.nrDOF(); ii++)
        for(int jj=0; jj<3; jj++) {
          if(pss[j]==&_p1) parallelAdd(grad._HThetaD(ii,jj+(_p1.getVertexId()[0]+c)*3),HThetaD(ii,jj));
          else if(pss[j]==&_p2) parallelAdd(grad._HThetaD(ii,jj+(_p2.getVertexId()[0]+c)*3),HThetaD(ii,jj));
        }
    }
  }
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::contractHAll(const ArticulatedBody& body,CollisionGradInfo<T>& grad,const MAll& m) const {
  if(_p1.jid()>=0)
    CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p1.jid(),_p1.jid(),body,grad,m._m11._Mww,m._m11._Mtw,m._m11._Mwt,m._m11._Mtt);
  if(_p2.jid()>=0)
    CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p2.jid(),_p2.jid(),body,grad,m._m22._Mww,m._m22._Mtw,m._m22._Mwt,m._m22._Mtt);
  if(_p1.jid()>=0 && _p2.jid()>=0) {
    CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p1.jid(),_p2.jid(),body,grad,m._m12._Mww,m._m12._Mtw,m._m12._Mwt,m._m12._Mtt);
    CCBarrierEnergy<T,PFunc,TH>::contractHTheta(_p2.jid(),_p1.jid(),body,grad,m._m21._Mww,m._m21._Mtw,m._m21._Mwt,m._m21._Mtt);
  }
}
//instance
template class CCBarrierMeshEnergy<FLOAT,Px>;
template class CCBarrierMeshEnergy<FLOAT,Logx>;
template class CCBarrierMeshEnergy<FLOAT,CLogx>;
template class CCBarrierMeshEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierMeshEnergy<FLOAT,Cubicx>;
}
