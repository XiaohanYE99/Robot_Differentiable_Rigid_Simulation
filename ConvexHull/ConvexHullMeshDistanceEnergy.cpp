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
bool CCBarrierMeshEnergy<T,PFunc,TH>::eval(T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<Mat3X4T>* DNDX,Vec* GTheta,MatT* HTheta,std::vector<ContactManifold>* ml) {
  if(E)
    *E=0;
  std::shared_ptr<MeshExact> c1=_p1.mesh();
  std::shared_ptr<MeshExact> c2=_p2.mesh();
  //intersection check
  typename GJKPolytope<T>::Point p1,p2;
  T dist=GJKPolytope<T>::distance(_p1,_p2,p1,p2);
  if(dist<=_d0)
    return false;
  //normal,project
  Vec3T n=Vec3T(p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]).normalized();
  Vec2T r1=GJKPolytope<T>::project(_p1,n);
  Vec2T r2=GJKPolytope<T>::project(_p2,n);
  if(r1[1]<r2[0]-_d0) {
    _x=Vec4TH((TH)-n[0],(TH)-n[1],(TH)-n[2],(TH)-(r1[1]+r2[0])/2);
  } else {
    if(r2[1]>=r1[0]-_d0) {
      return false;
    }
    _x=Vec4TH((TH)n[0],(TH)n[1],(TH)n[2],(TH)(r2[1]+r1[0])/2);
  }
  //SVM margin
  TH alpha=TH(_d0Half+1)/TH(dist/2+1);
  _x*=alpha;
  //energy computation
  if(_useLRI) {
    if(!evalBsh(c1,c2,E,body,grad,ml))
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
  /*std::shared_ptr<MeshExact> c1=_p1.mesh();
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
  }*/
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
/*template <typename T,typename PFunc,typename TH>
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
    if(!e2.eval(&E2,&body,&info2,NULL,&GTheta2,NULL,NULL))
      continue;
    DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    break;
  }
}*/
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::debugGradient(const ArticulatedBody& body,int JID,int JID2,T x0,T d0,bool output) {
  T E,E2,E3,coef;
  MatT HTheta;
  PFunc barrier;
  Vec GTheta,GTheta2;
  GJKPolytope<T> p,p2;
  CollisionGradInfo<T> info,info2;
  barrier._x0=100;
  //DEFINE_NUMERIC_DELTA_T(T)
  std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(body.joint(JID)._mesh);
  std::shared_ptr<MeshExact> local2=std::dynamic_pointer_cast<MeshExact>(body.joint(JID2)._mesh);
  while(true) {
    coef=1e-10;//rand()/(T)RAND_MAX;
    Vec x,dx;
    x.setZero(12);
    /*x[0]=1.49672087e+00;
    x[1]=5.03350648e-01;
    x[2]=-1.47205678e-01;
    x[3]=1.58870264+0;
    x[4]=9.14689821e-03;
    x[5]=1.90744384e-02;
    x[6]=1.74295909e+0;
    x[7]=1.67147605e-01;
    x[8]=1.69666133e-01;
    x[9]=1.57076642e+00;
    x[10]=4.02036761e-05;
    x[11]=3.94697528e-01;*/
    x[0]=1.44468235;
    x[1]=0.45520312;
    x[2]=-0.14815516;
    x[3]=1.51494954;
    x[4]=0.03536977;
    x[5]=-0.33706095;
    x[6]=1.74049476;
    x[7]=0.16173514;
    x[8]=-0.07856184;
    x[9]=1.57002575;
    x[10]=0.00222843;
    x[11]=0.45259831;
    /*x[0]=1.38986822;
    x[1]=0.42652917;
    x[2]=-0.14815935;
    x[3]=1.51728559;
    x[4]=0.03234518;
    x[5]=-0.35279355;
    x[6]=1.70305556;
    x[7]=0.12299526;
    x[8]=-0.07876093;
    x[9]=1.56989732;
    x[10]=0.00246981;
    x[11]=0.48944889;*/
    //Vec dx=Vec::Random(body.nrDOF());
    {
      std::ifstream is("pos.dat",std::ios::binary);
      readBinaryData(x,is);
    }
    {
        std::ifstream is1("dx.dat",std::ios::binary);
        readBinaryData(dx,is1);
    }
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
    if(!e.eval(&E,&body,&info,NULL,&GTheta,&HTheta)){
      std::cout<<"#";
      continue;
    }
      
    if(E==0)
      continue;
    //E2
    /*e._useBVH=false;
    e.eval(&E2,&body,&info,NULL,NULL,NULL);
    DEBUG_GRADIENT("ERef",E,E-E2)
    //E2
    e._useBVH=true;
    e.eval(&E3,&body,&info,NULL,NULL,NULL);
    DEBUG_GRADIENT("ERef2",E,E-E3)*/
    //evaluate again
    info2.reset(body,x+dx);
    T DELTA=1.0;
    p=GJKPolytope<T>(JID,body,info2);
    p2=GJKPolytope<T>(JID2,body,info2);
    CCBarrierMeshEnergy<T,PFunc,TH> e2(p,p2,barrier,d0,&info2,coef);
    e2.setOutput(output);
    if(!e2.eval(&E2,&body,&info2,NULL,&GTheta2,NULL,NULL))
      continue;
    DEBUG_GRADIENT("dE",GTheta.dot(dx),GTheta.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("dG",(HTheta*dx).norm(),(HTheta*dx-(GTheta2-GTheta)/DELTA).norm())
    std::cout<<HTheta.norm()<<std::endl;
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
bool CCBarrierMeshEnergy<T,PFunc,TH>::evalBsh(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,T* E,const ArticulatedBody* body,CollisionGradInfo<T>* grad,std::vector<ContactManifold>* ml,bool backward) const {
  MAll m;
  GAll g;
  bool flag=true;
  Mat3X4T DTG1,DTG2;
  const auto& bvh1=_p1.getBVH();
  const auto& bvh2=_p2.getBVH();
  //T P1=0,P2=0,P=0;
  int id1=bvh1.size()-1;
  int id2=bvh2.size()-1;

  ComputePotential(c1,c2,id1,id2,E,&DTG1,&DTG2,m,g,grad,*body,ml,&flag);
  if(!flag) return false;
  if(body && grad) {
    if(_p1.jid()>=0)
      for(int r=0; r<3; r++)
        for(int c=0; c<4; c++)
          parallelAdd(grad->_DTG(r,c+_p1.jid()*4),DTG1(r,c));
    if(_p2.jid()>=0)
      for(int r=0; r<3; r++)
        for(int c=0; c<4; c++)
          parallelAdd(grad->_DTG(r,c+_p2.jid()*4),DTG2(r,c));
  }
  if(body && grad && !backward)
    contractHAll(*body,*grad,m);
  return true;
}
template <typename T,typename PFunc,typename TH>
bool CCBarrierMeshEnergy<T,PFunc,TH>::ComputePotential(std::shared_ptr<MeshExact> c1,std::shared_ptr<MeshExact> c2,int id1, int id2, T* P,
     Mat3X4T* DTG1, Mat3X4T* DTG2, MAll& m,GAll& g,CollisionGradInfo<T>* grad,
     const ArticulatedBody& body,std::vector<ContactManifold>* ml,bool* flag) const {
  if(!(*flag)) return false;
  clearMAll(m);
  clearGAll(g);
  DTG1->setZero();
  DTG2->setZero();
  *P=0;
  MAll subm;
  GAll subg;
  GAll tmpg;
  Mat3X4T subDTG1,subDTG2;
  Vec3T D1,D2;
  Mat3T Rxi,Rxj;
  Rxi.setZero();
  Rxj.setZero();
  T subP;
  const auto& bvh1=_p1.getBVH();
  const auto& bvh2=_p2.getBVH();
  Vec3T x1=bvh1[id1]._bb.center();
  Vec3T x2=bvh2[id2]._bb.center();
  Vec3T deltax=x1-x2;
  T d1=bvh1[id1]._bb.rad()+bvh2[id2]._bb.rad();
  T d2=(1.0+_ep)*d1;
  //T d=deltax.norm()+_d0;
  T dist=(deltax.norm()-d1)/(d2-d1);
  //dist=0;
  T alpha=12*_coef;//*bvh1[id1]._num*bvh2[id2]._num;
  T Q=sqrt(deltax.norm());//-_d0;
  T phi=0,dphi=0,ddphi=0;
  T P2=alpha*pow((1+1.0/Q),2);
  T D=-2.0/(Q*Q)-2.0/(Q*Q*Q);
  T DD=4.0/(pow(Q,3))+6/(pow(Q,4));
  Vec3T DQ;
  Mat3T h,DDP,DDPhi;
  Mat6T H;
  h.setZero();
  H.setZero();
  DDP.setZero();
  DDPhi.setZero();
  Vec3T vl1=c1->getBVH()[id1]._bb.center();
  Vec3T vl2=c2->getBVH()[id2]._bb.center();
  //if(_p1.jid()==-1) std::cout<<x1;
  if(grad){
    if(_p1.jid()>=0) Rxi=cross<T>(ROTI(grad->_info._TM,_p1.jid())*vl1);
    if(_p2.jid()>=0) Rxj=cross<T>(ROTI(grad->_info._TM,_p2.jid())*vl2);
  }
  if(dist>1) {
    *P=P2;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,computeDTG<T>(alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,computeDTG<T>(-alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl2));
      contractGAll(g,Rxi,Rxj,alpha*D*deltax/(2*pow(deltax.norm(),1.5)));
      DQ=deltax/(2*pow(deltax.norm(),1.5));
      DDP=DD*DQ*DQ.transpose()+D*(Mat3T::Identity()/(2*pow(deltax.norm(),1.5))-3*deltax*deltax.transpose()/(4*pow(deltax.norm(),3.5)));
      h=alpha*DDP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
    }
    return *flag;
  }
  else if(dist>0) {
    phi=6*pow(dist,5)-15*pow(dist,4)+10*pow(dist,3);
    dphi=30*pow(dist,4)-60*pow(dist,3)+30*pow(dist,2);
    ddphi=120*pow(dist,3)-180*pow(dist,2)+60*dist;
  }
  if(bvh1[id1]._cell>=0 && bvh2[id2]._cell>=0) {
    *P+=phi*P2;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,P2*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,P2*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      parallelAdd<T,3,4>(*DTG1,0,0,computeDTG<T>(phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,computeDTG<T>(-phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl2));
      contractGAll(g,Rxi,Rxj,P2*dphi*deltax/((d2-d1)*deltax.norm())+phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)));
      DQ=deltax/(2*pow(deltax.norm(),1.5));
      DDP=DD*DQ*DQ.transpose()+D*(Mat3T::Identity()/(2*pow(deltax.norm(),1.5))-3*deltax*deltax.transpose()/(4*pow(deltax.norm(),3.5)));
      DDPhi=ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1);
      h=alpha*(DDPhi*P2/alpha+dphi*deltax/((d2-d1)*deltax.norm())*D*DQ.transpose()
              +D*DQ*dphi*deltax.transpose()/((d2-d1)*deltax.norm())+DDP*phi);
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
    }
    
    std::vector<Eigen::Matrix<double,3,1>> vssA,vssB;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    iss.push_back(Eigen::Matrix<int,3,1>(0,1,2));
    for(auto &pair : c1->getBVH()[id1]._bb._points) vssA.push_back(pair.second.template cast<double>());
    for(auto &pair : c2->getBVH()[id2]._bb._points) vssB.push_back(pair.second.template cast<double>()); 
    std::shared_ptr<MeshExact> m1(new MeshExact(vssA,iss,false,false));
    std::shared_ptr<MeshExact> m2(new MeshExact(vssB,iss,false,false));
    GJKPolytope<T>mA,mB;
    mA=GJKPolytope<T>(_p1.jid(),m1,*_grad);
    mB=GJKPolytope<T>(_p2.jid(),m2,*_grad);
    CCBarrierConvexEnergy<T,PFunc> cc(mA,mB,_p,_d0,_grad,_coef);
    ContactManifold c;
    T val=0;
    clearMAll(subm);
    clearGAll(subg);
    subDTG1.setZero();
    subDTG2.setZero();
    if(!cc.evalLRI(&val,&body,grad,&c._DNDX,NULL,NULL,&subDTG1,&subDTG2,subm,subg,_x)) return *flag=false;
    *P+=(1-phi)*val;
    if(ml){
      c._jidA=_p1.jid();
      c._jidB=_p2.jid();
      c._sA=m1;
      c._sB=m2;
      c._x=cc.getX();
      ml->push_back(c);
    }
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,val*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,val*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,val*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*val;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }
    return *flag;
  }
  else if(bvh1[id1]._cell>=0) {
    ComputePotential(c1,c2,id1,bvh2[id2]._l,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }

    ComputePotential(c1,c2,id1,bvh2[id2]._r,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }
  }
  else if(bvh2[id2]._cell>=0) {
    ComputePotential(c1,c2,bvh1[id1]._l,id2,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);      
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }

    ComputePotential(c1,c2,bvh1[id1]._r,id2,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }
  }
  else {
    ComputePotential(c1,c2,bvh1[id1]._l,bvh2[id2]._l,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }

    ComputePotential(c1,c2,bvh1[id1]._l,bvh2[id2]._r,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }

    ComputePotential(c1,c2,bvh1[id1]._r,bvh2[id2]._l,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }

    ComputePotential(c1,c2,bvh1[id1]._r,bvh2[id2]._r,&subP,&subDTG1,&subDTG2,subm,subg,grad,body,ml,flag);
    *P+=(1-phi)*subP;
    if(grad){
      parallelAdd<T,3,4>(*DTG1,0,0,(1-phi)*subDTG1);
      parallelAdd<T,3,4>(*DTG2,0,0,(1-phi)*subDTG2);
      parallelAdd<T,3,4>(*DTG1,0,0,subP*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl1));
      parallelAdd<T,3,4>(*DTG2,0,0,subP*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl2));
      addGAll(g,subg,1-phi);
      contractGAll(g,Rxi,Rxj,subP*-dphi*deltax/((d2-d1)*deltax.norm()));
      DDPhi=-(ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1));
      h=DDPhi*subP;
      H.setZero();
      H.template block<3,3>(0,0)=h;
      H.template block<3,3>(0,3)=-h;
      H.template block<3,3>(3,0)=-h;
      H.template block<3,3>(3,3)=h;
      contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
      clearGAll(tmpg);
      contractGAll(tmpg,Rxi,Rxj,-dphi*deltax/((d2-d1)*deltax.norm()));
      mergeGAll(tmpg,subg,m);
      mergeGAll(subg,tmpg,m);
      addMAll(m,subm,1-phi);
    }
  }
  *P+=phi*P2;
  if(grad){
    parallelAdd<T,3,4>(*DTG1,0,0,P2*computeDTG<T>(dphi*deltax/((d2-d1)*deltax.norm()),vl1));
    parallelAdd<T,3,4>(*DTG2,0,0,P2*computeDTG<T>(-dphi*deltax/((d2-d1)*deltax.norm()),vl2));
    parallelAdd<T,3,4>(*DTG1,0,0,computeDTG<T>(phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl1));
    parallelAdd<T,3,4>(*DTG2,0,0,computeDTG<T>(-phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)),vl2));
    contractGAll(g,Rxi,Rxj,P2*dphi*deltax/((d2-d1)*deltax.norm())+phi*alpha*D*deltax/(2*pow(deltax.norm(),1.5)));
    DQ=deltax/(2*pow(deltax.norm(),1.5));
    DDP=DD*DQ*DQ.transpose()+D*(Mat3T::Identity()/(2*pow(deltax.norm(),1.5))-3*deltax*deltax.transpose()/(4*pow(deltax.norm(),3.5)));
    DDPhi=ddphi*(deltax/deltax.norm())*(deltax/deltax.norm()).transpose()/(pow(d2-d1,2))+dphi*(Mat3T::Identity()/deltax.norm()-deltax*deltax.transpose()/pow(deltax.norm(),3))/(d2-d1);
    h=alpha*(DDPhi*P2/alpha+dphi*deltax/((d2-d1)*deltax.norm())*D*DQ.transpose()
            +D*DQ*dphi*deltax.transpose()/((d2-d1)*deltax.norm())+DDP*phi);
    H.setZero();
    H.template block<3,3>(0,0)=h;
    H.template block<3,3>(0,3)=-h;
    H.template block<3,3>(3,0)=-h;
    H.template block<3,3>(3,3)=h;
    contractMAll(m,Rxi,Rxj,Rxi,Rxj,H);
  }
  /*Vec dx=Vec::Random(3);
  Vec3T g=dphi*deltax/deltax.norm()/(d2-d1)*pow((1+1.0/Q),2)+deltax/(2*pow(deltax.norm(),1.5))*phi;
  deltax+=dx*1e-8;
  dist=(deltax.norm()-d1)/(d2-d1);
  dphi=30*pow(dist,4)-60*pow(dist,3)+30*pow(dist,2);
  phi=6*pow(dist,5)-15*pow(dist,4)+10*pow(dist,3);
  Q=sqrt(deltax.norm())+_d0;
  Vec3T g2=dphi*deltax/deltax.norm()/(d2-d1)*pow((1+1.0/Q),2)+deltax/(2*pow(deltax.norm(),1.5))*phi;
  std::cout<<(H*dx).norm()<<" "<<((H*dx)-(g2-g)/1e-8).norm()<<std::endl;*/
  return *flag;
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
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::contractGAll(GAll& g,Mat3T Rxi,Mat3T Rxj,Vec3T L)const{
  g._g1._w+=Rxi*L;
  g._g1._t+=L;
  g._g2._w-=Rxj*L;
  g._g2._t-=L;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::clearGAll(GAll& g) const {
  g._g1._w.setZero();
  g._g1._t.setZero();
  g._g2._w.setZero();
  g._g2._t.setZero();
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::addGAll(GAll& g1,GAll& g2,T alpha) const {
  g1._g1._w+=g2._g1._w*alpha;
  g1._g1._t+=g2._g1._t*alpha;
  g1._g2._w+=g2._g2._w*alpha;
  g1._g2._t+=g2._g2._t*alpha;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::contractMAll(MAll& m,Mat3T Rxi,Mat3T Rxj,Mat3T rxi,Mat3T rxj,Mat6T H) const {
  m._m11._Mww+=Rxi*H.template block<3,3>(0,0)*rxi.transpose();
  m._m11._Mwt+=Rxi*H.template block<3,3>(0,0);
  m._m11._Mtw+=H.template block<3,3>(0,0)*rxi.transpose();
  m._m11._Mtt+=H.template block<3,3>(0,0);
  m._m12._Mww+=Rxi*H.template block<3,3>(0,3)*rxj.transpose();
  m._m12._Mwt+=Rxi*H.template block<3,3>(0,3);
  m._m12._Mtw+=H.template block<3,3>(0,3)*rxj.transpose();
  m._m12._Mtt+=H.template block<3,3>(0,3);
  m._m21._Mww+=Rxj*H.template block<3,3>(3,0)*rxi.transpose();
  m._m21._Mwt+=Rxj*H.template block<3,3>(3,0);
  m._m21._Mtw+=H.template block<3,3>(3,0)*rxi.transpose();
  m._m21._Mtt+=H.template block<3,3>(3,0);
  m._m22._Mww+=Rxj*H.template block<3,3>(3,3)*rxj.transpose();
  m._m22._Mwt+=Rxj*H.template block<3,3>(3,3);
  m._m22._Mtw+=H.template block<3,3>(3,3)*rxj.transpose();
  m._m22._Mtt+=H.template block<3,3>(3,3);
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::addMAll(MAll& m1,MAll& m2,T alpha) const {
  m1._m11._Mww+=m2._m11._Mww*alpha;
  m1._m11._Mwt+=m2._m11._Mwt*alpha;
  m1._m11._Mtw+=m2._m11._Mtw*alpha;
  m1._m11._Mtt+=m2._m11._Mtt*alpha;
  m1._m12._Mww+=m2._m12._Mww*alpha;
  m1._m12._Mwt+=m2._m12._Mwt*alpha;
  m1._m12._Mtw+=m2._m12._Mtw*alpha;
  m1._m12._Mtt+=m2._m12._Mtt*alpha;
  m1._m21._Mww+=m2._m21._Mww*alpha;
  m1._m21._Mwt+=m2._m21._Mwt*alpha;
  m1._m21._Mtw+=m2._m21._Mtw*alpha;
  m1._m21._Mtt+=m2._m21._Mtt*alpha;
  m1._m22._Mww+=m2._m22._Mww*alpha;
  m1._m22._Mwt+=m2._m22._Mwt*alpha;
  m1._m22._Mtw+=m2._m22._Mtw*alpha;
  m1._m22._Mtt+=m2._m22._Mtt*alpha;
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::clearMAll(MAll& m) const {
  m._m11._Mww.setZero();
  m._m11._Mwt.setZero();
  m._m11._Mtw.setZero();
  m._m11._Mtt.setZero();
  m._m12._Mww.setZero();
  m._m12._Mwt.setZero();
  m._m12._Mtw.setZero();
  m._m12._Mtt.setZero();
  m._m21._Mww.setZero();
  m._m21._Mwt.setZero();
  m._m21._Mtw.setZero();
  m._m21._Mtt.setZero();
  m._m22._Mww.setZero();
  m._m22._Mwt.setZero();
  m._m22._Mtw.setZero();
  m._m22._Mtt.setZero();
}
template <typename T,typename PFunc,typename TH>
void CCBarrierMeshEnergy<T,PFunc,TH>::mergeGAll(GAll& g1,GAll& g2,MAll& m) const {
  m._m11._Mww+=g1._g1._w*g2._g1._w.transpose();
  m._m11._Mwt+=g1._g1._w*g2._g1._t.transpose();
  m._m11._Mtw+=g1._g1._t*g2._g1._w.transpose();
  m._m11._Mtt+=g1._g1._t*g2._g1._t.transpose();
  m._m12._Mww+=g1._g1._w*g2._g2._w.transpose();
  m._m12._Mwt+=g1._g1._w*g2._g2._t.transpose();
  m._m12._Mtw+=g1._g1._t*g2._g2._w.transpose();
  m._m12._Mtt+=g1._g1._t*g2._g2._t.transpose();
  m._m21._Mww+=g1._g2._w*g2._g1._w.transpose();
  m._m21._Mwt+=g1._g2._w*g2._g1._t.transpose();
  m._m21._Mtw+=g1._g2._t*g2._g1._w.transpose();
  m._m21._Mtt+=g1._g2._t*g2._g1._t.transpose();
  m._m22._Mww+=g1._g2._w*g2._g2._w.transpose();
  m._m22._Mwt+=g1._g2._w*g2._g2._t.transpose();
  m._m22._Mtw+=g1._g2._t*g2._g2._w.transpose();
  m._m22._Mtt+=g1._g2._t*g2._g2._t.transpose();
}

//instance
template class CCBarrierMeshEnergy<FLOAT,Px>;
template class CCBarrierMeshEnergy<FLOAT,Logx>;
template class CCBarrierMeshEnergy<FLOAT,CLogx>;
template class CCBarrierMeshEnergy<FLOAT,InvQuadraticx>;
template class CCBarrierMeshEnergy<FLOAT,Cubicx>;
}
