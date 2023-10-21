#include "ConvHullPBDSimulator.h"
#include "PBDMatrixSolver.h"
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
#define USE_CRBA 0
ConvHullPBDSimulator::ConvHullPBDSimulator(T dt):Simulator(dt),_gTol(1e-7f),_alpha(1e-6f),_epsV(1e-1f),_output(true),_JTJ(true),_crossTerm(false),_maxIt(1e4) {}
void ConvHullPBDSimulator::setArticulatedBody(std::shared_ptr<ArticulatedBody> body) {
  Simulator::setArticulatedBody(body);
  _JRCF.setZero(3,4*_body->nrJ());
  _pos.reset(*body,Vec::Zero(body->nrDOF()));
  _lastPos.reset(*body,Vec::Zero(body->nrDOF()));
  setCrossTerm(_crossTerm);
}
void ConvHullPBDSimulator::step() {
  if(_body->nrDOF()==0)
    return;
  //generate contacts
  
  //detectCurrentContact();
  //update kinematic state
  GradInfo newPos(*_body,setKinematic(_pos._xM,_t+_dt)),newPos2;
  //normal solve
  Vec D,DE,DE2;
  MatT DDE,DDE2;
  T alphaMax=1e6f,alphaMin=1e-6f;
  CollisionGradInfo<T> grad,grad2;//(*_body,setKinematic(_pos._xM,_t+_dt));
  grad.reset(*_body,newPos._xM);
  //std::cout<<"e:"<<std::endl;
  T e=energy(grad,DE,DDE,true),e2=0;
  mask(_diag,DE,DDE);
  //_alpha=alphaMin;
  for(int iter=0; iter<_maxIt;) {
    //update configuration
    
    update(newPos,newPos2,D,DE,DDE,_alpha);
    //CollisionGradInfo<T> grad2(*_body,setKinematic(newPos2._xM,_t+_dt));
    grad2.reset(*_body,newPos2._xM);
    //detectCurrentContact();
    
    e2=energy(grad2,DE2,DDE2);
    //std::cout<<"e:"<<e<<"e2:"<<e2<<std::endl;
    mask(_diag,DE2,DDE2);
    //iteration update
    //std::cout<<"ispenetrated: "<<_ispenetrated<<std::endl;
    if((e2<e || DE2.cwiseAbs().maxCoeff()<_gTol)){//
      _alpha=std::max<T>(_alpha*.8f,alphaMin);
      newPos=newPos2;
      e=e2;
      DE=DE2;
      DDE=DDE2;
      iter++;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " gNorm=" << DE.cwiseAbs().maxCoeff() << " alpha=" << _alpha << std::endl;
      //termination: gradient tolerance
      if(DE2.cwiseAbs().maxCoeff()<_gTol)
        break;
    } else {
      _alpha=std::min<T>(_alpha*1.5f,alphaMax);
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " E2=" << e2 << " alpha=" << _alpha << " gNorm=" << DE2.cwiseAbs().maxCoeff()<<std::endl;
      //termination: numerical issue
      if(_alpha>=alphaMax)
        break;
    }
  }
  //update
  _lastPos=_pos;
  _pos=newPos;
  _t+=_dt;
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::pos() const {
  return _pos._xM;
}
void ConvHullPBDSimulator::setPos(const Vec& pos) {
  _pos.reset(*_body,setKinematic(pos,_t));
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::vel() const {
  return (_pos._xM-_lastPos._xM)/_dt;
}
void ConvHullPBDSimulator::setVel(const Vec& vel) {
  _lastPos.reset(*_body,setKinematic(_pos._xM-vel*_dt,_t-_dt));
}
void ConvHullPBDSimulator::setGravity(const Vec3T& g) {
  _JRCF.setZero(3,_body->nrJ()*4);
  for(int i=0; i<_body->nrJ(); i++)
    if(_body->joint(i)._M>0) {
      ROTI(_JRCF,i)=-g*_body->joint(i)._MC.transpose().template cast<T>();
      CTRI(_JRCF,i)=-_body->joint(i)._M*g;
    }
}
void ConvHullPBDSimulator::detectCurrentContact() {
  _manifolds.clear();
  detectContact(_pos._TM);
}
void ConvHullPBDSimulator::simpleCollisionHandle(){
  int nrJ=_body->nrJ();
  CollisionGradInfo<T> grad(*_body,setKinematic(_pos._xM,_t+_dt));
  for(int k1=0; k1<nrJ; k1++) {
    for(int k2=k1+1;k2<nrJ;k2++){
      GJKPolytope<T>p1(k1,*_body,grad);
      GJKPolytope<T>p2(k2,*_body,grad);
      T dist=GJKPolytope<T>::distance(p1,p2,NULL,NULL);
    }
  }
}
void ConvHullPBDSimulator::debugEnergy(T scale) {
  DEFINE_NUMERIC_DELTA_T(T)
  //generate random pose
  GradInfo newPos(*_body,Vec::Random(_body->nrDOF())*scale);
  _pos.reset(*_body,Vec::Random(_body->nrDOF())*scale);
  _lastPos.reset(*_body,Vec::Random(_body->nrDOF())*scale);

  //generate random contact
  _manifolds.clear();
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    ContactPoint p;
    p._ptA=Vec3T::Random().template cast<GEOMETRY_SCALAR>();
    p._ptB=Vec3T::Random().template cast<GEOMETRY_SCALAR>();
    p._nA2B=Vec3T::Random().normalized().template cast<GEOMETRY_SCALAR>();
    if(p.depth()<0)
      p._nA2B*=-1;
    ContactManifold m;
    m._points.push_back(p);
    //add dynamic-static contact
    m._jidA=k;
    m._jidB=-1;
    _manifolds.push_back(m);
    //add static-dynamic contact
    m._jidA=-1;
    m._jidB=k;
    _manifolds.push_back(m);
    //add dynamic-dynamic contact
    m._jidA=rand()%nrJ;
    m._jidB=k;
    if(m._jidA!=m._jidB)
      _manifolds.push_back(m);
  }
  computeLocalContactPos(newPos._TM);

  //generate random drag
  _drags.clear();
  for(int k=0; k<nrJ; k++) {
    DragEnergy d;
    d._jid=k;
    d._k=rand()/(T)RAND_MAX;
    d._pt=Vec3T::Random();
    d._ptL=Vec3T::Random();
    _drags.push_back(d);
  }

  //generate random PD target and joint limit
  Vec P=Vec::Random(_body->nrDOF());
  Vec D=Vec::Random(_body->nrDOF());
  for(int k=0; k<nrJ; k++) {
    Joint& J=_body->joint(k);
    J._limits.row(2).setRandom();

    PhysicsParameter& p=_params[k];
    p._kp=rand()/(T)RAND_MAX;
    p._kd=rand()/(T)RAND_MAX;
    p._tarP=[&](T,int n) {
      return P.segment(J._offDOF,n);
    };
    p._tarD=[&](T,int n) {
      return D.segment(J._offDOF,n);
    };
  }

  //debug DE/DDE
  MatT DDE,DDE2;
  Vec DE,DE2,dx=Vec::Random(_body->nrDOF());
  GradInfo newPos2(*_body,newPos._xM+dx*DELTA);
  CollisionGradInfo<T> grad,grad2;
  grad.reset(*_body,newPos._xM);
  grad2.reset(*_body,newPos2._xM);
  T e=energy(grad,DE,DDE);
  T e2=energy(grad2,DE2,DDE2);
  //std::cout<<DE.cwiseAbs().maxCoeff()<<" "<<DELTA<<std::endl;
  DEBUG_GRADIENT("DE",DE.dot(dx),DE.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("DDE",(DDE*dx).norm(),(DDE*dx-(DE2-DE)/DELTA).norm())

  //matrix solver
  MatT HInvH;
  setJTJ(false);
  setCrossTerm(false);
  energy(grad,DE,DDE,true);
  PBDMatrixSolverEigen solEigen(_body);
  solEigen.compute(DDE);
  HInvH=solEigen.solve(DDE);
  DEBUG_GRADIENT("HInvH-Eigen",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
  PBDMatrixSolverCRBA solCRBA(_body);
  solCRBA.compute(DDE);
  HInvH=solCRBA.solve(DDE);
  DEBUG_GRADIENT("HInvH-CRBA",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
  PBDMatrixSolverABA solABA(_body);
  solABA.compute(Vec(newPos._xM),_MRR,_MRt,_MtR,_Mtt,_diag);
  HInvH=solABA.solve(DDE);
  DEBUG_GRADIENT("HInvH-ABA",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
}
void ConvHullPBDSimulator::setOutput(bool output) {
  _output=output;
}
void ConvHullPBDSimulator::setJTJ(bool JTJ) {
  _JTJ=JTJ;
}
void ConvHullPBDSimulator::setCrossTerm(bool cross) {
  _crossTerm=cross;
  if(_crossTerm)
    _sol.reset(new PBDMatrixSolverEigen(_body));
  else if(USE_CRBA)
    _sol.reset(new PBDMatrixSolverCRBA(_body));
  else {
    _sol.reset(new PBDMatrixSolverABA(_body));
    _JTJ=true;
  }
}
//helper
void ConvHullPBDSimulator::update(const GradInfo& newPos,GradInfo& newPos2,Vec& D,const Vec& DE,const MatT& DDE,T alpha) const {
  MatT DDER;//=DDE;
  DDER.diagonal().array()+=_alpha;
  //compute
  if(std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol))
    std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol)->compute(Vec(newPos._xM),_MRR,_MRt,_MtR,_Mtt,Vec(_diag+Vec::Constant(newPos._xM.size(),_alpha)));
  else _sol->compute(DDER);
  newPos2.reset(*_body,newPos._xM-_sol->solve(MatT(DE)));
  //D=-_sol->solve(MatT(DE));
  //newPos2.reset(*_body,newPos._xM-(1.0/_alpha)*DE);
}
void ConvHullPBDSimulator::computeLocalContactPos(const Mat3XT& t) {
  Simulator::computeLocalContactPos(t);
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      p._tangentBound=0;
      if(m._jidA>=0)
        p._ptALast=ROTI(_pos._TM,m._jidA).template cast<GEOMETRY_SCALAR>()*p._ptAL+CTRI(_pos._TM,m._jidA).template cast<GEOMETRY_SCALAR>();
      else p._ptALast=p._ptA;
      if(m._jidB>=0)
        p._ptBLast=ROTI(_pos._TM,m._jidB).template cast<GEOMETRY_SCALAR>()*p._ptBL+CTRI(_pos._TM,m._jidB).template cast<GEOMETRY_SCALAR>();
      else p._ptBLast=p._ptB;
    }
}
void ConvHullPBDSimulator::mask(Vec& diag,Vec& DE,MatT& DDE) const {
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    int nrDJ=J.nrDOF();
    if(_params[k]._isKinematic) {
      diag.segment(J._offDOF,nrDJ).setConstant(std::numeric_limits<double>::infinity());
      DE.segment(J._offDOF,nrDJ).setZero();
      DDE.block(J._offDOF,0,nrDJ,DDE.cols()).setZero();
      DDE.block(0,J._offDOF,DDE.cols(),nrDJ).setZero();
      DDE.diagonal().segment(J._offDOF,nrDJ).setOnes();
    }
  }
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::energy(CollisionGradInfo<T>& grad,Vec& DE,MatT& DDE,bool updateTangentBound) {
  const GradInfo newPos=grad._info;
  Vec tmp;
  Mat3X4T A;
  Vec3T P,ptA;
  Mat3T PPT,ptARC,H;
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  T coef=1.0/(_dt*_dt),E=0,_E=0,val;
  int nrJ=_body->nrJ();
  int nrD=_body->nrDOF();
  DE.setZero(nrD);
  G.setZero(3,4*nrJ);
  Vec _DE;
  MatT _DDE;
  _DE.setZero(nrD);
  _DDE.setZero(nrD,nrD);
  _MRR.setZero(3,3*nrJ);
  _MRt.setZero(3,3*nrJ);
  _MtR.setZero(3,3*nrJ);
  _Mtt.setZero(3,3*nrJ);
  _diag.setZero(nrD);
  DDE.setZero(nrD,nrD);
  //contact
  //GJKPolytope<T> m1,m2;
  std::shared_ptr<MeshExact> mesh;

  _barrier._x0=0.1;
  _ispenetrated=false;
  //bool hasCollisionEnergy=false;
  //std::cout<<"E: "<<E<<std::endl;
  //bool over=false;
  T DIST=10.0;
  for(int k1=0; k1<nrJ; k1++) {
    //if(over==true) break;
    GJKPolytope<T> m1(k1,*_body,grad);
    /*for(int k2=k1+1;k2<nrJ;k2++){
      //if(over==true) break;
      GJKPolytope<T> m2(k2,*_body,grad);
      T dist=GJKPolytope<T>::distance(m1,m2,NULL,NULL);
      if(dist<=_barrier._x0 && dist>_d0){
        CCBarrierConvexEnergy<T,Px> cc(m1,m2,_barrier,_d0,&grad,5e-4);
        
        if(!cc.eval(&_E,_body.get(),&grad,&DE,&DDE)){
          _ispenetrated=true;
        }
        E+=_E;
        //std::cout<<"distself: "<<dist<<" e: "<<E<<std::endl;
        //over=true;
      }
    }*/
    for(auto& m:_shapes){
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,3,1>> iss;
      m->getMesh(vss,iss);
      mesh.reset(new MeshExact(vss,iss,false));
      GJKPolytope<T> m2(mesh);
      T dist=GJKPolytope<T>::distance(m1,m2,NULL,NULL);
      DIST=std::min(DIST,dist);
      if(dist<=2*_barrier._x0){
        //std::cout<<k1<<" "<<dist<<std::endl;
        CCBarrierConvexEnergy<T,Px> cc(m1,m2,_barrier,_d0,&grad,5e-6,false);
        //T Q=_E;
        cc.initialize(NULL,_body.get());
        if(!cc.eval(&_E,_body.get(),&grad,&DE,&DDE)){
          _ispenetrated=true;
        }
        
        E+=_E;
        //over=true;
      }
    }
  }
  std::cout<<DIST<<std::endl;
  DE.setZero(_body->nrDOF());
  grad._info.DTG(*_body,mapM(grad._DTG),mapV(DE));


  /*DDE=grad._HTheta;
  grad._info.toolB(*_body,mapM(grad._DTG),[&](int r,int c,T val) {
    (DDE)(r,c)+=val;
  });*/
    
  for(int k=0; k<nrJ; k++) {
    //dynamic
    const Joint& J=_body->joint(k);
    _MRR.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(newPos._TM,k)*J._MCCT.template cast<T>()*ROTI(newPos._TM,k).transpose())*coef;
    _MRt.template block<3,3>(0,k*3)+=cross<T>(ROTI(newPos._TM,k)*J._MC.template cast<T>())*coef;
    _MtR.template block<3,3>(0,k*3)-=cross<T>(ROTI(newPos._TM,k)*J._MC.template cast<T>())*coef;
    _Mtt.template block<3,3>(0,k*3)+=Mat3T::Identity()*J._M*coef;
    PPT=J._MCCT.template cast<T>();
    P=J._MC.template cast<T>();
    nrD=J.nrDOF();
    //kinematic force
    A=TRANSI(newPos._TM,k)-2*TRANSI(_pos._TM,k)+TRANSI(_lastPos._TM,k);
    E+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace()*coef/2;
    ROTI(G,k)+=(ROT(A)*PPT+CTR(A)*P.transpose())*coef;
    CTRI(G,k)+=(CTR(A)*J._M+ROT(A)*P)*coef;
    //external force
    E+=(TRANSI(newPos._TM,k)*TRANSI(_JRCF,k).transpose()).trace();
    TRANSI(G,k)+=TRANSI(_JRCF,k);
    //P controller
    if(_params[k]._kp>0) {
      tmp=newPos._xM.segment(J._offDOF,nrD)-_params[k]._tarP(_t,nrD);
      E+=tmp.squaredNorm()*_params[k]._kp/2;
      DE.segment(J._offDOF,nrD)+=tmp*_params[k]._kp;
      _diag.segment(J._offDOF,nrD).array()+=_params[k]._kp;
    }
    //D controller
    if(_params[k]._kd>0) {
      tmp=(newPos._xM-_pos._xM).segment(J._offDOF,nrD)/_dt-_params[k]._tarD(_t,nrD);
      E+=tmp.squaredNorm()*_params[k]._kd/2;
      DE.segment(J._offDOF,nrD)+=tmp*_params[k]._kd/_dt;
      _diag.segment(J._offDOF,nrD).array()+=_params[k]._kd/_dt/_dt;
    }
    //joint limit
    for(int c=0; c<J._limits.cols(); c++)
      if(isfinite(J._limits(2,c)) && J._limits(2,c)>0) {
        if(isfinite(J._limits(0,c)) && newPos._xM[J._offDOF+c]<J._limits(0,c)) {
          //lower limited
          val=newPos._xM[J._offDOF+c]-J._limits(0,c);
          E+=val*val*J._limits(2,c)/2;
          DE[J._offDOF+c]+=val*J._limits(2,c);
          _diag[J._offDOF+c]+=J._limits(2,c);
        } else if(isfinite(J._limits(1,c)) && newPos._xM[J._offDOF+c]>J._limits(1,c)) {
          //upper limited
          val=newPos._xM[J._offDOF+c]-J._limits(1,c);
          E+=val*val*J._limits(2,c)/2;
          DE[J._offDOF+c]+=val*J._limits(2,c);
          _diag[J._offDOF+c]+=J._limits(2,c);
        }
      }
  }
  
  
  //drags
  for(const auto& d:_drags) {
    ptA=ROTI(newPos._TM,d._jid)*d._ptL+CTRI(newPos._TM,d._jid);
    TRANSI(G,d._jid)+=d._k*(ptA-d._pt)*Vec4T(d._ptL[0],d._ptL[1],d._ptL[2],1).transpose();
    ptARC=cross<T>(ROTI(newPos._TM,d._jid)*d._ptL);
    _MRR.template block<3,3>(0,d._jid*3)+=ptARC*ptARC.transpose()*d._k;
    _MRt.template block<3,3>(0,d._jid*3)+=H=ptARC*d._k;
    _MtR.template block<3,3>(0,d._jid*3)+=H.transpose();
    _Mtt.template block<3,3>(0,d._jid*3)+=Mat3T::Identity()*d._k;
    E+=(ptA-d._pt).squaredNorm()*d._k/2;
  }
  
  //gradient
  newPos.DTG(*_body,mapM(GB=G),mapV(DE));
  //hessian
  if(_JTJ) {
    newPos.toolA(*_body,newPos,mapM(MRR=_MRR),mapM(MRt=_MRt),mapM(MtR=_MtR),mapM(Mtt=_Mtt),[&](int r,int c,T val) {
      DDE(r,c)+=val;
    });
  } else {
    newPos.toolAB(*_body,mapM(MRR=_MRR),mapM(MRt=_MRt),mapM(MtR=_MtR),mapM(Mtt=_Mtt),mapM(GB=G),[&](int r,int c,T val) {
      DDE(r,c)+=val;
    });
  }
  DDE.diagonal()+=_diag;

  return E;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::contactEnergy
(const ContactManifold& m,ContactPoint& p,
 const GradInfo& newPos,Vec& DE,MatT& DDE,Mat3XT& G,
 Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool updateTangentBound) const {
  T E=0,val;
  Mat3T ptARC,ptBRC,HRR,HRt,HtR,H;
  //energy = kc/2*|max(p.depth(),0)|^2
  H.setZero();
  p._fA.setZero();
  p._fB.setZero();
  //compute energy/gradient/hessian: normal
  E+=normalEnergy(newPos,m,p,G,H);
  //compute energy/gradient/hessian: friction
  E+=tangentEnergy(newPos,m,p,G,H,updateTangentBound);
  //fill-in non-zero gradient/hessian
  if(!H.isZero()) {
    if(m._jidA>=0) {
      ptARC=cross<T>(ROTI(newPos._TM,m._jidA)*p._ptAL.template cast<T>());
      MRR.template block<3,3>(0,m._jidA*3)+=ptARC*H*ptARC.transpose();
      MRt.template block<3,3>(0,m._jidA*3)+=HRt=ptARC*H;
      MtR.template block<3,3>(0,m._jidA*3)+=HRt.transpose();
      Mtt.template block<3,3>(0,m._jidA*3)+=H;
    }
    if(m._jidB>=0) {
      ptBRC=cross<T>(ROTI(newPos._TM,m._jidB)*p._ptBL.template cast<T>());
      MRR.template block<3,3>(0,m._jidB*3)+=ptBRC*H*ptBRC.transpose();
      MtR.template block<3,3>(0,m._jidB*3)+=HtR=H*ptBRC.transpose();
      MRt.template block<3,3>(0,m._jidB*3)+=HtR.transpose();
      Mtt.template block<3,3>(0,m._jidB*3)+=H;
    }
    if(m._jidA>=0 && m._jidB>=0 && _crossTerm) {
      HRR=ptARC*H*ptBRC.transpose();
      newPos.JRCSparse(*_body,m._jidA,[&](int r,const Vec3T& JRA) {
        newPos.JRCSparse(*_body,m._jidB,[&](int c,const Vec3T& JRB) {
          DDE(r,c)-=val=JRA.dot(HRR*JRB);
          DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          DDE(r,c)-=val=JRA.dot(HRt*JtB);
          DDE(c,r)-=val;
        });
      },[&](int r,const Vec3T& JtA) {
        newPos.JRCSparse(*_body,m._jidB,[&](int c,const Vec3T& JRB) {
          DDE(r,c)-=val=JtA.dot(HtR*JRB);
          DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          DDE(r,c)-=val=JtA.dot(H*JtB);
          DDE(c,r)-=val;
        });
      });
    }
  }
  return E;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::normalEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H) const {
  Vec3T nA2B=p._nA2B.template cast<T>();
  Vec3T ptA=p._ptA.template cast<T>();
  Vec3T ptB=p._ptB.template cast<T>();
  if(m._jidA>=0)
    ptA=ROTI(newPos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(newPos._TM,m._jidA);
  if(m._jidB>=0)
    ptB=ROTI(newPos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(newPos._TM,m._jidB);
  T kc=std::max<T>(m._jidA>=0?_params[m._jidA]._kc:0,m._jidB>=0?_params[m._jidB]._kc:0);
  T depth=std::max<T>((ptA-ptB).dot(nA2B),0);
  if(depth>0 && kc>0) {
    H+=nA2B*nA2B.transpose()*kc;
    if(m._jidA>=0) {
      TRANSI(G,m._jidA)+=kc*depth*nA2B*Vec4T((T)p._ptAL[0],(T)p._ptAL[1],(T)p._ptAL[2],1).transpose();
      p._fA-=(kc*depth*nA2B).template cast<GEOMETRY_SCALAR>();
    }
    if(m._jidB>=0) {
      TRANSI(G,m._jidB)-=kc*depth*nA2B*Vec4T((T)p._ptBL[0],(T)p._ptBL[1],(T)p._ptBL[2],1).transpose();
      p._fB+=(kc*depth*nA2B).template cast<GEOMETRY_SCALAR>();
    }
    return kc/2*depth*depth;
  } else return 0;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::tangentEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H,bool updateTangentBound) const {
  if(updateTangentBound) {
    T fri=std::max<T>(m._jidA>=0?_params[m._jidA]._friction:0,m._jidB>=0?_params[m._jidB]._friction:0);
    p._tangentBound=(GEOMETRY_SCALAR)(fri*std::max<T>(p._fA.template cast<T>().norm(),p._fB.template cast<T>().norm()));
  }
  if(p._tangentBound>0) {
    Mat3X2T tA2B=p._tA2B.template cast<T>();
    Vec3T ptA=p._ptA.template cast<T>();
    Vec3T ptB=p._ptB.template cast<T>();
    Vec3T ptALast=p._ptALast.template cast<T>();
    Vec3T ptBLast=p._ptBLast.template cast<T>();
    if(m._jidA>=0)
      ptA=ROTI(newPos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(newPos._TM,m._jidA);
    if(m._jidB>=0)
      ptB=ROTI(newPos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(newPos._TM,m._jidB);
    //tangent velocity
    Vec2T dPT=tA2B.transpose()*((ptA-ptALast)-(ptB-ptBLast));
    T dPTLen=dPT.norm(),dPTLen2=dPTLen*dPTLen,dPTLen3=dPTLen2*dPTLen;
    T epsVDt=_epsV*_dt,invEpsVDt=1/epsVDt,invEpsVDt2=invEpsVDt*invEpsVDt,f0,tb=(T)p._tangentBound;
    Vec3T f1,dPTP=tA2B*dPT;
    if(dPTLen<epsVDt) {
      f0=-dPTLen3*invEpsVDt2/3+dPTLen2*invEpsVDt+epsVDt/3;
      //gradient
      f1=(-dPTLen*invEpsVDt2+2*invEpsVDt)*dPTP;
      if(m._jidA>=0) {
        TRANSI(G,m._jidA)+=f1*Vec4T((T)p._ptAL[0],(T)p._ptAL[1],(T)p._ptAL[2],1).transpose()*tb;
        p._fA-=f1.template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
      }
      if(m._jidB>=0) {
        TRANSI(G,m._jidB)-=f1*Vec4T((T)p._ptBL[0],(T)p._ptBL[1],(T)p._ptBL[2],1).transpose()*tb;
        p._fB+=f1.template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
      }
      //hessian
      H+=(-dPTLen*invEpsVDt2+2*invEpsVDt)*tA2B*tA2B.transpose()*tb;
      if(!_JTJ)
        H-=invEpsVDt2*dPTP*dPTP.transpose()/std::max<T>(dPTP.norm(),Epsilon<T>::defaultEps())*tb;
    } else {
      f0=dPTLen;
      //gradient
      f1=dPTP/f0;
      if(m._jidA>=0) {
        TRANSI(G,m._jidA)+=f1*Vec4T((T)p._ptAL[0],(T)p._ptAL[1],(T)p._ptAL[2],1).transpose()*tb;
        p._fA-=f1.template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
      }
      if(m._jidB>=0) {
        TRANSI(G,m._jidB)-=f1*Vec4T((T)p._ptBL[0],(T)p._ptBL[1],(T)p._ptBL[2],1).transpose()*tb;
        p._fB+=f1.template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
      }
      //hessian
      H+=tA2B*tA2B.transpose()/f0*tb;
      if(!_JTJ)
        H-=dPTP*dPTP.transpose()/(f0*f0*f0)*tb;
    }
    return f0*tb;
  } else return 0;
}
}
//ATBBFTWcVCbNWZdRmAH8BxchGeqgE185B444