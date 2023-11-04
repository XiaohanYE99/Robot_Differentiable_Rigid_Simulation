#include "ConvHullPBDSimulator.h"
#include "PBDMatrixSolver.h"
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
#define USE_CRBA 0
ConvHullPBDSimulator::ConvHullPBDSimulator(T dt):Simulator(dt),_gTol(1e-4f),_alpha(1e-6f),_epsV(1e-1f),_output(true),_JTJ(true),_crossTerm(false),_maxIt(1e4) {
  _barrier._x0=0.1;
}
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
  //update kinematic state
  GradInfo newPos(*_body,setKinematic(_pos._xM,_t+_dt)),newPos2;
  //normal solve
  Vec D,DE,DE2;
  MatT DDE,DDE2;
  T alphaMax=1e8f,alphaMin=1e-8f;
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
    if((e2<e || DE2.cwiseAbs().maxCoeff()<_gTol) && (isfinite(e2))) { //
      _alpha=std::max<T>(_alpha*.8f,alphaMin);
      newPos=newPos2;
      e=e2;
      DE=DE2;
      DDE=DDE2;
      iter++;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e  << " alpha=" << _alpha << " gNorm=" << DE.cwiseAbs().maxCoeff()<<std::endl;
      //termination: gradient tolerance
      if(DE2.cwiseAbs().maxCoeff()<_gTol)
        break;
    } else {
      _alpha=std::min<T>(_alpha*1.5f,alphaMax);
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e <<  " alpha=" << _alpha << " gNorm=" << DE2.cwiseAbs().maxCoeff()<<std::endl;
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
void ConvHullPBDSimulator::addShape(std::shared_ptr<ShapeExact> shape) {
  _shapes.push_back(shape);
  std::shared_ptr<MeshExact> mesh;
  for(auto& m:_shapes) {
    std::vector<Eigen::Matrix<T,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    m->getMesh(vss,iss);
    mesh.reset(new MeshExact(vss,iss,false));
    _obs.push_back(std::shared_ptr<GJKPolytope<T>>(new GJKPolytope<T>(mesh)));
    auto m2=_obs.back();
    //std::cout << "- " << m2->getBB().minCorner().transpose() << " " << m2->getBB().maxCorner().transpose() << std::endl;
  }
}
void ConvHullPBDSimulator::detectCurrentContact(){

}
void ConvHullPBDSimulator::detectCCDContact(const Mat3XT& t){
  _manifolds.clear();
  detectContact(t);
}
void ConvHullPBDSimulator::detectContact(const Mat3XT& t) {
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  T x0=(2*_barrier._x0)/(1-_barrier._x0);
  _contact->generateManifolds(_manifolds,t.template cast<GEOMETRY_SCALAR>(),x0);
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

void ConvHullPBDSimulator::debugEnergy(T scale) {
  DEFINE_NUMERIC_DELTA_T(T)
  //generate random pose
  GradInfo newPos(*_body,Vec::Random(_body->nrDOF())*scale);
  _pos.reset(*_body,Vec::Random(_body->nrDOF())*scale);
  _lastPos.reset(*_body,Vec::Random(_body->nrDOF())*scale);


  int nrJ=_body->nrJ();

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
  /*MatT HInvH;
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
  DEBUG_GRADIENT("HInvH-ABA",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())*/
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
  //debugBVHenergy(grad,updateTangentBound);
  const GradInfo newPos=grad._info;
  detectCCDContact(newPos._TM);
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
  _MRR.setZero(3,3*nrJ);
  _MRt.setZero(3,3*nrJ);
  _MtR.setZero(3,3*nrJ);
  _Mtt.setZero(3,3*nrJ);
  _diag.setZero(nrD);
  DDE.setZero(nrD,nrD);
  //contact
  //GJKPolytope<T> m1,m2;
  T coefBarrier=1e-2;
  std::shared_ptr<MeshExact> m1,m2;
  std::cout<<_manifolds.size()<<std::endl;
  for(auto& m:_manifolds){
    //std::cout<<m._jidA<<"  "<<m._jidB<<std::endl;
    std::vector<Eigen::Matrix<double,3,1>> vss1,vss2;
    std::vector<Eigen::Matrix<int,3,1>> iss1,iss2;
    if(m._sA){
      if(m._jidA<0){
        m._sA->getMesh(vss1,iss1);
        m1.reset(new MeshExact(vss1,iss1,false));
        m._pA.reset(new GJKPolytope<T>(m1));
      }
      else m._pA.reset(new GJKPolytope<T>(m._jidA,*_body,grad));
    }
    if(m._sB){
      if(m._jidB<0){
        m._sB->getMesh(vss2,iss2);
        m2.reset(new MeshExact(vss2,iss2,false));
        m._pB.reset(new GJKPolytope<T>(m2));
      }
      else m._pB.reset(new GJKPolytope<T>(m._jidB,*_body,grad));
    }
    if(!m._pA || !m._pB) continue;
    CCBarrierConvexEnergy<T,CLogx> cc(*m._pA,*m._pB,_barrier,0,&grad,coefBarrier);
    //cc.initialize(NULL,_body.get());
    if(!cc.eval(&_E,_body.get(),&grad,NULL,NULL))
      return std::numeric_limits<T>::infinity();
    E+=_E;
  }

  /*for(int k1=0; k1<nrJ; k1++) {
    std::shared_ptr<GJKPolytope<T>> m1=grad._polytopes[k1];
    if(!m1) continue;
    //std::cout << k1 << " " << m1->getBB().minCorner().transpose() << " " << m1->getBB().maxCorner().transpose();
    //m1->writeVTK("shape"+std::to_string(k1)+".vtk",NULL);
    for(int k2=k1+1; k2<nrJ; k2++) {
      if(_body->joint(k2)._parent==-1) continue;
      if(k2==_body->joint(k1)._parent || k1==_body->joint(k2)._parent) continue;
      std::shared_ptr<GJKPolytope<T>> m2=grad._polytopes[k2];
      //T dist=GJKPolytope<T>::distance(*m1,*m2,NULL,NULL);
      if(!m2) continue;
      CCBarrierConvexEnergy<T,CLogx> cc(*m1,*m2,_barrier,0,&grad,coefBarrier);
      //cc.initialize(NULL,_body.get());
      if(!cc.eval(&_E,_body.get(),&grad,NULL,NULL))
        return std::numeric_limits<T>::infinity();
      E+=_E;
    }
    for(auto& m2:_obs) {
      //T dist=GJKPolytope<T>::distance(*m1,*m2,NULL,NULL);
      //std::cout << " " << dist;
      CCBarrierConvexEnergy<T,CLogx> cc(*m1,*m2,_barrier,0,&grad,coefBarrier);
      //cc.initialize(NULL,_body.get());
      if(!cc.eval(&_E,_body.get(),&grad,NULL,NULL))
        return std::numeric_limits<T>::infinity();
      E+=_E;
    }
  }*/
  Mat3XT tmpDTG;
  DE.setZero(_body->nrDOF());
  grad._info.DTG(*_body,mapM(tmpDTG=grad._DTG),mapV(DE));

  DDE=grad._HTheta;
  grad._info.toolB(*_body,mapM(tmpDTG=grad._DTG),[&](int r,int c,T val) {
    (DDE)(r,c)+=val;
  });

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
void ConvHullPBDSimulator::debugBVHenergy(CollisionGradInfo<T>& grad,bool updateTangentBound){
  const GradInfo newPos=grad._info;
  CollisionGradInfo<T>grad1;
  grad1.reset(*_body,newPos._xM);
  Vec DE,DE1;
  MatT DDE,DDE1;
  detectCCDContact(newPos._TM);
  T E1=0,E=0,_E=0;
  int nrJ=_body->nrJ();
  int nrD=_body->nrDOF();
  DE.setZero(nrD);
  DDE.setZero(nrD,nrD);
  DE1.setZero(nrD);
  DDE1.setZero(nrD,nrD);
  //contact
  //GJKPolytope<T> m1,m2;
  T coefBarrier=1e-2;
  std::shared_ptr<MeshExact> m1,m2;
  for(auto& m:_manifolds){
    //std::cout<<m._jidA<<"  "<<m._jidB<<std::endl;
    std::vector<Eigen::Matrix<double,3,1>> vss1,vss2;
    std::vector<Eigen::Matrix<int,3,1>> iss1,iss2;
    if(m._sA){
      if(m._jidA<0){
        m._sA->getMesh(vss1,iss1);
        m1.reset(new MeshExact(vss1,iss1,false));
        m._pA.reset(new GJKPolytope<T>(m1));
      }
      else m._pA.reset(new GJKPolytope<T>(m._jidA,*_body,grad1));
    }
    if(m._sB){
      if(m._jidB<0){
        m._sB->getMesh(vss2,iss2);
        m2.reset(new MeshExact(vss2,iss2,false));
        m._pB.reset(new GJKPolytope<T>(m2));
      }
      else m._pB.reset(new GJKPolytope<T>(m._jidB,*_body,grad1));
    }
    if(!m._pA || !m._pB) continue;
    CCBarrierConvexEnergy<T,CLogx> cc(*m._pA,*m._pB,_barrier,0,&grad1,coefBarrier);
    //cc.initialize(NULL,_body.get());
    cc.eval(&_E,_body.get(),&grad1,NULL,NULL);
    E1+=_E;
  }
  Mat3XT tmpDTG1;
  DE1.setZero(_body->nrDOF());
  grad1._info.DTG(*_body,mapM(tmpDTG1=grad1._DTG),mapV(DE1));

  DDE1=grad1._HTheta;
  grad1._info.toolB(*_body,mapM(tmpDTG1=grad1._DTG),[&](int r,int c,T val) {
    (DDE1)(r,c)+=val;
  });

  for(int k1=0; k1<nrJ; k1++) {
    std::shared_ptr<GJKPolytope<T>> m1=grad._polytopes[k1];
    if(!m1) continue;
    //std::cout << k1 << " " << m1->getBB().minCorner().transpose() << " " << m1->getBB().maxCorner().transpose();
    //m1->writeVTK("shape"+std::to_string(k1)+".vtk",NULL);
    for(int k2=k1+1; k2<nrJ; k2++) {
      if(_body->joint(k2)._parent==-1) continue;
      if(k2==_body->joint(k1)._parent || k1==_body->joint(k2)._parent) continue;
      std::shared_ptr<GJKPolytope<T>> m2=grad._polytopes[k2];
      //T dist=GJKPolytope<T>::distance(*m1,*m2,NULL,NULL);
      if(!m2) continue;
      CCBarrierConvexEnergy<T,CLogx> cc(*m1,*m2,_barrier,0,&grad,coefBarrier);
      //cc.initialize(NULL,_body.get());
      cc.eval(&_E,_body.get(),&grad,NULL,NULL);
      E+=_E;
    }
    for(auto& m2:_obs) {
      //T dist=GJKPolytope<T>::distance(*m1,*m2,NULL,NULL);
      //std::cout << " " << dist;
      CCBarrierConvexEnergy<T,CLogx> cc(*m1,*m2,_barrier,0,&grad,coefBarrier);
      //cc.initialize(NULL,_body.get());
      !cc.eval(&_E,_body.get(),&grad,NULL,NULL);
      E+=_E;
    }
  }
  Mat3XT tmpDTG;
  DE.setZero(_body->nrDOF());
  grad._info.DTG(*_body,mapM(tmpDTG=grad._DTG),mapV(DE));

  DDE=grad._HTheta;
  grad._info.toolB(*_body,mapM(tmpDTG=grad._DTG),[&](int r,int c,T val) {
    (DDE)(r,c)+=val;
  });
  std::cout<<"E1: "<<E1<<"   "<<"error: "<<E1-E<<std::endl;
  std::cout<<"DE1: "<<DE.cwiseAbs().maxCoeff()<<"error: "<<(DE1-DE).cwiseAbs().maxCoeff()<<std::endl;
  std::cout<<"DDE1: "<<DDE.norm()<<"   "<<"error: "<<(DDE1-DDE).norm()<<std::endl;
}
}
//ATBBnHDkBNsAuRvr4EUGTmyZCytYB0E4A4FC
//git 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
