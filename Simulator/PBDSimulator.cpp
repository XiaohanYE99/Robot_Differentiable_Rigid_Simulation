#include "PBDSimulator.h"
#include "PBDMatrixSolver.h"
#include "JointLimit.h"
#include "SoftJoint.h"
#include <Utils/RotationUtils.h>

namespace PHYSICSMOTION {
//SolverState
PBDSimulator::SolverState::SolverState() {}
PBDSimulator::SolverState::SolverState(const ArticulatedBody& body,const Vec& x) {
  reset(body,x);
}
void PBDSimulator::SolverState::reset(const ArticulatedBody& body,const Vec& x) {
  _pos.reset(body,x);
}
//PBDSimulator
PBDSimulator::PBDSimulator(T dt):Simulator(dt),_gTol(1e-7f),_epsV(1e-1f),_JTJ(true),_crossTerm(false),_maxIt(1e4) {}
void PBDSimulator::setArticulatedBody(std::shared_ptr<ArticulatedBody> body) {
  Simulator::setArticulatedBody(body);
  setGravity(Vec3T::Zero());
  _pos.reset(*body,Vec::Zero(body->nrDOF()));
  _lastPos.reset(*body,Vec::Zero(body->nrDOF()));
  setCrossTerm(_crossTerm);
}
void PBDSimulator::step() {
  if(_body->nrDOF()==0)
    return;
  //generate contacts
  detectCurrentContact();
  //update kinematic state
  SolverState state(*_body,setKinematic(_pos._xM,_t+_dt)),state2;
  //integrate all constraints
  solveBody(state,state2);
  //update
  _lastPos=_pos;
  _pos=state._pos;
  _t+=_dt;
}
PBDSimulator::Vec PBDSimulator::pos() const {
  return _pos._xM;
}
void PBDSimulator::setPos(const Vec& pos) {
  _pos.reset(*_body,setKinematic(pos,_t));
}
PBDSimulator::Vec PBDSimulator::vel() const {
  return (_pos._xM-_lastPos._xM)/_dt;
}
void PBDSimulator::setVel(const Vec& vel) {
  _lastPos.reset(*_body,setKinematic(_pos._xM-vel*_dt,_t-_dt));
}
void PBDSimulator::detectCurrentContact() {
  _manifolds.clear();
  detectContact(_pos._TM);
}
void PBDSimulator::debugEnergy(T scale) {
  DEFINE_NUMERIC_DELTA_T(T)
  //generate random pose
  SolverState state(*_body,Vec::Random(_body->nrDOF())*scale);
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
  computeLocalContactPos(state._pos._TM);

  //generate random drag
  _joints.clear();
  for(int k=0; k<nrJ; k++) {
    SoftJoint j;
    j.setRandom(*_body);
    _joints.push_back(j);
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
  SolverState state2(*_body,state._pos._xM+dx*DELTA);
  T e=energy(state,true);
  T e2=energy(state2);
  DEBUG_GRADIENT("DE",state._DE.dot(dx),state._DE.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("DDE",(state._DDE*dx).norm(),(state._DDE*dx-(state2._DE-state._DE)/DELTA).norm())

  //matrix solver
  MatT HInvH;
  setJTJ(false);
  setCrossTerm(false);
  energy(state,true);
  PBDMatrixSolverEigen solEigen(_body);
  solEigen.compute(state._DDE);
  HInvH=solEigen.solve(state._DDE);
  DEBUG_GRADIENT("HInvH-Eigen",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
  PBDMatrixSolverCRBA solCRBA(_body);
  solCRBA.compute(state._DDE);
  HInvH=solCRBA.solve(state._DDE);
  DEBUG_GRADIENT("HInvH-CRBA",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
  PBDMatrixSolverABA solABA(_body);
  solABA.compute(Vec(state._pos._xM),state._MRR,state._MRt,state._MtR,state._Mtt,state._diag);
  HInvH=solABA.solve(state._DDE);
  DEBUG_GRADIENT("HInvH-ABA",HInvH.norm(),(HInvH-MatT::Identity(DDE.rows(),DDE.cols())).norm())
}
void PBDSimulator::setJTJ(bool JTJ) {
  _JTJ=JTJ;
}
void PBDSimulator::setCrossTerm(bool cross,bool CRBA) {
  _crossTerm=cross;
  if(_crossTerm) {
    _sol.reset(new PBDMatrixSolverEigen(_body));
    if(_output)
      std::cout << "Choosing Eigen's native matrix solver!" << std::endl;
  } else if(CRBA) {
    _sol.reset(new PBDMatrixSolverCRBA(_body));
    if(_output)
      std::cout << "Choosing CRBA matrix solver!" << std::endl;
  } else {
    _sol.reset(new PBDMatrixSolverABA(_body));
    if(_output)
      std::cout << "Choosing ABA matrix solver!" << std::endl;
    _JTJ=true;
  }
}
//helper
void PBDSimulator::solveBody(SolverState& state,SolverState& state2) {
  T alphaMax=1e6f,alphaMin=1e-6f,alpha=alphaMin;
  T e=energy(state,true),e2=0;
  mask(&(state._diag),&(state._DE),&(state._DDE));
  for(int iter=0; iter<_maxIt;) {
    //update configuration
    update(state,state2,alpha);
    e2=energy(state2);
    mask(&(state2._diag),&(state2._DE),&(state2._DDE));
    //iteration update
    if(e2<e || state2._DE.cwiseAbs().maxCoeff()<_gTol) {
      alpha=std::max<T>(alpha*.5f,alphaMin);
      state=state2;
      e=e2;
      iter++;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " gNorm=" << state2._DE.cwiseAbs().maxCoeff() << " alpha=" << alpha << std::endl;
      //termination: gradient tolerance
      if(state2._DE.cwiseAbs().maxCoeff()<_gTol)
        break;
    } else {
      alpha=std::min<T>(alpha*5.f,alphaMax);
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " E2=" << e2 << " alpha=" << alpha << std::endl;
      //termination: numerical issue
      if(alpha>=alphaMax)
        break;
    }
  }
}
void PBDSimulator::update(const SolverState& state,SolverState& state2,T alpha) {
  //compute
  if(std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol)) {
    _DDER=state._diag;
    _DDER.diagonal().array()+=alpha;
    std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol)->compute(Vec(state._pos._xM),state._MRR,state._MRt,state._MtR,state._Mtt,_DDER);
  } else {
    _DDER=state._DDE;
    _DDER.diagonal().array()+=alpha;
    _sol->compute(_DDER);
  }
  //update
  state2.reset(*_body,state._pos._xM-_sol->solve(MatT(state._DE)));
}
void PBDSimulator::computeLocalContactPos(const Mat3XT& t) {
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
PBDSimulator::T PBDSimulator::energy(SolverState& state,bool updateTangentBound) {
  T E=0,damping;
  int nrJ=_body->nrJ();
  int nrD=_body->nrDOF();
  state._DE.setZero(nrD);
  _G.setZero(3,4*nrJ);
  state._MRR.setZero(3,3*nrJ);
  state._MRt.setZero(3,3*nrJ);
  state._MtR.setZero(3,3*nrJ);
  state._Mtt.setZero(3,3*nrJ);
  state._diag.setZero(nrD,nrD);
  state._DDE.setZero(nrD,nrD);
  for(int k=0; k<nrJ; k++) {
    //inertial
    const auto& J=_body->joint(k);
    nrD=J.nrDOF();
    damping=0;
    if(J._damping.size()>0 && J._damping.maxCoeff()>0)
      damping=J._damping.maxCoeff();
    E+=energyInertial(state._pos,_pos,_lastPos,k,damping,
                      J._M,J._MC.template cast<T>(),J._MCCT.template cast<T>(),
                      &_G,state._MRR,state._MRt,state._MtR,state._Mtt);
    //PD controller
    E+=energyPDController(mapV2CV(state._pos._xM),mapV2CV(_pos._xM),k,J,nrD,&(state._DE),&(state._diag));
    //joint limit
    E+=JointLimit::energy(mapV2CV(state._pos._xM),J,nrD,&(state._DE),&(state._diag));
  }
  //contact
  for(auto& m:_manifolds)
    for(auto& p:m._points)
      E+=contactEnergy(m,p,state,updateTangentBound);
  //joints
  for(const auto& joint:_joints)
    E+=joint.energy(*_body,state._pos,&(state._DE),state._DDE,_G,state._MRR,state._MRt,state._MtR,state._Mtt,_crossTerm);
  //gradient
  state._pos.DTG(*_body,mapM(_GB=_G),mapV(state._DE));
  //hessian
  if(!std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol)) {
    if(_JTJ) {
      state._pos.toolA(*_body,state._pos,mapM(state._MRR),mapM(state._MRt),mapM(state._MtR),mapM(state._Mtt),[&](int r,int c,T val) {
        state._DDE(r,c)+=val;
      });
    } else {
      state._pos.toolAB(*_body,mapM(state._MRR),mapM(state._MRt),mapM(state._MtR),mapM(state._Mtt),mapM(_GB=_G),[&](int r,int c,T val) {
        state._DDE(r,c)+=val;
      });
    }
    for(int k=0; k<nrJ; k++) {
      const Joint& J=_body->joint(k);
      nrD=J.nrDOF();
      state._DDE.template block(J._offDOF,J._offDOF,nrD,nrD)+=state._diag.template block(J._offDOF,J._offDOF,nrD,nrD);
    }
  }
  return E;
}
PBDSimulator::T PBDSimulator::contactEnergy(const ContactManifold& m,ContactPoint& p,SolverState& state,bool updateTangentBound) {
  T E=0,val;
  Mat3T ptARC,ptBRC,HRR,HRt,HtR,H;
  //energy = kc/2*|max(p.depth(),0)|^2
  H.setZero();
  p._fA.setZero();
  p._fB.setZero();
  //compute energy/gradient/hessian: normal
  E+=normalEnergy(state._pos,m,p,_G,H);
  //compute energy/gradient/hessian: friction
  E+=tangentEnergy(state._pos,m,p,_G,H,updateTangentBound);
  //fill-in non-zero gradient/hessian
  if(!H.isZero()) {
    if(m._jidA>=0) {
      ptARC=cross<T>(ROTI(state._pos._TM,m._jidA)*p._ptAL.template cast<T>());
      state._MRR.template block<3,3>(0,m._jidA*3)+=ptARC*H*ptARC.transpose();
      state._MRt.template block<3,3>(0,m._jidA*3)+=HRt=ptARC*H;
      state._MtR.template block<3,3>(0,m._jidA*3)+=HRt.transpose();
      state._Mtt.template block<3,3>(0,m._jidA*3)+=H;
    }
    if(m._jidB>=0) {
      ptBRC=cross<T>(ROTI(state._pos._TM,m._jidB)*p._ptBL.template cast<T>());
      state._MRR.template block<3,3>(0,m._jidB*3)+=ptBRC*H*ptBRC.transpose();
      state._MtR.template block<3,3>(0,m._jidB*3)+=HtR=H*ptBRC.transpose();
      state._MRt.template block<3,3>(0,m._jidB*3)+=HtR.transpose();
      state._Mtt.template block<3,3>(0,m._jidB*3)+=H;
    }
    if(m._jidA>=0 && m._jidB>=0 && _crossTerm) {
      HRR=ptARC*H*ptBRC.transpose();
      state._pos.JRCSparse(*_body,m._jidA,[&](int r,const Vec3T& JRA) {
        state._pos.JRCSparse(*_body,m._jidB,[&](int c,const Vec3T& JRB) {
          state._DDE(r,c)-=val=JRA.dot(HRR*JRB);
          state._DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          state._DDE(r,c)-=val=JRA.dot(HRt*JtB);
          state._DDE(c,r)-=val;
        });
      },[&](int r,const Vec3T& JtA) {
        state._pos.JRCSparse(*_body,m._jidB,[&](int c,const Vec3T& JRB) {
          state._DDE(r,c)-=val=JtA.dot(HtR*JRB);
          state._DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          state._DDE(r,c)-=val=JtA.dot(H*JtB);
          state._DDE(c,r)-=val;
        });
      });
    }
  }
  return E;
}
PBDSimulator::T PBDSimulator::normalEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H) const {
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
PBDSimulator::T PBDSimulator::tangentEnergy(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,Mat3XT& G,Mat3T& H,bool updateTangentBound) const {
  if(updateTangentBound) {
    T fri=std::max<T>(m._jidA>=0?_params[m._jidA]._friction:0,m._jidB>=0?_params[m._jidB]._friction:0);
    p._tangentBound=(GEOMETRY_SCALAR)(fri*std::max<T>(p._fA.template cast<T>().norm(),p._fB.template cast<T>().norm()));
    //std::cout << p._tangentBound << std::endl;
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
