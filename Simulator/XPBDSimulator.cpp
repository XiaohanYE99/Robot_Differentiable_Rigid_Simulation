#include "XPBDSimulator.h"
#include "PBDMatrixSolver.h"
#include "JointLimit.h"
#include "SoftJoint.h"
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
XPBDSimulator::XPBDSimulator(T dt):PBDSimulator(dt) {
  _maxIt=5;
  _gTol=1e-2f;
  _epsDenom=1e-3f;
}
void XPBDSimulator::step() {
  if(_body->nrDOF()==0)
    return;
  //generate contacts
  detectCurrentContact();
  //update kinematic state
  SolverState state(*_body,setKinematic(_pos._xM,_t+_dt)),state2;
  //integrate all constraints
  solveBody(state,state2);
  //project all constraints
  projectBody(state);
  //update
  _lastPos=_pos;
  _pos=state._pos;
  _t+=_dt;
}
void XPBDSimulator::debugEnergy(T scale) {
  DEFINE_NUMERIC_DELTA_T(T)
  //generate random pose
  SolverState state(*_body,Vec::Random(_body->nrDOF())*scale);
  _pos.reset(*_body,Vec::Random(_body->nrDOF())*scale);
  _lastPos.reset(*_body,Vec::Random(_body->nrDOF())*scale);

  //debug DE/DDE
  MatT DDE,DDE2;
  Vec DE,DE2,dx=Vec::Random(_body->nrDOF());
  SolverState state2(*_body,state._pos._xM+dx*DELTA);
  T e=energy(state);
  T e2=energy(state2);
  DEBUG_GRADIENT("DE",state._DE.dot(dx),state._DE.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("DDE",(state._DDE*dx).norm(),(state._DDE*dx-(state2._DE-state._DE)/DELTA).norm())
}
void XPBDSimulator::setCrossTerm(bool cross,bool CRBA) {
  _crossTerm=cross;
  if(_crossTerm) {
    if(_output)
      std::cout << "XPBD does not support cross term!" << std::endl;
    _crossTerm=false;
  }
  if(CRBA) {
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
XPBDSimulator::T XPBDSimulator::energy(SolverState& state,bool updateTangentBound) {
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
  }
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
void XPBDSimulator::projectBody(SolverState& state) {
  _Css.clear();
  int nrD=_body->nrDOF();
  //contact
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      p._fA.setZero();
      p._fB.setZero();
    }
  //position project
  for(int i=0,off=0,off0; i<_maxIt; i++,off=0) {
    //project joint limit
    for(int k=0; k<_body->nrJ(); k++) {
      const Joint& J=_body->joint(k);
      off0=off;
      JointLimit::constraint(mapV2CV(state._pos._xM),J,J.nrDOF(),nrD,_Css,off);
      project(state._pos,off0,off,NULL,false); //only update at last
    }
    state.reset(*_body,state._pos._xM);
    //project contact
    for(auto& m:_manifolds)
      for(auto& p:m._points) {
        T fri=std::max<T>(m._jidA>=0?_params[m._jidA]._friction:0,m._jidB>=0?_params[m._jidB]._friction:0);
        Vec2T fT=Vec2T::Zero();
        T fN=0,dfN=0,lmt=0;
        //normal
        off0=off;
        contactNormal(state._pos,m,p,off);
        project(state._pos,off0,off,NULL,true);
        for(int k=off0; k<off; k++) {
          fN+=_Css[k]._lambda;
          dfN+=_Css[k]._dLambda;
        }
        //tangent
        if(fri>0) {
          off0=off;
          lmt=dfN*fri;
          for(int d=0; d<2; d++)
            contactTangent(state._pos,m,p,d,off);
          project(state._pos,off0,off,&lmt,true);
          for(int k=off0; k<off; k++)
            fT[k-off0]=_Css[k]._lambda;
        }
        //update force
        p._fA=-(p._nA2B*(GEOMETRY_SCALAR)fN+p._tA2B*fT.template cast<GEOMETRY_SCALAR>());
        p._fB=-p._fA;
      }
    //project drag
    for(auto& joint:_joints) {
      off0=off;
      joint.constraint(*_body,state._pos,nrD,_Css,off);
      project(state._pos,off0,off,NULL,true);
    }
  }
}
void XPBDSimulator::contactNormal(const GradInfo& newPos,const ContactManifold& m,const ContactPoint& p,int& off) {
  //allocate constraint
  if((int)_Css.size()<=off) {
    _Css.push_back(XPBDConstraint());
    _Css[off]._JC.setZero(_body->nrDOF());
  }
  //build position constraint
  XPBDConstraint& C=_Css[off++];
  Mat3X4T GK;
  Vec3T nA2B=p._nA2B.template cast<T>();
  Vec3T ptA=p._ptA.template cast<T>(),ptALast=ptA;
  Vec3T ptB=p._ptB.template cast<T>(),ptBLast=ptB;
  C._JC.setZero();
  if(m._jidA>=0) {
    ptA=ROTI(newPos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(newPos._TM,m._jidA);
    ptALast=ROTI(_pos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(_pos._TM,m._jidA);
    ROT(GK)=nA2B*p._ptAL.template cast<T>().transpose();
    CTR(GK)=nA2B;
    newPos.DTG(m._jidA,*_body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  if(m._jidB>=0) {
    ptB=ROTI(newPos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(newPos._TM,m._jidB);
    ptBLast=ROTI(_pos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(_pos._TM,m._jidB);
    ROT(GK)=-nA2B*p._ptBL.template cast<T>().transpose();
    CTR(GK)=-nA2B;
    newPos.DTG(m._jidB,*_body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  C._C=std::max<T>((ptA-ptB).dot(nA2B),0);
  C._alpha=std::max<T>(m._jidA>=0?_params[m._jidA]._kc:0,m._jidB>=0?_params[m._jidB]._kc:0);
  //build velocity constraint
  T C0=std::max<T>((ptALast-ptBLast).dot(nA2B),0);
  if(C0==0)
    return;
  //allocate constraint
  if((int)_Css.size()<=off) {
    _Css.push_back(XPBDConstraint());
    _Css[off]._JC.setZero(_body->nrDOF());
  }
  //build constraint
  XPBDConstraint& CV=_Css[off++];
  CV._C=((ptA-ptALast)-(ptB-ptBLast)).dot(nA2B);
  CV._JC=_Css[off-2]._JC;
  CV._alpha=_Css[off-2]._alpha;
}
void XPBDSimulator::contactTangent(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,int d,int& off) {
  //allocate constraint
  if((int)_Css.size()<=off) {
    _Css.push_back(XPBDConstraint());
    _Css[off]._JC.setZero(_body->nrDOF());
  }
  //build constraint
  XPBDConstraint& C=_Css[off++];
  Mat3X4T GK;
  Vec3T tA2B=p._tA2B.col(d).template cast<T>();
  Vec3T ptA=p._ptA.template cast<T>(),ptALast=ptA;
  Vec3T ptB=p._ptB.template cast<T>(),ptBLast=ptB;
  C._JC.setZero();
  if(m._jidA>=0) {
    ptA=ROTI(newPos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(newPos._TM,m._jidA);
    ptALast=ROTI(_pos._TM,m._jidA)*p._ptAL.template cast<T>()+CTRI(_pos._TM,m._jidA);
    ROT(GK)=tA2B*p._ptAL.template cast<T>().transpose();
    CTR(GK)=tA2B;
    newPos.DTG(m._jidA,*_body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  if(m._jidB>=0) {
    ptB=ROTI(newPos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(newPos._TM,m._jidB);
    ptBLast=ROTI(_pos._TM,m._jidB)*p._ptBL.template cast<T>()+CTRI(_pos._TM,m._jidB);
    ROT(GK)=-tA2B*p._ptBL.template cast<T>().transpose();
    CTR(GK)=-tA2B;
    newPos.DTG(m._jidB,*_body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  C._C=((ptA-ptALast)-(ptB-ptBLast)).dot(tA2B);
  C._alpha=std::max<T>(m._jidA>=0?_params[m._jidA]._kc:0,m._jidB>=0?_params[m._jidB]._kc:0);
}
template <int N>
void XPBDSimulator::project(GradInfo& newPos,int beg,int end,T* limit,bool updateKinematics) {
  Eigen::Matrix<T,-1,1,0,N,1> RHS;
  Eigen::Matrix<T,-1,-1,0,N,N> JTInvMJ;
  RHS.setZero(end-beg);
  JTInvMJ.setZero(end-beg,end-beg);
  ASSERT_MSG(end-beg<=4,"We can only couple at most 4 constraints!")
  if(end<=beg)
    return;
  for(int k=beg; k<end; k++) {
    XPBDConstraint& C=_Css[k];
    mask(NULL,&(C._JC),NULL);
    C._invMJC=_sol->solve(MatT(C._JC));
    RHS[k-beg]=C._C-C._lambda/C._alpha;
    for(int k2=k; k2<end; k2++) {
      const XPBDConstraint& C2=_Css[k2];
      JTInvMJ(k2-beg,k-beg)=JTInvMJ(k-beg,k2-beg)=C2._JC.dot(C._invMJC);
    }
    JTInvMJ(k-beg,k-beg)+=1/C._alpha;
  }
  //compute step size
  RHS=JTInvMJ.inverse()*RHS;
  //handle limit
  if(limit) {
    for(int k=beg; k<end; k++)
      if(*limit<=0)
        RHS[k-beg]=0;
      else RHS[k-beg]=std::max<T>(std::min<T>(RHS[k-beg],*limit),-*limit);
  }
  //update position/lambda
  for(int k=beg; k<end; k++) {
    XPBDConstraint& C=_Css[k];
    C._dLambda=RHS[k-beg];
    newPos._xM-=C._invMJC*C._dLambda;
    C._lambda+=C._dLambda;
  }
  //update position
  if(updateKinematics)
    newPos.reset(*_body,newPos._xM);
}
}
