#include "XPBDSimulator.h"
#include "PBDMatrixSolver.h"
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>

namespace PHYSICSMOTION {
XPBDSimulator::XPBDSimulator(T dt):PBDSimulator(dt) {
  _maxIt=8;
}
void XPBDSimulator::setArticulatedBody(std::shared_ptr<ArticulatedBody> body) {
  PBDSimulator::setArticulatedBody(body);
  _lambdaLimit.setZero(body->nrDOF());
}
void XPBDSimulator::debugEnergy(T scale) {
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
  Vec x,r,DE,DE2;
  _lambdaLimit=Vec::Random(_body->nrDOF());
  T eRef=energy(newPos,DE,DDE,true);
  T e=PBDSimulator::energy(newPos,DE2,DDE2,true);
  GaussSeidel(newPos,x,r,DE,DDE,0,1e-8f,1e4,2);
  schurComplement(DE,DDE);
  DEBUG_GRADIENT("E",e,e-eRef)
  DEBUG_GRADIENT("DE",DE.norm(),(DE-DE2).norm())
  DEBUG_GRADIENT("DDE",DDE.norm(),(DDE-DDE2).norm())
}
//helper
void XPBDSimulator::update(const GradInfo& newPos,GradInfo& newPos2,const Vec& DE,const MatT& DDE,T alpha) const {
  Vec x,r;
  int nrD=_body->nrDOF();
  GaussSeidel(newPos,x,r,DE,DDE,alpha,1e-4f,100,0);
  newPos2.reset(*_body,newPos._xM+x.segment(0,nrD));
}
void XPBDSimulator::schurComplement(Vec& DE,MatT& DDE) const {
  int nrD=_body->nrDOF(),nrDC=DE.size(),nrC=nrDC-nrD;
  DE=(DE.segment(0,nrD)-DDE.block(0,nrD,nrD,nrC)*DDE.block(nrD,nrD,nrC,nrC).inverse()*DE.segment(nrD,nrC)).eval();
  DDE=(DDE.block(0,0,nrD,nrD)-DDE.block(0,nrD,nrD,nrC)*DDE.block(nrD,nrD,nrC,nrC).inverse()*DDE.block(nrD,0,nrC,nrD)).eval();
}
void XPBDSimulator::updateLambda(const GradInfo& newPos,Vec& lambdaLimit,VecCM dx) {
  int nrD=_body->nrDOF(),nrDC=nrD;
  jointLimitConstraint(newPos,lambdaLimit,mapM((MatT*)NULL),mapV((Vec*)NULL),dx,nrDC);
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      normalConstraint(newPos,m,p,mapM((MatT*)NULL),mapV((Vec*)NULL),dx,nrDC);
      tangentConstraint(newPos,m,p,mapM((MatT*)NULL),mapV((Vec*)NULL),dx,nrDC,false);
    }
  for(auto& d:_drags)
    dragConstraint(newPos,d,mapM((MatT*)NULL),mapV((Vec*)NULL),dx,nrDC);
}
void XPBDSimulator::GaussSeidel(const GradInfo& newPos,Vec& x,Vec& r,const Vec& DE,MatT DDE,T alpha,T tol,int maxIter,int debug) const {
  T dbi;
  DEFINE_NUMERIC_DELTA_T(T)
  int nrD=_body->nrDOF(),nrC=DE.size()-nrD;
  x.setZero(nrD+nrC);
  r.setZero(nrD+nrC);
  //Use one pass of guass-seidel to solve: DDE*x+DE=0.
  //We use temporary variable: r=DDE*x+DX,
  //and the goal of gauss-seidel is to modify x by dx, so we have: DDE*dx+r=0.
  DDE.diagonal().segment(0,nrD).array()+=alpha;
  if(std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol))
    std::dynamic_pointer_cast<PBDMatrixSolverABA>(_sol)->compute(Vec(newPos._xM),_MRR,_MRt,_MtR,_Mtt,Vec(_diag+Vec::Constant(nrD,alpha)));
  else _sol->compute(MatT(DDE.block(0,0,nrD,nrD)));
  Eigen::Block<MatT> C=DDE.block(nrD,0,nrC,nrD);
  Eigen::Block<MatT> H=DDE.block(nrD,nrD,nrC,nrC);
  x.segment(0,nrD)=-_sol->solve(MatT(DE.segment(0,nrD)));
  r.segment(nrD,nrC)=C*x.segment(0,nrD)+DE.segment(nrD,nrC);
  MatT invMCT=_sol->solve(MatT(C.transpose())),CInvMCT=C*invMCT;
  if(debug>=2) {
    DEBUG_GRADIENT("||DDE*x+DE-r||-before",r.norm(),(DDE*x+DE-r).norm())
  }
  //We first write out the KKT system:
  //    dxa=da
  //    dxb=db
  //    DDEa=(M -C^T)(da)
  //    DDEb=(C   H )(db),
  //so we have:
  //    0=Mda-C^Tdb+ra
  //    0=Cda+ H db+rb.
  //Apply Gauss elimination, we have:
  //    da=M^{-1}(C^Tdb-ra)
  //    0=CM^{-1}(C^Tdb-ra)+Hdb+rb.
  //Rearrange and we have:
  //    0=(CM^{-1}C^T+H)db+(rb-CM^{-1}ra).
  //Now if we modify the ith entry of db, then:
  //    dbi=(CiM^{-1}*ra-rbi)/(CiM^{-1}Ci^T+Hii)
  //       =-rbi/(CiM^{-1}Ci^T+Hii),
  //where we use the fact that ra=0, to make sure this is always true, we need to update x by:
  //    da=M^{-1}Ci^Tdbi,
  //and update rb as:
  //    rb+=CM^{-1}Ci^Tdbi+ei*Hii*dbi.
  for(int iter=0; iter<maxIter; iter++) {
    for(int off=nrD,i=0; i<nrC; i++,off++) {
      dbi=-r[off]/(CInvMCT(i,i)+H(i,i));
      //update x
      x[off]+=dbi;
      //x.segment(0,nrD)+=invMCT.col(i)*dbi;    //we can delay this update to the end
      //update r
      r.segment(nrD,nrC)+=CInvMCT.col(i)*dbi;
      r[off]+=H(i,i)*dbi;
    }
    if(debug>=1)
      std::cout << "Gauss-Seidel iter=" << iter << " residual=" << r.segment(nrD,nrC).cwiseAbs().maxCoeff() << std::endl;
    if(r.cwiseAbs().maxCoeff()<tol)
      break;
  }
  x.segment(0,nrD)+=invMCT*x.segment(nrD,nrC);  //this is the delayed update, more efficient
  if(debug>=2) {
    DEBUG_GRADIENT("||DDE*x+DE-r||-after",r.norm(),(DDE*x+DE-r).norm())
  }
}
XPBDSimulator::T XPBDSimulator::energy(const GradInfo& newPos,Vec& DE,MatT& DDE,bool updateTangentBound) {
  //count constraint
  int nrD=_body->nrDOF(),nrDC=nrD;
  jointLimitConstraint(newPos,_lambdaLimit,mapM((MatT*)NULL),mapV((Vec*)NULL),mapCV((const Vec*)NULL),nrDC);
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      normalConstraint(newPos,m,p,mapM((MatT*)NULL),mapV((Vec*)NULL),mapCV((const Vec*)NULL),nrDC);
      tangentConstraint(newPos,m,p,mapM((MatT*)NULL),mapV((Vec*)NULL),mapCV((const Vec*)NULL),nrDC,updateTangentBound);
    }
  for(auto& d:_drags)
    dragConstraint(newPos,d,mapM((MatT*)NULL),mapV((Vec*)NULL),mapCV((const Vec*)NULL),nrDC);
  //resize
  T E=0;
  DE.setZero(nrDC);
  DDE.setZero(nrDC,nrDC);
  nrDC=nrD;
  //fill-in data
  E+=dynamics(newPos,mapM(DDE),mapV(DE));
  E+=jointLimitConstraint(newPos,_lambdaLimit,mapM(DDE),mapV(DE),mapCV((const Vec*)NULL),nrDC);
  for(auto& m:_manifolds)
    for(auto& p:m._points) {
      E+=normalConstraint(newPos,m,p,mapM(DDE),mapV(DE),mapCV((const Vec*)NULL),nrDC);
      E+=tangentConstraint(newPos,m,p,mapM(DDE),mapV(DE),mapCV((const Vec*)NULL),nrDC,updateTangentBound);
    }
  for(auto& d:_drags)
    E+=dragConstraint(newPos,d,mapM(DDE),mapV(DE),mapCV((const Vec*)NULL),nrDC);
  return E;
}
XPBDSimulator::T XPBDSimulator::dynamics(const GradInfo& newPos,MatTM M,VecM g) {
  Vec tmp;
  Vec3T P;
  Mat3X4T A;
  Mat3T PPT,H;
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  T coef=1.0/(_dt*_dt),E=0;
  int nrJ=_body->nrJ();
  int nrD=_body->nrDOF();
  G.setZero(3,4*nrJ);
  _MRR.setZero(3,3*nrJ);
  _MRt.setZero(3,3*nrJ);
  _MtR.setZero(3,3*nrJ);
  _Mtt.setZero(3,3*nrJ);
  _diag.setZero(nrD);
  //dynamic
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    _MRR.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(newPos._TM,k)*J._MCCT.template cast<T>()*ROTI(newPos._TM,k).transpose())*coef;
    _MRt.template block<3,3>(0,k*3)+=cross<T>(ROTI(newPos._TM,k)*J._MC.template cast<T>())*coef;
    _MtR.template block<3,3>(0,k*3)-=cross<T>(ROTI(newPos._TM,k)*J._MC.template cast<T>())*coef;
    _Mtt.template block<3,3>(0,k*3)+=Mat3T::Identity()*J._M*coef;
  }
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
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
      g.segment(J._offDOF,nrD)+=tmp*_params[k]._kp;
      _diag.segment(J._offDOF,nrD).array()+=_params[k]._kp;
    }
    //D controller
    if(_params[k]._kd>0) {
      tmp=(newPos._xM-_pos._xM).segment(J._offDOF,nrD)/_dt-_params[k]._tarD(_t,nrD);
      E+=tmp.squaredNorm()*_params[k]._kd/2;
      g.segment(J._offDOF,nrD)+=tmp*_params[k]._kd/_dt;
      _diag.segment(J._offDOF,nrD).array()+=_params[k]._kd/_dt/_dt;
    }
  }
  //gradient
  newPos.DTG(*_body,mapM(GB=G),g);
  //hessian
  newPos.toolA(*_body,newPos,mapM(MRR=_MRR),mapM(MRt=_MRt),mapM(MtR=_MtR),mapM(Mtt=_Mtt),[&](int r,int c,T val) {
    M(r,c)+=val;
  });
  M.diagonal().segment(0,_diag.size())+=_diag;
  return E;
}
XPBDSimulator::T XPBDSimulator::jointLimitConstraint(const GradInfo& newPos,Vec& lambdaLimit,MatTM H,VecM g,VecCM dx,int& off) const {
  T E=0,val;
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    for(int c=0; c<J._limits.cols(); c++)
      if(isfinite(J._limits(2,c)) && J._limits(2,c)>0) {
        if(isfinite(J._limits(0,c)) && newPos._xM[J._offDOF+c]<J._limits(0,c)) {
          //lower limited
          val=newPos._xM[J._offDOF+c]-J._limits(0,c);
          if(H.data() && g.data()) {
            H(off,J._offDOF+c)+=1;
            H(J._offDOF+c,off)-=1;
            g[J._offDOF+c]-=lambdaLimit[J._offDOF+c];
            H(off,off)+=1/J._limits(2,c);
            g[off]+=val+lambdaLimit[J._offDOF+c]/J._limits(2,c);
            E+=val*val*J._limits(2,c)/2;
          } else if(dx.data())
            lambdaLimit[J._offDOF+c]+=dx[off];
          off++;
        } else if(isfinite(J._limits(1,c)) && newPos._xM[J._offDOF+c]>J._limits(1,c)) {
          //upper limited
          val=newPos._xM[J._offDOF+c]-J._limits(1,c);
          if(H.data() && g.data()) {
            H(off,J._offDOF+c)+=1;
            H(J._offDOF+c,off)-=1;
            g[J._offDOF+c]-=lambdaLimit[J._offDOF+c];
            H(off,off)+=1/J._limits(2,c);
            g[off]+=val+lambdaLimit[J._offDOF+c]/J._limits(2,c);
            E+=val*val*J._limits(2,c)/2;
          } else if(dx.data())
            lambdaLimit[J._offDOF+c]+=dx[off];
          off++;
        } else {
          //clear lambda
          if(dx.data())
            lambdaLimit[J._offDOF+c]=0;
        }
      }
  }
  return E;
}
XPBDSimulator::T XPBDSimulator::normalConstraint(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,MatTM H,VecM g,VecCM dx,int& off) const {
  T E=0;
  Mat3X4T G;
  p._fA.setZero();
  p._fB.setZero();
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
    if(m._jidA>=0)
      p._fA-=(kc*depth*nA2B).template cast<GEOMETRY_SCALAR>();
    if(m._jidB>=0)
      p._fB+=(kc*depth*nA2B).template cast<GEOMETRY_SCALAR>();
    if(H.data() && g.data()) {
      if(m._jidA>=0) {
        G=nA2B*Vec4T((T)p._ptAL[0],(T)p._ptAL[1],(T)p._ptAL[2],1).transpose();
        newPos.DTG(m._jidA,*_body,G,[&](int c,T val) {
          H(off,c)+=val;
          H(c,off)-=val;
          g[c]-=val*(T)p._lambda[0];
        });
      }
      if(m._jidB>=0) {
        G=-nA2B*Vec4T((T)p._ptBL[0],(T)p._ptBL[1],(T)p._ptBL[2],1).transpose();
        newPos.DTG(m._jidB,*_body,G,[&](int c,T val) {
          H(off,c)+=val;
          H(c,off)-=val;
          g[c]-=val*(T)p._lambda[0];
        });
      }
      H(off,off)+=1/kc;
      g[off]+=depth+(T)p._lambda[0]/kc;
      E+=kc/2*depth*depth;
    } else if(dx.data())
      p._lambda[0]+=(GEOMETRY_SCALAR)dx[off];
    off++;
  }
  return E;
}
XPBDSimulator::T XPBDSimulator::tangentConstraint(const GradInfo& newPos,const ContactManifold& m,ContactPoint& p,MatTM H,VecM g,VecCM dx,int& off,bool updateTangentBound) const {
  T E=0;
  Mat3X4T G;
  Mat2T H2,invH2;
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
    Vec2T dPT=tA2B.transpose()*((ptA-ptALast)-(ptB-ptBLast)),f1;
    T dPTLen=dPT.norm(),dPTLen2=dPTLen*dPTLen,dPTLen3=dPTLen2*dPTLen;
    T epsVDt=_epsV*_dt,invEpsVDt=1/epsVDt,invEpsVDt2=invEpsVDt*invEpsVDt,f0,tb=(T)p._tangentBound;
    if(dPTLen<epsVDt) {
      f0=-dPTLen3*invEpsVDt2/3+dPTLen2*invEpsVDt+epsVDt/3;
      f1=(-dPTLen*invEpsVDt2+2*invEpsVDt)*dPT;
      H2=Mat2T::Identity()*(-dPTLen*invEpsVDt2+2*invEpsVDt)*tb;
      invH2=Mat2T::Identity()/(-dPTLen*invEpsVDt2+2*invEpsVDt)/tb;
    } else {
      f0=dPTLen;
      f1=dPT/f0;
      H2=Mat2T::Identity()/f0*tb;
      invH2=Mat2T::Identity()*f0/tb;
    }
    if(m._jidA>=0)
      p._fA-=(tA2B*f1).template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
    if(m._jidB>=0)
      p._fB+=(tA2B*f1).template cast<GEOMETRY_SCALAR>()*(GEOMETRY_SCALAR)tb;
    if(H.data() && g.data()) {
      for(int d=0; d<2; d++) {
        if(m._jidA>=0) {
          G=tA2B.col(d)*Vec4T((T)p._ptAL[0],(T)p._ptAL[1],(T)p._ptAL[2],1).transpose();
          newPos.DTG(m._jidA,*_body,G,[&](int c,T val) {
            H(off+d,c)+=val;
            H(c,off+d)-=val;
            g[c]-=val*(T)p._lambda[d+1];
          });
        }
        if(m._jidB>=0) {
          G=-tA2B.col(d)*Vec4T((T)p._ptBL[0],(T)p._ptBL[1],(T)p._ptBL[2],1).transpose();
          newPos.DTG(m._jidB,*_body,G,[&](int c,T val) {
            H(off+d,c)+=val;
            H(c,off+d)-=val;
            g[c]-=val*(T)p._lambda[d+1];
          });
        }
      }
      H.template block<2,2>(off,off)+=invH2;
      g.template segment<2>(off)+=invH2*(p._lambda.template segment<2>(1).template cast<T>()+f1*tb);
      E+=f0*tb;
    } else if(dx.data())
      p._lambda.template segment<2>(1)+=dx.template segment<2>(off).template cast<GEOMETRY_SCALAR>();
    off+=2;
  }
  return E;
}
XPBDSimulator::T XPBDSimulator::dragConstraint(const GradInfo& newPos,DragEnergy& drag,MatTM H,VecM g,VecCM dx,int& off) const {
  T E=0,val;
  Mat3X4T G;
  Vec3T pt=ROTI(newPos._TM,drag._jid)*drag._ptL+CTRI(newPos._TM,drag._jid);
  for(int d=0; d<3; d++) {
    val=pt[d]-drag._pt[d];
    if(H.data() && g.data()) {
      G.setZero();
      G.row(d)=Vec4T(drag._ptL[0],drag._ptL[1],drag._ptL[2],1).transpose();
      newPos.DTG(drag._jid,*_body,G,[&](int c,T val) {
        H(off,c)+=val;
        H(c,off)-=val;
        g[c]-=val*drag._lambda[d];
      });
      H(off,off)+=1/drag._k;
      g[off]+=val+drag._lambda[d]/drag._k;
      E+=val*val*drag._k/2;
    } else if(dx.data())
      drag._lambda[d]+=dx[off];
    off++;
  }
  return E;
}
}
