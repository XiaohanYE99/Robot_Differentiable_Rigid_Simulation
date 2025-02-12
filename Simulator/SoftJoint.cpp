#include "SoftJoint.h"
#include <Utils/RotationUtils.h>

namespace PHYSICSMOTION {
SoftJoint::SoftJoint() {
  _linearLimit.setConstant(std::numeric_limits<double>::infinity());
  _angularLimit.setConstant(std::numeric_limits<double>::infinity());
  _transA.setIdentity();
  _transB.setIdentity();
}
SoftJoint::Vec3T SoftJoint::posA(const GradInfo& pos) const {
  if(_jidA==-1)
    return CTR(_transA);
  else return ROTI(pos._TM,_jidA)*CTR(_transA)+CTRI(pos._TM,_jidA);
}
SoftJoint::Vec3T SoftJoint::posB(const GradInfo& pos) const {
  if(_jidB==-1)
    return CTR(_transB);
  else return ROTI(pos._TM,_jidB)*CTR(_transB)+CTRI(pos._TM,_jidB);
}
SoftJoint::Vec3T SoftJoint::rotA(const GradInfo& pos,int d) const {
  if(_jidA==-1)
    return rotAL(d);
  else return ROTI(pos._TM,_jidA)*rotAL(d);
}
SoftJoint::Vec3T SoftJoint::rotB(const GradInfo& pos,int d) const {
  if(_jidB==-1)
    return rotBL(d);
  else return ROTI(pos._TM,_jidB)*rotBL(d);
}
SoftJoint::Vec3T SoftJoint::rotAL(int d) const {
  return ROT(_transA)*Vec3T::Unit(d);
}
SoftJoint::Vec3T SoftJoint::rotBL(int d) const {
  return ROT(_transB)*Vec3T::Unit(d);
}
void SoftJoint::setRandom(const ArticulatedBody& body) {
  //randomize joint id
  while(true) {
    _jidA=rand()%(body.nrJ()+1)-1;
    _jidB=rand()%(body.nrJ()+1)-1;
    if(_jidA!=_jidB)
      break;
  }
  //randomize transformation
  ROT(_transA)=expWGradV<T,Vec3T>(Vec3T::Random());
  CTR(_transA).setRandom();
  ROT(_transB)=expWGradV<T,Vec3T>(Vec3T::Random());
  CTR(_transB).setRandom();
  //randomize linear limit
  _linearLimit.setConstant(std::numeric_limits<double>::infinity());
  _linearLimit.row(1).setConstant(rand()/(double)RAND_MAX);
  _linearLimit.row(2).setConstant(rand()/(double)RAND_MAX);
  //randomize angular limit
  _angularLimit.setConstant(std::numeric_limits<double>::infinity());
  for(int c=0; c<3; c++) {
    _angularLimit(1,c)=rand()/(double)RAND_MAX;
    _angularLimit(2,c)=rand()/(double)RAND_MAX;
  }
}
//energy model
SoftJoint::T SoftJoint::energyLinear(const ArticulatedBody& body,const GradInfo& newPos,int d,Vec* DE,MatT& DDE,Mat3XT& G,
                                     Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm) const {
  Vec3T D=Vec3T::Unit(d);
  Mat3T ptARC,ptBRC,H,HRR,HRt,HtR;
  T dpt=(posA(newPos)-posB(newPos)).dot(D),val;
  if(isfinite(_linearLimit(0,d)) && dpt<_linearLimit(0,d))
    dpt-=_linearLimit(0,d);
  else if(isfinite(_linearLimit(1,d)) && dpt>_linearLimit(1,d))
    dpt-=_linearLimit(1,d);
  else dpt=0;
  if(!isfinite(_linearLimit(2,d)) || dpt==0)
    return 0;
  if(DE) {
    H=D*D.transpose()*_linearLimit(2,d);
    if(_jidA>=0) {
      ROTI(G,_jidA)+=D*CTR(_transA).transpose()*dpt*_linearLimit(2,d);
      CTRI(G,_jidA)+=D*dpt*_linearLimit(2,d);
      ptARC=cross<T>(ROTI(newPos._TM,_jidA)*CTR(_transA));
      MRR.template block<3,3>(0,_jidA*3)+=ptARC*H*ptARC.transpose();
      MRt.template block<3,3>(0,_jidA*3)+=HRt=ptARC*H;
      MtR.template block<3,3>(0,_jidA*3)+=HRt.transpose();
      Mtt.template block<3,3>(0,_jidA*3)+=H;
    }
    if(_jidB>=0) {
      ROTI(G,_jidB)-=D*CTR(_transB).transpose()*dpt*_linearLimit(2,d);
      CTRI(G,_jidB)-=D*dpt*_linearLimit(2,d);
      ptBRC=cross<T>(ROTI(newPos._TM,_jidB)*CTR(_transB));
      MRR.template block<3,3>(0,_jidB*3)+=ptBRC*H*ptBRC.transpose();
      MtR.template block<3,3>(0,_jidB*3)+=HtR=H*ptBRC.transpose();
      MRt.template block<3,3>(0,_jidB*3)+=HtR.transpose();
      Mtt.template block<3,3>(0,_jidB*3)+=H;
    }
    if(_jidA>=0 && _jidB>=0 && crossTerm) {
      HRR=ptARC*H*ptBRC.transpose();
      newPos.JRCSparse(body,_jidA,[&](int r,const Vec3T& JRA) {
        newPos.JRCSparse(body,_jidB,[&](int c,const Vec3T& JRB) {
          DDE(r,c)-=val=JRA.dot(HRR*JRB);
          DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          DDE(r,c)-=val=JRA.dot(HRt*JtB);
          DDE(c,r)-=val;
        });
      },[&](int r,const Vec3T& JtA) {
        newPos.JRCSparse(body,_jidB,[&](int c,const Vec3T& JRB) {
          DDE(r,c)-=val=JtA.dot(HtR*JRB);
          DDE(c,r)-=val;
        },[&](int c,const Vec3T& JtB) {
          DDE(r,c)-=val=JtA.dot(H*JtB);
          DDE(c,r)-=val;
        });
      });
    }
  }
  return dpt*dpt*_linearLimit(2,d)/2;
}
SoftJoint::T SoftJoint::energyBall(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,Mat3XT& G,
                                   Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm,T tol) const {
  //build stiffness and coefficient
  T k=0,lmt=0,szSqr,sz,val,E=0;
  Vec3T coeff,D;
  coeff.setZero(3);
  for(int c=0; c<3; c++) {
    //ball
    if(!isfinite(_linearLimit(0,c)) && isfinite(_linearLimit(1,c)) && isfinite(_linearLimit(2,c)) && _linearLimit(2,c)>0 && _linearLimit(1,c)>=0) {
      ASSERT_MSG(k==0 || k==_linearLimit(2,c),"We only support ball-shaped translational joint limit!");
      ASSERT_MSG(lmt==0 || lmt==_linearLimit(1,c),"We only support ball-shaped translational joint limit!");
      k=_linearLimit(2,c);
      lmt=_linearLimit(1,c);
      coeff[c]=1;
    }
    //lock
    if(isfinite(_linearLimit(0,c)) && isfinite(_linearLimit(1,c)) && isfinite(_linearLimit(2,c)) && _linearLimit(2,c)>0 && _linearLimit(1,c)==0) {
      ASSERT_MSG(_linearLimit(0,c)==0,"Asymmetric translational joint limit not supported!");
      E+=energyLinear(body,newPos,c,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm);
    }
  }
  Mat3T ptARC,ptBRC,H,HRR,HRt,HtR;
  Vec3T dpt=posA(newPos)-posB(newPos);
  szSqr=(dpt.array()*coeff.array()).matrix().squaredNorm()+tol;
  sz=sqrt(szSqr);
  if(sz>lmt) {
    if(DE) {
      D=(dpt.array()*coeff.array()*coeff.array()).matrix()/sz;
      H=k*D*D.transpose();
      H+=(sz-lmt)*k*(coeff.array()*coeff.array()/sz).matrix().asDiagonal();
      H-=(sz-lmt)*k*(D*D.transpose())/sz;
      if(_jidA>=0) {
        ROTI(G,_jidA)+=D*CTR(_transA).transpose()*(sz-lmt)*k;
        CTRI(G,_jidA)+=D*(sz-lmt)*k;
        ptARC=cross<T>(ROTI(newPos._TM,_jidA)*CTR(_transA));
        MRR.template block<3,3>(0,_jidA*3)+=ptARC*H*ptARC.transpose();
        MRt.template block<3,3>(0,_jidA*3)+=HRt=ptARC*H;
        MtR.template block<3,3>(0,_jidA*3)+=HRt.transpose();
        Mtt.template block<3,3>(0,_jidA*3)+=H;
      }
      if(_jidB>=0) {
        ROTI(G,_jidB)-=D*CTR(_transB).transpose()*(sz-lmt)*k;
        CTRI(G,_jidB)-=D*(sz-lmt)*k;
        ptBRC=cross<T>(ROTI(newPos._TM,_jidB)*CTR(_transB));
        MRR.template block<3,3>(0,_jidB*3)+=ptBRC*H*ptBRC.transpose();
        MtR.template block<3,3>(0,_jidB*3)+=HtR=H*ptBRC.transpose();
        MRt.template block<3,3>(0,_jidB*3)+=HtR.transpose();
        Mtt.template block<3,3>(0,_jidB*3)+=H;
      }
      if(_jidA>=0 && _jidB>=0 && crossTerm) {
        HRR=ptARC*H*ptBRC.transpose();
        newPos.JRCSparse(body,_jidA,[&](int r,const Vec3T& JRA) {
          newPos.JRCSparse(body,_jidB,[&](int c,const Vec3T& JRB) {
            DDE(r,c)-=val=JRA.dot(HRR*JRB);
            DDE(c,r)-=val;
          },[&](int c,const Vec3T& JtB) {
            DDE(r,c)-=val=JRA.dot(HRt*JtB);
            DDE(c,r)-=val;
          });
        },[&](int r,const Vec3T& JtA) {
          newPos.JRCSparse(body,_jidB,[&](int c,const Vec3T& JRB) {
            DDE(r,c)-=val=JtA.dot(HtR*JRB);
            DDE(c,r)-=val;
          },[&](int c,const Vec3T& JtB) {
            DDE(r,c)-=val=JtA.dot(H*JtB);
            DDE(c,r)-=val;
          });
        });
      }
    }
    E+=(sz-lmt)*(sz-lmt)*k/2;
  }
  return E;
}
SoftJoint::T SoftJoint::energyAngular(const ArticulatedBody& body,const GradInfo& newPos,int d,Vec* DE,MatT& DDE,Mat3XT& G,
                                      Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm) const {
  ASSERT_MSG(!isfinite(_angularLimit(0,d)),"We do not support angular lower bound!")
  if(!isfinite(_angularLimit(1,d)) || _angularLimit(1,d)<0)
    return 0;
  if(!isfinite(_angularLimit(2,d)) || _angularLimit(2,d)<=0)
    return 0;
  Mat3T HRR;
  Vec3T dA=rotA(newPos,d),dB=rotB(newPos,d),dACdB;
  T sz=dA.dot(dB),lmt=cos(_angularLimit(1,d)),val;
  if(sz<lmt) {
    if(DE) {
      dACdB=dA.cross(dB);
      HRR=dACdB*dACdB.transpose();
      if(_jidA>=0) {
        ROTI(G,_jidA)+=dB*rotAL(d).transpose()*(sz-lmt)*_angularLimit(2,d);
        MRR.template block<3,3>(0,_jidA*3)+=HRR*_angularLimit(2,d);
      }
      if(_jidB>=0) {
        ROTI(G,_jidB)+=dA*rotBL(d).transpose()*(sz-lmt)*_angularLimit(2,d);
        MRR.template block<3,3>(0,_jidB*3)+=HRR*_angularLimit(2,d);
      }
      if(_jidA>=0 && _jidB>=0 && crossTerm) {
        HRR+=cross(dA)*cross(dB)*(sz-lmt);
        HRR*=_angularLimit(2,d);
        newPos.JRSparse(body,_jidA,[&](int r,const Vec3T& JRA) {
          newPos.JRSparse(body,_jidB,[&](int c,const Vec3T& JRB) {
            DDE(r,c)-=val=JRA.dot(HRR*JRB);
            DDE(c,r)-=val;
          });
        });
      }
    }
    return (sz-lmt)*(sz-lmt)*_angularLimit(2,d)/2;
  } else return 0;
}
SoftJoint::T SoftJoint::energy(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,Mat3XT& G,
                               Mat3XT& MRR,Mat3XT& MRt,Mat3XT& MtR,Mat3XT& Mtt,bool crossTerm,T tol) const {
  return energyBall(body,newPos,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm,tol)+
         energyAngular(body,newPos,0,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm)+
         energyAngular(body,newPos,1,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm)+
         energyAngular(body,newPos,2,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm);
}
SoftJoint::T SoftJoint::energy(const ArticulatedBody& body,const GradInfo& newPos,Vec* DE,MatT& DDE,bool JTJ,bool crossTerm,T tol) const {
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  Mat3XT MRRB,MRtB,MtRB,MttB;
  G.setZero(3,body.nrJ()*4);
  MRR.setZero(3,body.nrJ()*3);
  MRt.setZero(3,body.nrJ()*3);
  MtR.setZero(3,body.nrJ()*3);
  Mtt.setZero(3,body.nrJ()*3);
  T ret=energy(body,newPos,DE,DDE,G,MRR,MRt,MtR,Mtt,crossTerm,tol);
  if(DE) {
    //gradient
    newPos.DTG(body,mapM(GB=G),mapV(DE));
    //hessian
    if(JTJ) {
      newPos.toolA(body,newPos,mapM(MRRB=MRR),mapM(MRtB=MRt),mapM(MtRB=MtR),mapM(MttB=Mtt),[&](int r,int c,T val) {
        DDE(r,c)+=val;
      });
    } else {
      newPos.toolAB(body,mapM(MRRB=MRR),mapM(MRtB=MRt),mapM(MtRB=MtR),mapM(MttB=Mtt),mapM(GB=G),[&](int r,int c,T val) {
        DDE(r,c)+=val;
      });
    }
  }
  return ret;
}
//constraint model
void SoftJoint::constraintLinear(const ArticulatedBody& body,const GradInfo& newPos,int d,int nrDOF,
                                 std::vector<Simulator::XPBDConstraint>& Css,int& off) const {
  Mat3X4T GK;
  Vec3T D=Vec3T::Unit(d);
  //allocate constraint
  if((int)Css.size()<=off) {
    Css.push_back(Simulator::XPBDConstraint());
    Css[off]._JC.setZero(nrDOF);
  }
  //build constraint
  Simulator::XPBDConstraint& C=Css[off++];
  C._JC.setZero();
  C._C=0;
  T dpt=(posA(newPos)-posB(newPos)).dot(D);
  if(isfinite(_linearLimit(0,d)) && dpt<_linearLimit(0,d))
    C._C=dpt-_linearLimit(0,d);
  else if(isfinite(_linearLimit(1,d)) && dpt>_linearLimit(1,d))
    C._C=dpt-_linearLimit(1,d);
  if(!isfinite(_linearLimit(2,d)))
    return;
  if(_jidA>=0) {
    ROT(GK)=D*CTR(_transA).transpose();
    CTR(GK)=D;
    newPos.DTG(_jidA,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  if(_jidB>=0) {
    ROT(GK)=-D*CTR(_transB).transpose();
    CTR(GK)=-D;
    newPos.DTG(_jidB,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  C._alpha=_linearLimit(2,d);
}
void SoftJoint::constraintBall(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,
                               std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol) const {
  //build stiffness and coefficient
  Mat3X4T GK;
  T k=0,lmt=0,szSqr,sz;
  Vec3T coeff,D,dpt;
  coeff.setZero(3);
  for(int c=0; c<3; c++) {
    //ball
    if(!isfinite(_linearLimit(0,c)) && isfinite(_linearLimit(1,c)) && isfinite(_linearLimit(2,c)) && _linearLimit(2,c)>0 && _linearLimit(1,c)>=0) {
      ASSERT_MSG(k==0 || k==_linearLimit(2,c),"We only support ball-shaped translational joint limit!");
      ASSERT_MSG(lmt==0 || lmt==_linearLimit(1,c),"We only support ball-shaped translational joint limit!");
      k=_linearLimit(2,c);
      lmt=_linearLimit(1,c);
      coeff[c]=1;
    }
    //lock
    if(isfinite(_linearLimit(0,c)) && isfinite(_linearLimit(1,c)) && isfinite(_linearLimit(2,c)) && _linearLimit(2,c)>0 && _linearLimit(1,c)==0) {
      ASSERT_MSG(_linearLimit(0,c)==0,"Asymmetric translational joint limit not supported!");
      constraintLinear(body,newPos,c,nrDOF,Css,off);
    }
  }
  if(coeff.isZero())
    return;
  //allocate constraint
  if((int)Css.size()<=off) {
    Css.push_back(Simulator::XPBDConstraint());
    Css[off]._JC.setZero(nrDOF);
  }
  //build constraint
  Simulator::XPBDConstraint& C=Css[off++];
  C._JC.setZero();
  C._C=0;
  dpt=posA(newPos)-posB(newPos);
  szSqr=(dpt.array()*coeff.array()).matrix().squaredNorm()+tol;
  sz=sqrt(szSqr);
  if(sz>lmt)
    C._C=sz-lmt;
  D=(dpt.array()*coeff.array()*coeff.array()).matrix()/sz;
  if(_jidA>=0) {
    ROT(GK)=D*CTR(_transA).transpose();
    CTR(GK)=D;
    newPos.DTG(_jidA,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  if(_jidB>=0) {
    ROT(GK)=-D*CTR(_transB).transpose();
    CTR(GK)=-D;
    newPos.DTG(_jidB,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  C._alpha=k;
}
void SoftJoint::constraintAngular(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,int d,
                                  std::vector<Simulator::XPBDConstraint>& Css,int& off) const {
  ASSERT_MSG(!isfinite(_angularLimit(0,d)),"We do not support angular lower bound!")
  if(!isfinite(_angularLimit(1,d)) || _angularLimit(1,d)<0)
    return;
  if(!isfinite(_angularLimit(2,d)) || _angularLimit(2,d)<=0)
    return;
  Mat3X4T GK=Mat3X4T::Zero();
  Vec3T dA=rotA(newPos,d);
  Vec3T dB=rotB(newPos,d);
  T sz=dA.dot(dB),lmt=cos(_angularLimit(1,d));
  //allocate constraint
  if((int)Css.size()<=off) {
    Css.push_back(Simulator::XPBDConstraint());
    Css[off]._JC.setZero(nrDOF);
  }
  //build constraint
  Simulator::XPBDConstraint& C=Css[off++];
  C._JC.setZero();
  C._C=0;
  if(sz<lmt)
    C._C=sz-lmt;
  if(_jidA>=0) {
    ROT(GK)=dB*rotAL(d).transpose();
    newPos.DTG(_jidA,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  if(_jidB>=0) {
    ROT(GK)=dA*rotBL(d).transpose();
    newPos.DTG(_jidB,body,GK,[&](int r,T val) {
      C._JC[r]+=val;
    });
  }
  C._alpha=_angularLimit(2,d);
}
void SoftJoint::constraint(const ArticulatedBody& body,const GradInfo& newPos,int nrDOF,
                           std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol) const {
  constraintBall(body,newPos,nrDOF,Css,off,tol);
  constraintAngular(body,newPos,nrDOF,0,Css,off);
  constraintAngular(body,newPos,nrDOF,1,Css,off);
  constraintAngular(body,newPos,nrDOF,2,Css,off);
}
//debug
void SoftJoint::debug(const ArticulatedBody& body) {
  //linear
  for(int d=0; d<4; d++)
    for(int pass=0; pass<3; pass++) {
      setRandom(body);
      if(pass==0) {
        if(_jidA!=-1 || _jidB==-1) {
          pass--;
          continue;
        }
      }
      if(pass==1) {
        if(_jidA==-1 || _jidB!=-1) {
          pass--;
          continue;
        }
      }
      if(pass==2) {
        if(_jidA==-1 || _jidB==-1) {
          pass--;
          continue;
        }
      }
      _angularLimit.setConstant(std::numeric_limits<double>::infinity());
      if(d<3)
        _linearLimit(0,d)=_linearLimit(1,d)=0;
      debugInner("BallLock"+std::to_string(d),body);
    }
  //angular
  for(int pass=0; pass<3; pass++) {
    setRandom(body);
    if(pass==0) {
      if(_jidA!=-1 || _jidB==-1) {
        pass--;
        continue;
      }
    }
    if(pass==1) {
      if(_jidA==-1 || _jidB!=-1) {
        pass--;
        continue;
      }
    }
    if(pass==2) {
      if(_jidA==-1 || _jidB==-1) {
        pass--;
        continue;
      }
    }
    _linearLimit.setConstant(std::numeric_limits<double>::infinity());
    debugInner("Angular",body);
  }
}
void SoftJoint::debugInner(const std::string& name,const ArticulatedBody& body) {
  DEFINE_NUMERIC_DELTA_T(T)
  while(true) {
    Vec DE0=Vec::Random(body.nrDOF()),DE=DE0,DE2=DE0;
    MatT DDE0=MatT::Random(body.nrDOF(),body.nrDOF()),DDE=DDE0,DDE2=DDE0;
    Vec x=Vec::Random(body.nrDOF());
    Vec dx=Vec::Random(body.nrDOF());
    GradInfo newPos(body,x);
    GradInfo newPos2(body,x+dx*DELTA);
    T E=energy(body,newPos,&DE,DDE,false,true);
    if(E==0)
      continue;
    T E2=energy(body,newPos2,&DE2,DDE2,false,true);
    DE-=DE0;
    DDE-=DDE0;
    DE2-=DE0;
    DEBUG_GRADIENT(name+"-G",DE.dot(dx),DE.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT(name+"-H",(DDE*dx).norm(),(DDE*dx-(DE2-DE)/DELTA).norm())
    //debug constraint
    int off=0;
    Vec DERef=Vec::Zero(x.size());
    std::vector<Simulator::XPBDConstraint> Css;
    constraint(body,newPos,body.nrDOF(),Css,off);
    for(int k=0; k<off; k++)
      DERef+=Css[k]._JC*Css[k]._C*Css[k]._alpha;
    DEBUG_GRADIENT(name+"-CG",DE.norm(),(DE-DERef).norm())
    //ensure this can be called
    energy(body,newPos,NULL,DDE,false,true);
    break;
  }
}
}
