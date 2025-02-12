#include "JointLimit.h"
#include <ConvexHull/Barrier.h>
#include <Utils/RotationUtils.h>

namespace PHYSICSMOTION {
//energy model
JointLimit::T JointLimit::energyExp(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE,T tol) {
  T E=0,val;
  //if user wants joint limit, X-axis must be locked!
  if(J._limits.array().isFinite().any()) {
    ASSERT_MSG(isfinite(J._limits(2,0)) && J._limits(2,0)>0 &&
               J._limits(0,0)==0 && J._limits(1,0)==0,
               "Expoential joint must be locked in X-axis!")
  } else return E;
  if(isfinite(J._limits(2,1)) && !isfinite(J._limits(2,2))) {
    ASSERT_MSG(isfinite(J._limits(0,1)) && -J._limits(0,1)==J._limits(1,1),"Dual cone joint limits must be symmetric in Y-axis!")
    Vec2T DY;
    Mat2T DDY;
    T y=getASinY(x[J._offDOF+1],x[J._offDOF+2],&DY,&DDY,tol);
    if(isfinite(J._limits(0,1)) && y<J._limits(0,1)) {
      val=y-J._limits(0,1);
      if(DE)
        DE->template segment<2>(J._offDOF+1)+=val*DY*J._limits(2,1);
      if(DDE)
        DDE->template block<2,2>(J._offDOF+1,J._offDOF+1)+=(DY*DY.transpose()+val*DDY)*J._limits(2,1);
      E+=val*val*J._limits(2,1)/2;
    } else if(isfinite(J._limits(1,1)) && y>J._limits(1,1)) {
      val=y-J._limits(1,1);
      if(DE)
        DE->template segment<2>(J._offDOF+1)+=val*DY*J._limits(2,1);
      if(DDE)
        DDE->template block<2,2>(J._offDOF+1,J._offDOF+1)+=(DY*DY.transpose()+val*DDY)*J._limits(2,1);
      E+=val*val*J._limits(2,1)/2;
    }
  } else if(!isfinite(J._limits(2,1)) && isfinite(J._limits(2,2))) {
    ASSERT_MSG(isfinite(J._limits(0,2)) && -J._limits(0,2)==J._limits(1,2),"Dual cone joint limits must be symmetric in Z-axis!")
    Vec2T DZ;
    Mat2T DDZ;
    T z=getASinZ(x[J._offDOF+1],x[J._offDOF+2],&DZ,&DDZ,tol);
    if(isfinite(J._limits(0,2)) && z<J._limits(0,2)) {
      val=z-J._limits(0,2);
      if(DE)
        DE->template segment<2>(J._offDOF+1)+=val*DZ*J._limits(2,2);
      if(DDE)
        DDE->template block<2,2>(J._offDOF+1,J._offDOF+1)+=(DZ*DZ.transpose()+val*DDZ)*J._limits(2,2);
      E+=val*val*J._limits(2,2)/2;
    } else if(isfinite(J._limits(1,2)) && z>J._limits(1,2)) {
      val=z-J._limits(1,2);
      if(DE)
        DE->template segment<2>(J._offDOF+1)+=val*DZ*J._limits(2,2);
      if(DDE)
        DDE->template block<2,2>(J._offDOF+1,J._offDOF+1)+=(DZ*DZ.transpose()+val*DDZ)*J._limits(2,2);
      E+=val*val*J._limits(2,2)/2;
    }
  } else if(isfinite(J._limits(2,1))) {
    ASSERT_MSG(J._limits(2,1)==J._limits(2,2),"Expoential joint does not support different limit stiffness in YZ-axis!")
    ASSERT_MSG(!isfinite(J._limits(0,1)),"Expoential joint limits must be symmetric in Y-axis!")
    ASSERT_MSG(!isfinite(J._limits(0,2)),"Expoential joint limits must be symmetric in Z-axis!")
    Vec2T DY,DZ,D;
    Mat2T DDY,DDZ,DD;
    T y=getSwingAngleY(x[J._offDOF+1],x[J._offDOF+2],&DY,&DDY,tol)/J._limits(1,1);
    T z=getSwingAngleZ(x[J._offDOF+1],x[J._offDOF+2],&DZ,&DDZ,tol)/J._limits(1,2);
    T szSqr=y*y+z*z+tol,sz=sqrt(szSqr);
    if(sz>1) {
      E+=(sz-1)*(sz-1)*J._limits(2,1)/2;
      D=(DY*(y/J._limits(1,1))+DZ*(z/J._limits(1,2)))/sz;
      if(DE)
        DE->template segment<2>(J._offDOF+1)+=D*(sz-1)*J._limits(2,1);
      if(DDE) {
        DD=(DDY*(y/J._limits(1,1))+DDZ*(z/J._limits(1,2)))/sz;
        DD+=(DY*DY.transpose())*(1-y*y/szSqr)/(J._limits(1,1)*J._limits(1,1)*sz);
        DD+=(DZ*DZ.transpose())*(1-z*z/szSqr)/(J._limits(1,2)*J._limits(1,2)*sz);
        DD-=(DY*DZ.transpose()+DZ*DY.transpose())*(y*z)/(J._limits(1,1)*J._limits(1,2)*sz*szSqr);
        DDE->template block<2,2>(J._offDOF+1,J._offDOF+1)+=(DD*(sz-1)+D*D.transpose())*J._limits(2,1);
      }
    }
  }
  return E;
}
JointLimit::T JointLimit::energySimple(VecCM x,const Joint& J,int c,int nrD,Vec* DE,MatT* DDE) {
  T E=0,val;
  if(isfinite(J._limits(2,c)) && J._limits(2,c)>0) {
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && J._limits(0,c)==J._limits(1,c)) {
      //this will be locked, no need to do anything
    } else if(isfinite(J._limits(0,c)) && x[J._offDOF+c]<J._limits(0,c)) {
      //lower limited
      val=x[J._offDOF+c]-J._limits(0,c);
      E+=val*val*J._limits(2,c)/2;
      if(DE)
        DE->coeffRef(J._offDOF+c)+=val*J._limits(2,c);
      if(DDE)
        DDE->diagonal()[J._offDOF+c]+=J._limits(2,c);
    } else if(isfinite(J._limits(1,c)) && x[J._offDOF+c]>J._limits(1,c)) {
      //upper limited
      val=x[J._offDOF+c]-J._limits(1,c);
      E+=val*val*J._limits(2,c)/2;
      if(DE)
        DE->coeffRef(J._offDOF+c)+=val*J._limits(2,c);
      if(DDE)
        DDE->diagonal()[J._offDOF+c]+=J._limits(2,c);
    }
  }
  return E;
}
template <typename Barrier>
bool JointLimit::energySimple(VecCM x,const Joint& J,int c,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p) {
  T D,DD;
  if(isfinite(J._limits(2,c)) && J._limits(2,c)>0) {
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && J._limits(0,c)==J._limits(1,c)) {
      //this will be locked, no need to do anything
    } else if(isfinite(J._limits(0,c))) {
      //lower limited
      if(x[J._offDOF+c]<=J._limits(0,c))
        return false;
      E+=p.template eval<T>(x[J._offDOF+c]-J._limits(0,c),DE?&D:NULL,DDE?&DD:NULL,0,J._limits(2,c));
      if(DE)
        DE->coeffRef(J._offDOF+c)+=D;
      if(DDE)
        DDE->diagonal()[J._offDOF+c]+=DD;
    } else if(isfinite(J._limits(1,c))) {
      //upper limited
      if(x[J._offDOF+c]>=J._limits(1,c))
        return false;
      E+=p.template eval<T>(J._limits(1,c)-x[J._offDOF+c],DE?&D:NULL,DDE?&DD:NULL,0,J._limits(2,c));
      if(DE)
        DE->coeffRef(J._offDOF+c)-=D;
      if(DDE)
        DDE->diagonal()[J._offDOF+c]-=DD;
    }
  }
  return true;
}
JointLimit::T JointLimit::energyBall(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE,T tol) {
  T k=0,lmt=0,szSqr,sz,E=0;
  Eigen::Matrix<T,-1,1,0,3,1> coeff,D;
  coeff.setZero(nrD);
  for(int c=0; c<nrD; c++) {
    //ball
    if(!isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(2,c)>0 && J._limits(1,c)>0) {
      ASSERT_MSG(k==0 || k==J._limits(2,c),"We only support ball-shaped translational joint limit!");
      ASSERT_MSG(lmt==0 || lmt==J._limits(1,c),"We only support ball-shaped translational joint limit!");
      k=J._limits(2,c);
      lmt=J._limits(1,c);
      coeff[c]=1;
    }
    //lock
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(0,c)<=J._limits(1,c))
      E+=energySimple(x,J,c,nrD,DE,DDE);
  }
  szSqr=(x.segment(J._offDOF,nrD).array()*coeff.array()).matrix().squaredNorm()+tol;
  sz=sqrt(szSqr);
  if(sz>lmt) {
    D=(x.segment(J._offDOF,nrD).array()*coeff.array()*coeff.array()).matrix()/sz;
    if(DE)
      DE->segment(J._offDOF,nrD)+=(sz-lmt)*k*D;
    if(DDE) {
      DDE->template block(J._offDOF,J._offDOF,nrD,nrD)+=k*D*D.transpose();
      DDE->template block(J._offDOF,J._offDOF,nrD,nrD)+=(sz-lmt)*k*(coeff.array()*coeff.array()/sz).matrix().asDiagonal();
      DDE->template block(J._offDOF,J._offDOF,nrD,nrD)-=(sz-lmt)*k*(D*D.transpose())/sz;
    }
    E+=(sz-lmt)*(sz-lmt)*k/2;
  }
  return E;
}
template <typename Barrier>
bool JointLimit::energyBall(VecCM x,const Joint& J,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p,T tol) {
  T k=0,lmt=0,szSqr,sz,d,dd;
  Eigen::Matrix<T,-1,1,0,3,1> coeff,D;
  coeff.setZero(nrD);
  for(int c=0; c<nrD; c++) {
    //ball
    if(!isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(2,c)>0 && J._limits(1,c)>0) {
      ASSERT_MSG(k==0 || k==J._limits(2,c),"We only support ball-shaped translational joint limit!");
      ASSERT_MSG(lmt==0 || lmt==J._limits(1,c),"We only support ball-shaped translational joint limit!");
      k=J._limits(2,c);
      lmt=J._limits(1,c);
      coeff[c]=1;
    }
    //lock
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(0,c)<=J._limits(1,c)) {
      if(!energySimple(x,J,c,nrD,E,DE,DDE,p))
        return false;
    }
  }
  szSqr=(x.segment(J._offDOF,nrD).array()*coeff.array()).matrix().squaredNorm()+tol;
  sz=sqrt(szSqr);
  if(sz>=lmt)
    return false;
  E+=p.template eval<T>(lmt-sz,DE?&d:NULL,DDE?&dd:NULL,0,k);
  D=(x.segment(J._offDOF,nrD).array()*coeff.array()*coeff.array()).matrix()/sz;
  if(DE)
    DE->segment(J._offDOF,nrD)-=d*D;
  if(DDE) {
    DDE->template block(J._offDOF,J._offDOF,nrD,nrD)+=dd*D*D.transpose();
    DDE->template block(J._offDOF,J._offDOF,nrD,nrD)-=d*(coeff.array()*coeff.array()/sz).matrix().asDiagonal();
    DDE->template block(J._offDOF,J._offDOF,nrD,nrD)+=d*(D*D.transpose())/sz;
  }
  return true;
}
JointLimit::T JointLimit::energy(VecCM x,const Joint& J,int nrD,Vec* DE,MatT* DDE) {
  if(!J.isRotational())
    return energyBall(x,J,nrD,DE,DDE);
  else if(J._typeJoint==Joint::ROT_3D_EXP)
    return energyExp(x,J,nrD,DE,DDE);
  else {
    T E=0;
    for(int c=0; c<nrD; c++)
      E+=energySimple(x,J,c,nrD,DE,DDE);
    return E;
  }
}
template <typename Barrier>
bool JointLimit::energy(VecCM x,const Joint& J,int nrD,T& E,Vec* DE,MatT* DDE,const Barrier& p) {
  if(!J.isRotational()) {
    if(!energyBall(x,J,nrD,E,DE,DDE,p))
      return false;
  } else if(J._typeJoint==Joint::ROT_3D_EXP) {
    E+=energyExp(x,J,nrD,DE,DDE);
  } else {
    for(int c=0; c<nrD; c++)
      if(!energySimple(x,J,c,nrD,E,DE,DDE,p))
        return false;
  }
  return true;
}
//constraint model
void JointLimit::constraintExp(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol) {
  //if user wants joint limit, X-axis must be locked!
  if(J._limits.array().isFinite().any()) {
    ASSERT_MSG(isfinite(J._limits(2,0)) && J._limits(2,0)>0 &&
               J._limits(0,0)==0 && J._limits(1,0)==0,
               "Expoential joint must be locked in X-axis!")
  } else return;
  if(isfinite(J._limits(2,1)) && !isfinite(J._limits(2,2))) {
    ASSERT_MSG(isfinite(J._limits(0,1)) && -J._limits(0,1)==J._limits(1,1),"Dual cone joint limits must be symmetric in Y-axis!")
    //allocate constraint
    if((int)Css.size()<=off) {
      Css.push_back(Simulator::XPBDConstraint());
      Css[off]._JC.setZero(nrDOF);
    }
    //build constraint
    Simulator::XPBDConstraint& C=Css[off++];
    C._JC.setZero();
    C._C=0;
    Vec2T DY;
    T y=getASinY(x[J._offDOF+1],x[J._offDOF+2],&DY,NULL,tol);
    if(isfinite(J._limits(0,1)) && y<J._limits(0,1))
      C._C=y-J._limits(0,1);
    else if(isfinite(J._limits(1,1)) && y>J._limits(1,1))
      C._C=y-J._limits(1,1);
    C._JC.template segment<2>(J._offDOF+1)=DY;
    C._alpha=J._limits(2,1);
  } else if(!isfinite(J._limits(2,1)) && isfinite(J._limits(2,2))) {
    ASSERT_MSG(isfinite(J._limits(0,2)) && -J._limits(0,2)==J._limits(1,2),"Dual cone joint limits must be symmetric in Z-axis!")
    //allocate constraint
    if((int)Css.size()<=off) {
      Css.push_back(Simulator::XPBDConstraint());
      Css[off]._JC.setZero(nrDOF);
    }
    //build constraint
    Simulator::XPBDConstraint& C=Css[off++];
    C._JC.setZero();
    C._C=0;
    Vec2T DZ;
    T z=getASinZ(x[J._offDOF+1],x[J._offDOF+2],&DZ,NULL,tol);
    if(isfinite(J._limits(0,2)) && z<J._limits(0,2))
      C._C=z-J._limits(0,2);
    else if(isfinite(J._limits(1,2)) && z>J._limits(1,2))
      C._C=z-J._limits(1,2);
    C._JC.template segment<2>(J._offDOF+1)=DZ;
    C._alpha=J._limits(2,2);
  } else if(isfinite(J._limits(2,1))) {
    ASSERT_MSG(J._limits(2,1)==J._limits(2,2),"Expoential joint does not support different limit stiffness in YZ-axis!")
    ASSERT_MSG(!isfinite(J._limits(0,1)),"Expoential joint limits must be symmetric in Y-axis!")
    ASSERT_MSG(!isfinite(J._limits(0,2)),"Expoential joint limits must be symmetric in Z-axis!")
    //allocate constraint
    if((int)Css.size()<=off) {
      Css.push_back(Simulator::XPBDConstraint());
      Css[off]._JC.setZero(nrDOF);
    }
    //build constraint
    Simulator::XPBDConstraint& C=Css[off++];
    C._JC.setZero();
    C._C=0;
    Vec2T DY,DZ,D;
    T y=getSwingAngleY(x[J._offDOF+1],x[J._offDOF+2],&DY,NULL,tol)/J._limits(1,1);
    T z=getSwingAngleZ(x[J._offDOF+1],x[J._offDOF+2],&DZ,NULL,tol)/J._limits(1,2);
    T szSqr=y*y+z*z+tol,sz=sqrt(szSqr);
    D=(DY*(y/J._limits(1,1))+DZ*(z/J._limits(1,2)))/sz;
    if(sz>1)
      C._C=sz-1;
    C._JC.template segment<2>(J._offDOF+1)=D;
    C._alpha=J._limits(2,1);
  }
}
void JointLimit::constraintSimple(VecCM x,const Joint& J,int c,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off) {
  if(isfinite(J._limits(2,c)) && J._limits(2,c)>0) {
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && J._limits(0,c)==J._limits(1,c)) {
      //this will be locked, no need to do anything
      return;
    }
    //allocate constraint
    if((int)Css.size()<=off) {
      Css.push_back(Simulator::XPBDConstraint());
      Css[off]._JC.setZero(nrDOF);
    }
    //build constraint
    Simulator::XPBDConstraint& C=Css[off++];
    C._JC.setZero();
    C._C=0;
    C._isLinear=true;
    T xc=x[J._offDOF+c];
    if(isfinite(J._limits(0,c)) && xc<J._limits(0,c)) {
      C._C=xc-J._limits(0,c);
      C._CL=J._limits(0,c);
    } else if(isfinite(J._limits(1,c)) && xc>J._limits(1,c)) {
      C._C=xc-J._limits(1,c);
      C._CU=J._limits(1,c);
    }
    C._JC.setUnit(C._JC.size(),J._offDOF+c);
    C._alpha=J._limits(2,c);
  }
}
void JointLimit::constraintBall(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off,T tol) {
  T k=0,lmt=0,szSqr,sz;
  Eigen::Matrix<T,-1,1,0,3,1> coeff;
  coeff.setZero(nrD);
  for(int c=0; c<nrD; c++) {
    //ball
    if(!isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(2,c)>0 && J._limits(1,c)>0) {
      ASSERT_MSG(k==0 || k==J._limits(2,c),"We only support ball-shaped translational joint limit!");
      ASSERT_MSG(lmt==0 || lmt==J._limits(1,c),"We only support ball-shaped translational joint limit!");
      k=J._limits(2,c);
      lmt=J._limits(1,c);
      coeff[c]=1;
    }
    //lock
    if(isfinite(J._limits(0,c)) && isfinite(J._limits(1,c)) && isfinite(J._limits(2,c)) && J._limits(2,c)>0 && J._limits(0,c)<=J._limits(1,c))
      constraintSimple(x,J,c,nrD,nrDOF,Css,off);
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
  szSqr=(x.segment(J._offDOF,nrD).array()*coeff.array()).matrix().squaredNorm()+tol;
  sz=sqrt(szSqr);
  if(sz>lmt)
    C._C=sz-lmt;
  C._JC.segment(J._offDOF,nrD)=(x.segment(J._offDOF,nrD).array()*coeff.array()*coeff.array()).matrix()/sz;
  C._alpha=k;
}
void JointLimit::constraint(VecCM x,const Joint& J,int nrD,int nrDOF,std::vector<Simulator::XPBDConstraint>& Css,int& off) {
  if(!J.isRotational())
    constraintBall(x,J,nrD,nrDOF,Css,off);
  else if(J._typeJoint==Joint::ROT_3D_EXP)
    constraintExp(x,J,nrD,nrDOF,Css,off);
  else {
    for(int c=0; c<nrD; c++)
      constraintSimple(x,J,c,nrD,nrDOF,Css,off);
  }
}
//debug
void JointLimit::debugSwing(T tol) {
  Mat2T DD;
  T val,val2;
  Vec2T x,dx,D,D2;
  DEFINE_NUMERIC_DELTA_T(T)
  x=Vec2T::Random();
  dx=Vec2T::Random();
  //ASinY
  val=getASinY(x[0],x[1],&D,&DD,tol);
  val2=getASinY(x[0]+dx[0]*DELTA,x[1]+dx[1]*DELTA,&D2,NULL,tol);
  DEBUG_GRADIENT("SinY-G",D.dot(dx),D.dot(dx)-(val2-val)/DELTA)
  DEBUG_GRADIENT("SinY-H",(DD*dx).norm(),(DD*dx-(D2-D)/DELTA).norm())
  //ASinZ
  val=getASinZ(x[0],x[1],&D,&DD,tol);
  val2=getASinZ(x[0]+dx[0]*DELTA,x[1]+dx[1]*DELTA,&D2,NULL,tol);
  DEBUG_GRADIENT("SinZ-G",D.dot(dx),D.dot(dx)-(val2-val)/DELTA)
  DEBUG_GRADIENT("SinZ-H",(DD*dx).norm(),(DD*dx-(D2-D)/DELTA).norm())
  //SwingY
  val=getSwingAngleY(x[0],x[1],&D,&DD,tol);
  val2=getSwingAngleY(x[0]+dx[0]*DELTA,x[1]+dx[1]*DELTA,&D2,NULL,tol);
  DEBUG_GRADIENT("SwingY-G",D.dot(dx),D.dot(dx)-(val2-val)/DELTA)
  DEBUG_GRADIENT("SwingY-H",(DD*dx).norm(),(DD*dx-(D2-D)/DELTA).norm())
  //SwingZ
  val=getSwingAngleZ(x[0],x[1],&D,&DD,tol);
  val2=getSwingAngleZ(x[0]+dx[0]*DELTA,x[1]+dx[1]*DELTA,&D2,NULL,tol);
  DEBUG_GRADIENT("SwingZ-G",D.dot(dx),D.dot(dx)-(val2-val)/DELTA)
  DEBUG_GRADIENT("SwingZ-H",(DD*dx).norm(),(DD*dx-(D2-D)/DELTA).norm())
}
void JointLimit::debugExp(int off,T tol) {
  Vec x;
  Joint J;
  J._offDOF=off;
  J._typeJoint=Joint::ROT_3D_EXP;
  J._limits.setZero(3,J.nrDOF());
  //ASinYZ
  for(int c=1; c<=2; c++) {
    J._limits.setConstant(std::numeric_limits<double>::infinity());
    J._limits.col(0).setZero();
    J._limits(2,0)=1;
    while(true) {
      J._limits.col(c).setRandom();
      if(J._limits(1,c)<=0)
        continue;
      if(J._limits(2,c)<=0)
        continue;
      J._limits(0,c)=-J._limits(1,c);
      break;
    }
    T angle[3];
    debug("ASin"+std::to_string(c)+"L",x,J,[&](const Vec& x) {
      angle[1]=getASinY(x[J._offDOF+1],x[J._offDOF+2],NULL,NULL,tol);
      angle[2]=getASinZ(x[J._offDOF+1],x[J._offDOF+2],NULL,NULL,tol);
      return angle[c]<J._limits(0,c);
    });
    debug("ASin"+std::to_string(c)+"U",x,J,[&](const Vec& x) {
      angle[1]=getASinY(x[J._offDOF+1],x[J._offDOF+2],NULL,NULL,tol);
      angle[2]=getASinZ(x[J._offDOF+1],x[J._offDOF+2],NULL,NULL,tol);
      return angle[c]>J._limits(1,c);
    });
  }
  //SwingYZ
  J._limits.setConstant(std::numeric_limits<double>::infinity());
  J._limits.col(0).setZero();
  J._limits(2,0)=1;
  J._limits.template block<1,2>(2,1).setConstant(rand()/(double)RAND_MAX);
  J._limits.template block<1,2>(1,1).setConstant(rand()/(double)RAND_MAX);
  debug("Swing",x,J,[&](const Vec& x) {
    return true;
  });
}
void JointLimit::debugSimple(int off) {
  Vec x;
  Joint J;
  J._offDOF=off;
  J._typeJoint=Joint::ROT_3D_XYZ;
  J._limits.setZero(3,J.nrDOF());
  for(int c=0; c<J.nrDOF(); c++)
    while(true) {
      J._limits.col(c).setRandom();
      if(J._limits(0,c)>J._limits(1,c))
        continue;
      if(J._limits(2,c)<=0)
        continue;
      break;
    }
  for(int c=0; c<J.nrDOF(); c++) {
    debug("SimpleL",x,J,[&](const Vec& x) {
      return x[J._offDOF+c]<J._limits(0,c);
    });
    debug("SimpleU",x,J,[&](const Vec& x) {
      return x[J._offDOF+c]>J._limits(1,c);
    });
  }
}
template <typename Barrier>
void JointLimit::debugSimple(int off,const Barrier& p) {
  Vec x;
  Joint J;
  J._offDOF=off;
  J._typeJoint=Joint::ROT_3D_XYZ;
  J._limits.setZero(3,J.nrDOF());
  for(int c=0; c<J.nrDOF(); c++)
    while(true) {
      J._limits.col(c).setRandom();
      if(J._limits(0,c)>J._limits(1,c))
        continue;
      if(J._limits(2,c)<=0)
        continue;
      break;
    }
  for(int c=0; c<J.nrDOF(); c++)
    debugHard("SimpleLU",x,J,p);
}
void JointLimit::debugBall(int off) {
  Vec x;
  Joint J;
  J._offDOF=off;
  J._typeJoint=Joint::TRANS_3D;
  J._limits.resize(3,J.nrDOF());
  for(int d=0; d<J.nrDOF()+1; d++) {
    J._limits.setConstant(std::numeric_limits<double>::infinity());
    J._limits.row(1).setConstant(rand()/(double)RAND_MAX);
    J._limits.row(2).setConstant(rand()/(double)RAND_MAX);
    if(d<J.nrDOF())
      J._limits(0,d)=J._limits(1,d)=0;
    debug("BallLock"+std::to_string(d),x,J,[&](const Vec& x) {
      return true;
    });
  }
}
template <typename Barrier>
void JointLimit::debugBall(int off,const Barrier& p) {
  Vec x;
  Joint J;
  J._offDOF=off;
  J._typeJoint=Joint::TRANS_3D;
  J._limits.resize(3,J.nrDOF());
  for(int d=0; d<J.nrDOF()+1; d++) {
    J._limits.setConstant(std::numeric_limits<double>::infinity());
    J._limits.row(1).setConstant(rand()/(double)RAND_MAX);
    J._limits.row(2).setConstant(rand()/(double)RAND_MAX);
    if(d<J.nrDOF())
      J._limits(0,d)=J._limits(1,d)=0;
    debugHard("BallLock"+std::to_string(d),x,J,p);
  }
}
void JointLimit::debug(const std::string& name,const Vec& x,const Joint& J,std::function<bool(const Vec&)> func) {
  DEFINE_NUMERIC_DELTA_T(T)
  while(true) {
    //debug energy
    Vec x=Vec::Random(J._offDOF+J.nrDOF());
    Vec dx=Vec::Random(J._offDOF+J.nrDOF()),x2=x+dx*DELTA;
    Vec DE0=Vec::Random(x.size()),DE=DE0,DE2=DE0;
    MatT DDE0=MatT::Random(x.size(),x.size()),DDE=DDE0;
    if(!func(x))
      continue;
    T E=energy(mapCV(x),J,J.nrDOF(),&DE,&DDE);
    T E2=energy(mapCV(x2),J,J.nrDOF(),&DE2,NULL);
    if(E==0)
      continue;
    DE-=DE0;
    DDE-=DDE0;
    DE2-=DE0;
    DEBUG_GRADIENT("JointLimit"+name+"-G",DE.dot(dx),DE.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("JointLimit"+name+"-H",(DDE*dx).norm(),(DDE*dx-(DE2-DE)/DELTA).norm())
    //debug constraint
    int off=0;
    Vec DERef=Vec::Zero(x.size());
    std::vector<Simulator::XPBDConstraint> Css;
    constraint(mapCV(x),J,J.nrDOF(),x.size(),Css,off);
    for(int k=0; k<off; k++)
      DERef+=Css[k]._JC*Css[k]._C*Css[k]._alpha;
    DEBUG_GRADIENT("JointLimit"+name+"-CG",DE.norm(),(DE-DERef).norm())
    break;
  }
}
template <typename Barrier>
void JointLimit::debugHard(const std::string& name,const Vec& x,const Joint& J,const Barrier& p) {
  DEFINE_NUMERIC_DELTA_T(T)
  while(true) {
    //debug energy
    Vec x=Vec::Random(J._offDOF+J.nrDOF());
    Vec dx=Vec::Random(J._offDOF+J.nrDOF()),x2=x+dx*DELTA;
    Vec DE0=Vec::Random(x.size()),DE=DE0,DE2=DE0;
    MatT DDE0=MatT::Random(x.size(),x.size()),DDE=DDE0;
    T E=x[0],E2=E;
    if(!energy(mapCV(x),J,J.nrDOF(),E,&DE,&DDE,p))
      continue;
    if(!energy(mapCV(x2),J,J.nrDOF(),E2,&DE2,NULL,p))
      continue;
    DE-=DE0;
    DDE-=DDE0;
    DE2-=DE0;
    DEBUG_GRADIENT("JointLimit"+name+"-HardG",DE.dot(dx),DE.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("JointLimit"+name+"-HardH",(DDE*dx).norm(),(DDE*dx-(DE2-DE)/DELTA).norm())
    break;
  }
}
//helper
JointLimit::T JointLimit::getASinY(T y,T z,Vec2T* DY,Mat2T* DDY,T tol) {
  T asy;
  T t = sqrt(y*y+z*z);
  if(t < tol) {
    T x0 = pow(y, 2);
    T x1 = (1.0/6.0)*x0 - 1;
    asy = -x1*z;
    if(DY) {
      (*DY)[0] = -1.0/3.0*y*z;
      (*DY)[1] = -1.0/5.0*x0*pow(z, 2) - x1;
    }
    if(DDY) {
      x0 = pow(y, 2.0);
      x1 = -2.0/5.0*y*pow(z, 2) - 1.0/3.0*y;
      (*DDY)(0,0) = z*((1.0/10.0)*x0 - 1.0/3.0);
      (*DDY)(0,1) = x1;
      (*DDY)(1,0) = x1;
      (*DDY)(1,1) = -2.0/5.0*x0*z;
    }
  } else {
    T x0 = pow(y, 2);
    T x1 = pow(z, 2);
    T x2 = x0 + x1;
    T x3 = sqrt(x2);
    T x4 = sin(x3);
    T x5 = cos(x3);
    T x6 = x3*x5;
    T x7 = 1/(pow(x2, 3.0/2.0)*sqrt((x0 + x1*pow(x5, 2))/x2));
    asy = asin(x4*z/x3);
    if(DY) {
      (*DY)[0] = x7*y*z*(-x4 + x6);
      (*DY)[1] = x7*(x0*x4 + x1*x6);
    }
    if(DDY) {
      x0 = pow(y, 4.0);
      x1 = 2*x0;
      x2 = pow(y, 2.0);
      x3 = pow(z, 2.0);
      x4 = x2*x3;
      x5 = pow(z, 4.0);
      x6 = x2 + x3;
      x7 = sqrt(x6);
      T x8 = cos(x7);
      T x9 = pow(x8, 2);
      T x10 = x5*x9;
      T x11 = pow(x6, 5.0/2.0)*x8;
      T x12 = x0*x3;
      T x13 = x3*x9;
      T x14 = sin(x7);
      T x15 = x14*pow(x6, 2);
      T x16 = x13 + x2;
      T x17 = 1/(x16*pow(x6, 9.0/2.0)*sqrt(x16/x6));
      T x18 = x17*z;
      T x19 = x4 + x5;
      T x20 = -x17*y*(x11*(-x0 + x10 + x19) + x15*(x0 - 2*x10 + x12 + x2*x5 - x4));
      (*DDY)(0,0) = -x18*(x11*(x1 - x10 + x4) + x15*(-x1 + x10 + x12 - x13*x2 + pow(y, 6)));
      (*DDY)(0,1) = x20;
      (*DDY)(1,0) = x20;
      (*DDY)(1,1) = x18*x2*(x11*(-pow(x14, 2)*x3 + 3*x2 + 3*x3) - x15*(3*x13 + x19 + x2*x9 + 2*x2));
    }
  }
  return asy;
}
JointLimit::T JointLimit::getASinZ(T y,T z,Vec2T* DZ,Mat2T* DDZ,T tol) {
  T asz;
  T t = sqrt(y*y+z*z);
  if(t < tol) {
    T x0 = pow(z, 2);
    T x1 = (1.0/6.0)*x0 - 1;
    asz = x1*y;
    if(DZ) {
      (*DZ)[0] = (1.0/5.0)*x0*pow(y, 2) + x1;
      (*DZ)[1] = (1.0/3.0)*y*z;
    }
    if(DDZ) {
      x0 = pow(z, 2.0);
      x1 = (2.0/5.0)*pow(y, 2)*z + (1.0/3.0)*z;
      (*DDZ)(0,0) = (2.0/5.0)*x0*y;
      (*DDZ)(0,1) = x1;
      (*DDZ)(1,0) = x1;
      (*DDZ)(1,1) = y*(1.0/3.0 - 1.0/10.0*x0);
    }
  } else {
    T x0 = pow(y, 2);
    T x1 = pow(z, 2);
    T x2 = x0 + x1;
    T x3 = sqrt(x2);
    T x4 = sin(x3);
    T x5 = cos(x3);
    T x6 = x3*x5;
    T x7 = 1/(pow(x2, 3.0/2.0)*sqrt((x0*pow(x5, 2) + x1)/x2));
    asz = -asin(x4*y/x3);
    if(DZ) {
      (*DZ)[0] = -x7*(x0*x6 + x1*x4);
      (*DZ)[1] = x7*y*z*(x4 - x6);
    }
    if(DDZ) {
      x0 = pow(z, 2.0);
      x1 = pow(y, 2.0);
      x2 = x0 + x1;
      x3 = sqrt(x2);
      x4 = sin(x3);
      x5 = cos(x3);
      x6 = pow(x2, 5.0/2.0)*x5;
      x7 = pow(x5, 2.0);
      T x8 = x1*x7;
      T x9 = x0*x1;
      T x10 = pow(y, 4);
      T x11 = x10 + x9;
      T x12 = pow(x2, 2)*x4;
      T x13 = x0 + x8;
      T x14 = 1/(x13*pow(x2, 9.0/2.0)*sqrt(x13/x2));
      T x15 = x14*y;
      T x16 = pow(z, 4);
      T x17 = x10*x7;
      T x18 = x1*x16;
      T x19 = x14*z*(x12*(x0*x10 + x16 - 2*x17 + x18 - x9) + x6*(x11 - x16 + x17));
      T x20 = 2*x16;
      (*DDZ)(0,0) = x0*x15*(x12*(x0*x7 + 2*x0 + x11 + 3*x8) + x6*(-3*x0 + x1*pow(x4, 2) - 3*x1));
      (*DDZ)(0,1) = x19;
      (*DDZ)(1,0) = x19;
      (*DDZ)(1,1) = x15*(x12*(-x0*x8 + x17 + x18 - x20 + pow(z, 6)) + x6*(-x17 + x20 + x9));
    }
  }
  return asz;
}
JointLimit::T JointLimit::getSwingAngleY(T y, T z,Vec2T* DY,Mat2T* DDY,T tol) {
  T result;
  T t = sqrt(y*y+z*z);
  if(t < tol) {
    if(DY) {
      (*DY)[0] = pow(z, 2)*(1.0/48.0 - 1.0/1280.0*pow(y, 2)) + 1;
      (*DY)[1] = (1.0/24.0)*y*z;
    }
    if(DDY) {
      T x0 = y*pow(z, 2);
      T x1 = z*(1.0/24.0 - 1.0/640.0*pow(y, 2));
      (*DDY)(0,0) = -1.0/640.0*x0;
      (*DDY)(0,1) = x1;
      (*DDY)(1,0) = x1;
      (*DDY)(1,1) = (1.0/160.0)*x0 + (1.0/24.0)*y;
    }
  } else {
    if(DY) {
      T x0 = pow(y, 2);
      T x1 = pow(z, 2);
      T x2 = sqrt(x0 + x1);
      T x3 = (1.0/2.0)*x2;
      T x4 = 2*sin(x3);
      T x5 = 2/(x2*(2*x0 + x1*cos(x3) + x1));
      (*DY)[0] = x5*(x0*x2 + x1*x4);
      (*DY)[1] = x5*y*z*(x2 - x4);
    }
    if(DDY) {
      T x0 = pow(y, 2);
      T x1 = pow(z, 2);
      T x2 = 2*x1;
      T x3 = x0 + x1;
      T x4 = sqrt(x3);
      T x5 = (1.0/2.0)*x4;
      T x6 = cos(x5);
      T x7 = x1*x6;
      T x8 = pow(x3, 5.0/2.0);
      T x9 = 2*x0;
      T x10 = x1 + x7 + x9;
      T x11 = 2*x10;
      T x12 = sin(x5);
      T x13 = x1*x12;
      T x14 = x0*x4 + 2*x13;
      T x15 = x11*x14*pow(x3, 2);
      T x16 = x13 - 8*x4;
      T x17 = x14*x8;
      T x18 = pow(x3, 7.0/2.0);
      T x19 = pow(x10, -2);
      T x20 = x19/x18;
      T x21 = -x13 + 4*x4*(x6 + 1);
      T x22 = x10*x9;
      T x23 = pow(x3, 3);
      T x24 = x23*(x6 - 1);
      T x25 = -2*x12 + x4;
      T x26 = -2*x10*x18*x25;
      T x27 = x25*x8;
      T x28 = x19/pow(x3, 4);
      T x29 = x10*x2;
      (*DDY)(0,0) = x20*y*(x11*x8*(3*x0 + x2 + x7) - x15 + x16*x17);
      (*DDY)(0,1) = x20*z*(2*x10*x8*(x0 + 4*x12*x4 + x7) - x15 - x17*x21);
      (*DDY)(1,0) = x28*z*(x0*x16*x23*x25 - x22*x24 - x22*x27 - x26);
      (*DDY)(1,1) = x28*y*(-x1*x21*x23*x25 - x24*x29 - x26 - x27*x29);
    }
  }
  {
    T x0 = sqrt(pow(y, 2) + pow(z, 2));
    T x1 = (1.0/2.0)*x0;
    result = 4*atan2(y*sin(x1), x0*(cos(x1) + 1));
    return result;
  }
}
JointLimit::T JointLimit::getSwingAngleZ(T y, T z,Vec2T* DZ,Mat2T* DDZ,T tol) {
  T result;
  T t = sqrt(y*y+z*z);
  if(t < tol) {
    if(DZ) {
      (*DZ)[0] = (1.0/24.0)*y*z;
      (*DZ)[1] = pow(y, 2)*(1.0/48.0 - 1.0/1280.0*pow(z, 2)) + 1;
    }
    if(DDZ) {
      T x0 = pow(y, 2)*z;
      T x1 = y*(1.0/24.0 - 1.0/640.0*pow(z, 2));
      (*DDZ)(0,0) = (1.0/160.0)*x0 + (1.0/24.0)*z;
      (*DDZ)(0,1) = x1;
      (*DDZ)(1,0) = x1;
      (*DDZ)(1,1) = -1.0/640.0*x0;
    }
  } else {
    if(DZ) {
      T x0 = pow(y, 2);
      T x1 = pow(z, 2);
      T x2 = sqrt(x0 + x1);
      T x3 = (1.0/2.0)*x2;
      T x4 = 2*sin(x3);
      T x5 = 2/(x2*(x0*cos(x3) + x0 + 2*x1));
      (*DZ)[0] = x5*y*z*(x2 - x4);
      (*DZ)[1] = x5*(x0*x4 + x1*x2);
    }
    if(DDZ) {
      T x0 = pow(y, 2);
      T x1 = pow(z, 2);
      T x2 = 2*x1;
      T x3 = x0 + x1;
      T x4 = sqrt(x3);
      T x5 = (1.0/2.0)*x4;
      T x6 = cos(x5);
      T x7 = x0*x6;
      T x8 = x0 + x2 + x7;
      T x9 = 2*x0;
      T x10 = x8*x9;
      T x11 = pow(x3, 3);
      T x12 = x11*(x6 - 1);
      T x13 = pow(x3, 7.0/2.0);
      T x14 = sin(x5);
      T x15 = 2*x14;
      T x16 = -x15 + x4;
      T x17 = -2*x13*x16*x8;
      T x18 = pow(x3, 5.0/2.0);
      T x19 = x16*x18;
      T x20 = x0*x14;
      T x21 = -x20 + 4*x4*(x6 + 1);
      T x22 = pow(x8, -2);
      T x23 = x22/pow(x3, 4);
      T x24 = x2*x8;
      T x25 = x20 - 8*x4;
      T x26 = x0*x15 + x1*x4;
      T x27 = 2*x8;
      T x28 = x26*x27*pow(x3, 2);
      T x29 = x18*x26;
      T x30 = x22/x13;
      (*DDZ)(0,0) = x23*z*(-x0*x11*x16*x21 - x10*x12 - x10*x19 - x17);
      (*DDZ)(0,1) = x23*y*(x1*x11*x16*x25 - x12*x24 - x17 - x19*x24);
      (*DDZ)(1,0) = x30*y*(2*x18*x8*(x1 + 4*x14*x4 + x7) - x21*x29 - x28);
      (*DDZ)(1,1) = x30*z*(x18*x27*(3*x1 + x7 + x9) + x25*x29 - x28);
    }
  }
  {
    T x0 = sqrt(pow(y, 2) + pow(z, 2));
    T x1 = (1.0/2.0)*x0;
    result = 4*atan2(z*sin(x1), x0*(cos(x1) + 1));
    return result;
  }
}
//instance
template void JointLimit::debugSimple(int off,const CLogx& p);
template void JointLimit::debugBall(int off,const CLogx& p);
}
