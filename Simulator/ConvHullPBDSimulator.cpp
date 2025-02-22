#include "ConvHullPBDSimulator.h"
#include "PBDMatrixSolver.h"
#include "JointLimit.h"
#include "SoftJoint.h"
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>
#include <Articulated/ArticulatedUtils.h>
#include <time.h>

namespace PHYSICSMOTION {
#define USE_CRBA 0
void moveMesh(Joint& J,const Eigen::Matrix<GEOMETRY_SCALAR,-1,1>& X) {
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(J._mesh);
  mesh->moveMesh(X);
  mesh->init(mesh->vss(),mesh->iss(),true);
}
void setMesh(Joint& J,const Eigen::Matrix<GEOMETRY_SCALAR,-1,1>& X) {
  std::shared_ptr<MeshExact> mesh=std::dynamic_pointer_cast<MeshExact>(J._mesh);
  mesh->setMesh(X);
  mesh->init(mesh->vss(),mesh->iss(),true);
}
ConvHullPBDSimulator::ConvHullPBDSimulator(T dt):Simulator(dt),_gTol(1e-4f),_alpha(1e-6f),_epsV(1e-1f),_coefBarrier(1e-16),_hardLimit(false),_maxIt(1e4) {
  _barrier._x0=10000;//0.01
  _sol.reset(new PBDMatrixSolverEigen(_body));
}
void ConvHullPBDSimulator::clearShape() {
  _shapes.clear();
  _obs.clear();
  _contact=NULL;
  _t=0;
}
void ConvHullPBDSimulator::addShape(std::shared_ptr<ShapeExact> shape) {
  _shapes.push_back(shape);
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  shape->getMesh(vss,iss);
  std::shared_ptr<MeshExact> mesh(new MeshExact(vss,iss,true));
  _obs.push_back(GJKPolytope<T>(mesh));
  _contact=NULL;
}
void ConvHullPBDSimulator::setCustomEnergy(std::shared_ptr<CustomPBDEnergy<T>> custom) {
  _custom=custom;
}
void ConvHullPBDSimulator::setArticulatedBody(std::shared_ptr<ArticulatedBody> body) {
  ArticulatedUtils(*body).tessellate(true); //tessellate the mesh
  Simulator::setArticulatedBody(body);
  //setGravity(Vec3T::Zero());
  _JRCF.setZero(3,4*_body->nrJ());
  _pos.reset(*body,Vec::Zero(body->nrDOF()));
  _lastPos.reset(*body,Vec::Zero(body->nrDOF()));
}
ConvHullPBDSimulator::Vec3T ConvHullPBDSimulator::getCentre(int k) {
  return _pos._centre[k];
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdPTarget() {
  return _pos._HThetaPTarget;
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdDTarget() {
  return _pos._HThetaDTarget;
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdL() {
  return _pos._HThetaL;
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdLL() {
  return _pos._HThetaLL;
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdD() {
  return _pos._HThetaD;
}
ConvHullPBDSimulator::MatT ConvHullPBDSimulator::getdPos() {
  return _pos._HPos;
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::getGlobalPoints(bool backward) {
  int nrJ=_body->nrJ();
  Vec ConvexPoints;
  Mat3T Mwt,Mtt,C;
  int nrDB=_body->nrDOF();
  ConvexPoints.setZero(_pos._nrVex*3);
  _pos._HPos.setZero(nrDB,_pos._nrVex*3);
  OMP_PARALLEL_FOR_
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
      for(int c=0; c<(int)local->vss().size(); c++) {
        Vec3T global=(ROT(TRANSI(_pos._info._TM,k))*local->vss()[c].template cast<T>()+CTR(TRANSI(_pos._info._TM,k)));
        for(int j=0; j<3; j++)
          ConvexPoints(j+(_pos._polytopes[k].getVertexId()[0]+c)*3)=global[j];
        if(backward){
          MatX3T HThetaX;
          HThetaX.setZero(nrDB,3);
          C=cross<T>(ROTI(_pos._info._TM,k)*global);
          Mwt=C;
          Mtt=Mat3T::Identity();
          _pos._info.JRCSparse(*_body,k,[&](int row,const Vec3T& JR) {
            parallelAdd<T,3>(HThetaX,row,0,JR.transpose()*Mwt);
          },[&](int row,const Vec3T& JC) {
            parallelAdd<T,3>(HThetaX,row,0,JC.transpose()*Mtt);
          });
          for(int i=0; i<nrDB; i++)
            for(int j=0; j<3; j++)
              parallelAdd(_pos._HPos(i,j+(_pos._polytopes[k].getVertexId()[0]+c)*3),HThetaX(i,j));
        }
      }
    }    
  }
  return ConvexPoints;
}
void ConvHullPBDSimulator::setclass(Eigen::Matrix<int,-1,1> c){
  for(int k=0; k<=_body->nrJ(); k++) {
    //control
    Joint& J=_body->joint(k);
    J._class=c[k];
  }
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::getConvexPoints() {
  int nrJ=_body->nrJ();
  Vec ConvexPoints;
  ConvexPoints.setZero(_pos._nrVex*3);
  OMP_PARALLEL_FOR_
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
      for(int c=0; c<(int)local->vss().size(); c++) {
        for(int j=0; j<3; j++)
          ConvexPoints(j+(_pos._polytopes[k].getVertexId()[0]+c)*3)=(local->vss()[c].template cast<T>())(j);
      }
    }
  }
  return ConvexPoints;
}
void ConvHullPBDSimulator::getJointPosGrad(GradInfo& grad,int k) {
  const Joint& J=_body->joint(k);
  Vec3T centre=Vec3T::Zero();
  Mat3T Mwt,Mtt,C;
  MatX3T HThetaX;
  int nrDB=_body->nrDOF();
  HThetaX.setZero(nrDB,3);
  if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
    int sz=(int)local->vss().size();
    for(int i=0; i<sz; i++) {
      centre+=local->vss()[i].template cast<T>()/sz;
    }
    C=cross<T>(ROTI(grad._info._TM,k)*centre);
    Mwt=C;
    Mtt=Mat3T::Identity();
    grad._info.JRCSparse(*_body,k,[&](int row,const Vec3T& JR) {
      parallelAdd<T,3>(HThetaX,row,0,JR.transpose()*Mwt);
    },[&](int row,const Vec3T& JC) {
      parallelAdd<T,3>(HThetaX,row,0,JC.transpose()*Mtt);
    });
    parallelAdd<T,-1,-1>(grad._HPos,0,0,HThetaX);
  } else std::cout<<"Invalid Joint !"<<std::endl;
}
void ConvHullPBDSimulator::reset() {
  setPos(Vec::Zero(_body->nrDOF()));
  setVel(Vec::Zero(_body->nrDOF()));
  setTime(0);
}
void ConvHullPBDSimulator::resetWithPos(Vec pos) {
  setPos(pos);
  setVel(Vec::Zero(_body->nrDOF()));
  setTime(0);
}
void ConvHullPBDSimulator::step() {
  if(_body->nrDOF()==0)
    return;
  //update kinematic state
  Vec D,DE;
  D=setKinematic(_lastPos._info._xM,_t-_dt);
  if(D!=_lastPos._info._xM) {
    std::cout << "Warning: kinematic state for _lastPos not set!" << std::endl;
    _lastPos.reset(*_body,D);
  }
  D=setKinematic(_pos._info._xM,_t);
  if(D!=_pos._info._xM) {
    std::cout << "Warning: kinematic state for _pos not set!" << std::endl;
    _pos.reset(*_body,D);
  }
  GradInfo newPos(*_body,setKinematic(_pos._info._xM,_t+_dt)),newPos2;
  detectLastContact();
  //ASSERT_MSG(detectLastContact(),"Penetrating initial guess!")
  //normal solve
  T nu=2,alphaMax=1e20,alphaMin=1e-6f;
  T e=energy(newPos,&DE),e2=0,rho=0;
  mask(NULL,&DE,&(newPos._HTheta));
  for(int iter=0; iter<_maxIt;) {
    //update configuration
    update(newPos,newPos2,D,DE,newPos._HTheta,_alpha);
    e2=energy(newPos2,NULL);
    //iteration update
    rho=(e2-e)/D.dot(DE+newPos._HTheta*D/2);
    if(isfinite(e2) && (e2<e && rho>0)) {
      _alpha=_alpha*std::max<T>(1./3.,1-(2*rho-1)*(2*rho-1)*(2*rho-1));
      _alpha=std::max<T>(std::min<T>(_alpha,alphaMax),alphaMin);
      nu=2;
      newPos=newPos2;
      e=energy(newPos,&DE);
      mask(NULL,&DE,&(newPos._HTheta));
      iter++;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " gNorm=" << DE.cwiseAbs().maxCoeff() << " alpha=" << _alpha << " nu=" << nu << " rho=" << rho << std::endl;
      //termination: gradient tolerance
      if(DE.cwiseAbs().maxCoeff()<_gTol)
        break;
    } else {
      _alpha=_alpha*nu;
      nu=nu*2;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " alpha=" << _alpha << " nu=" << nu << std::endl;
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
void ConvHullPBDSimulator::backward() {
  backward(_pos);
}
void ConvHullPBDSimulator::backward(GradInfo& grad,bool debug) {
  Vec DE;
  energy(grad,&DE); //we do not rely on assumption that energy(*,*) has been called
  detectContact(grad._info._TM);
  T coef=1.0/(_dt*_dt),rho;
  Mat3XT GB,MRR,MRt,MtR,Mtt;
  Mat3XT MRRL,MRtL,MtRL,MttL;
  std::shared_ptr<MeshExact> local;
  int nrJ=_body->nrJ();
  int nrDB=_body->nrDOF();
  _MRR.setZero(3,3*nrJ);
  _MRt.setZero(3,3*nrJ);
  _MtR.setZero(3,3*nrJ);
  _Mtt.setZero(3,3*nrJ);
  _MRRL.setZero(3,3*nrJ);
  _MRtL.setZero(3,3*nrJ);
  _MtRL.setZero(3,3*nrJ);
  _MttL.setZero(3,3*nrJ);
  grad._HThetaD.setZero(nrDB,grad._nrVex*3);
  grad._HThetaPTarget.setZero(nrDB,nrDB);
  grad._HThetaDTarget.setZero(nrDB,nrDB);
  grad._HThetaL.setZero(nrDB,nrDB);
  grad._HThetaLL.setZero(nrDB,nrDB);
  grad._HPos.setZero(nrDB,3);
  OMP_PARALLEL_FOR_
  for(int k=0; k<nrJ; k++) {
    //dynamic
    const Joint& J=_body->joint(k);
    int nrD=J.nrDOF();
    //update MC,MCCT
    Mat3T PPT=Mat3T::Zero();
    Vec3T P=Vec3T::Zero();
    if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
      rho=J._M/(1.0*local->vss().size());
      for(int i=0; i<(int)local->vss().size(); i++) {
        P+=local->vss()[i].template cast<T>()*rho;
        PPT+=(local->vss()[i]*(local->vss()[i]).transpose()).template cast<T>()*rho;
      }
    }
    _MRRL.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(grad._info._TM,k)*PPT*ROTI(_lastPos._info._TM,k).transpose())*coef;
    _MRtL.template block<3,3>(0,k*3)+=cross<T>(ROTI(grad._info._TM,k)*P)*coef;
    _MtRL.template block<3,3>(0,k*3)-=cross<T>(ROTI(_lastPos._info._TM,k)*P)*coef;
    _MttL.template block<3,3>(0,k*3)+=Mat3T::Identity()*J._M*coef;
    _MRR.template block<3,3>(0,k*3)+=2*invDoubleCrossMatTrace<T>(ROTI(grad._info._TM,k)*PPT*ROTI(_pos._info._TM,k).transpose())*coef;
    _MRt.template block<3,3>(0,k*3)-=2*cross<T>(ROTI(grad._info._TM,k)*P)*coef;
    _MtR.template block<3,3>(0,k*3)+=2*cross<T>(ROTI(_pos._info._TM,k)*P)*coef;
    _Mtt.template block<3,3>(0,k*3)-=2*Mat3T::Identity()*J._M*coef;
    Mat3X4T A=TRANSI(grad._info._TM,k)-2*TRANSI(_pos._info._TM,k)+TRANSI(_lastPos._info._TM,k);
    Mat3T Mwt,Mtt,C;
    MatX3T HThetaX;
    if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
      T rho=J._M/(1.0*local->vss().size());
      for(int c=0; c<(int)local->vss().size(); c++) {
        HThetaX.setZero(nrDB,3);
        C=cross<T>(ROTI(grad._info._TM,k)*local->vss()[c].template cast<T>());
        Mwt=C*coef*rho;
        Mtt=Mat3T::Identity()*coef*rho;
        Mwt=Mwt*ROT(A)-cross<T>((CTR(A)*rho+ROT(A)*rho*local->vss()[c].template cast<T>())*coef)*ROTI(grad._info._TM,k);
        Mtt=Mtt*ROT(A);
        grad._info.JRCSparse(*_body,k,[&](int row,const Vec3T& JR) {
          parallelAdd<T,3>(HThetaX,row,0,JR.transpose()*Mwt);
        },[&](int row,const Vec3T& JC) {
          parallelAdd<T,3>(HThetaX,row,0,JC.transpose()*Mtt);
        });
        for(int i=0; i<_body->nrDOF(); i++)
          for(int j=0; j<3; j++)
            parallelAdd(grad._HThetaD(i,j+(grad._polytopes[k].getVertexId()[0]+c)*3),HThetaX(i,j));
      }
    }
    //P controller
    Eigen::Matrix<T,-1,1,0,3,1> coef;
    if(J._control.size()>=nrD)
      coef=J._control.segment(0,nrD).template cast<T>();
    else coef.setOnes(nrD);
    if(_params[k]._kp>0) {
      grad._HThetaPTarget.diagonal().segment(J._offDOF,nrD).array()-=coef.array()*_params[k]._kp;
    }
    //D controller
    if(_params[k]._kd>0) {
      grad._HThetaDTarget.diagonal().segment(J._offDOF,nrD).array()-=coef.array()*_params[k]._kd/_dt;
      grad._HThetaL.diagonal().segment(J._offDOF,nrD).array()-=coef.array()*_params[k]._kd/_dt/_dt;
    }
    //joint limit
  }
  //pos grad
  //getJointPosGrad(grad,1);
  //contact
  normalEnergy(grad,NULL,true);
  tangentEnergy(grad,NULL,true);
  grad._info.toolA(*_body,_pos._info,mapM(MRR=_MRR),mapM(MRt=_MRt),mapM(MtR=_MtR),mapM(Mtt=_Mtt),[&](int r,int c,T val) {
    grad._HThetaL(r,c)+=val;
  });
  grad._info.toolA(*_body,_lastPos._info,mapM(MRRL=_MRRL),mapM(MRtL=_MRtL),mapM(MtRL=_MtRL),mapM(MttL=_MttL),[&](int r,int c,T val) {
    grad._HThetaLL(r,c)+=val;
  });
  //custom
  if(_custom)
    _custom->backward(grad,_pos,_lastPos,_dt);
  if(!debug) {
    MatT diag=MatT::Zero(_body->nrDOF(),_body->nrDOF());
    MatT HTheta=grad._HTheta;
    MatT HThetaPTarget=grad._HThetaPTarget;
    MatT HThetaDTarget=grad._HThetaDTarget;
    MatT HThetaD=grad._HThetaD;
    MatT HThetaL=grad._HThetaL;
    MatT HThetaLL=grad._HThetaLL;
    MatT HThetaDesign=grad._HTheta+grad._HThetaL+grad._HThetaLL;
    //mask out design parameters
    mask(&diag,NULL,&HTheta);
    mask(NULL,NULL,&HThetaPTarget);
    mask(NULL,NULL,&HThetaDTarget);
    mask(NULL,NULL,NULL,&HThetaD);
    mask(NULL,NULL,&HThetaL);
    mask(NULL,NULL,&HThetaLL);
    for(int i=0; i<diag.rows(); i++)
      if(!isfinite(diag(i,i)))
        diag(i,i)=1;
    //solve for non-design derivatives
    _sol->compute(MatT(HTheta+diag));
    grad._HThetaPTarget=-_sol->solve(HThetaPTarget);
    grad._HThetaDTarget=-_sol->solve(HThetaDTarget);
    grad._HThetaD=-_sol->solve(HThetaD);
    grad._HThetaL=-_sol->solve(HThetaL);
    grad._HThetaLL=-_sol->solve(HThetaLL);
    //compute design derivatives
    maskNonDesign(HThetaDesign);
    grad._HThetaDesign=-_sol->solve(HThetaDesign);
    markDesign(grad._HThetaDesign);
  } else {
    grad._HThetaDesign=grad._HTheta+grad._HThetaL+grad._HThetaLL;
    maskNonDesign(grad._HThetaDesign);
  }
}
void ConvHullPBDSimulator::stepwithbackward(bool debug) {
  if(_body->nrDOF()==0)
    return;
  //std::cout<<_count<<std::endl;
  //if(_count==500) exit(0);
  //auto start = std::chrono::high_resolution_clock::now();
  //update kinematic state
  GradInfo newPos(*_body,setKinematic(_pos._info._xM,_t+_dt)),newPos2;
  detectLastContact();
  //ASSERT_MSG(detectLastContact(),"Penetrating initial guess!")
  //normal solve
  Vec D,DE;
  D=setKinematic(_lastPos._info._xM,_t-_dt);
  if(D!=_lastPos._info._xM) {
    std::cout << "Warning: kinematic state for _lastPos not set!" << std::endl;
    _lastPos.reset(*_body,D);
  }
  D=setKinematic(_pos._info._xM,_t);
  if(D!=_pos._info._xM) {
    std::cout << "Warning: kinematic state for _pos not set!" << std::endl;
    _pos.reset(*_body,D);
  }
  T nu=2,alphaMax=1e20,alphaMin=1e-6f;
  T e=energy(newPos,&DE),e2=0,rho=0;
  mask(NULL,&DE,&(newPos._HTheta));
  for(int iter=0; iter<_maxIt;) {
    //update configuration
    update(newPos,newPos2,D,DE,newPos._HTheta,_alpha);
    e2=energy(newPos2,NULL);
    //iteration update
    rho=(e2-e)/D.dot(DE+newPos._HTheta*D/2);
    if(isfinite(e2) && (e2<e && rho>0)) {
      _alpha=_alpha*std::max<T>(1./3.,1-(2*rho-1)*(2*rho-1)*(2*rho-1));
      _alpha=std::max<T>(std::min<T>(_alpha,alphaMax),alphaMin);
      nu=2;
      newPos=newPos2;
      e=energy(newPos,&DE);
      mask(NULL,&DE,&(newPos._HTheta));
      iter++;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " gNorm=" << DE.cwiseAbs().maxCoeff() << " alpha=" << _alpha << " nu=" << nu << " rho=" << rho << std::endl;
      //termination: gradient tolerance
      if(DE.cwiseAbs().maxCoeff()<_gTol)
        break;
    } else {
      _alpha=_alpha*nu;
      nu=nu*2;
      if(_output)
        std::cout << "Iter=" << iter << " E=" << e << " alpha=" << _alpha << " nu=" << nu << std::endl;
      //termination: numerical issue
      if(_alpha>=alphaMax) {
        std::cout << "Failed!" << std::endl;
        _alpha=alphaMin;
        break;
      }
    }
  }
  backward(newPos,debug);
  //update
  _lastPos=_pos;
  _pos=newPos;
  _t+=_dt;
}
void ConvHullPBDSimulator::detectCurrentContact() {
  FUNCTION_NOT_IMPLEMENTED
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::gTol() const {
  return _gTol;
}
void ConvHullPBDSimulator::setGTol(T gTol) {
  _gTol=gTol;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::x0() const {
  return _barrier._x0;
}
void ConvHullPBDSimulator::setX0(T x0) {
  _barrier._x0=(double)x0;
}
bool ConvHullPBDSimulator::hardLimit() const {
  return _hardLimit;
}
void ConvHullPBDSimulator::setHardLimit(bool hardLimit) {
  _hardLimit=hardLimit;
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::pos() const {
  return _pos._info._xM;
}
void ConvHullPBDSimulator::setPos(const Vec& pos) {
  _pos.reset(*_body,setKinematic(pos,_t));
}
ConvHullPBDSimulator::Vec ConvHullPBDSimulator::vel() const {
  return (_pos._info._xM-_lastPos._info._xM)/_dt;
}
void ConvHullPBDSimulator::setVel(const Vec& vel) {
  _lastPos.reset(*_body,setKinematic(_pos._info._xM-vel*_dt,_t-_dt));
}
void ConvHullPBDSimulator::setCoefBarrier(T coefBarrier) {
  _coefBarrier=coefBarrier;
}
void ConvHullPBDSimulator::debugEnergy(T scale,const T* customDelta) {
  DEFINE_NUMERIC_DELTA_T(T)
  if(customDelta)
    DELTA=*customDelta;
  while(true) {
    //generate random pose
    setDesign(Vec::Random(_body->nrDOF())*scale);
    GradInfo newPos,newPos2,pos,lastPos;
    _pos.reset(*_body,setKinematic(Vec::Random(_body->nrDOF())*scale,_t));
    _lastPos.reset(*_body,setKinematic(Vec::Random(_body->nrDOF())*scale,_t));
    newPos.reset(*_body,setKinematic(Vec::Random(_body->nrDOF())*scale,_t));
    pos=_pos;
    lastPos=_lastPos;

    if(!detectLastContact())
      continue;
    pos.reset(*_body,_pos._info._xM);
    int nrJ=_body->nrJ();

    //generate random joint
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
    }
    Vec DE,DE2,dx=Vec::Random(_body->nrDOF()),dc,X,design=_design;
    setPD(P,D);

    //debug DE/DDE
    T e=energy(newPos,&DE);
    mask(NULL,&DE,NULL);
    backward(newPos,true);
    newPos2.reset(*_body,setKinematic(newPos._info._xM+dx*DELTA,_t));
    T e2=energy(newPos2,&DE2);
    mask(NULL,&DE2,&(newPos._HTheta));
    DEBUG_GRADIENT("DE",DE.dot(dx),DE.dot(dx)-(e2-e)/DELTA)
    DEBUG_GRADIENT("DDE",(newPos._HTheta*dx).norm(),(newPos._HTheta*dx-(DE2-DE)/DELTA).norm())

    //debug DDE-L
    _pos.reset(*_body,setKinematic(pos._info._xM+dx*DELTA,_t));
    _lastPos=lastPos;
    newPos2=newPos;
    if(!detectLastContact())
      continue;
    energy(newPos2,&DE2);
    mask(NULL,&DE2,&(newPos._HThetaL));
    DEBUG_GRADIENT("DDE-L",(newPos._HThetaL*dx).norm(),(newPos._HThetaL*dx-(DE2-DE)/DELTA).norm())

    //debug DDE-L
    _pos=pos;
    _lastPos.reset(*_body,setKinematic(lastPos._info._xM+dx*DELTA,_t));
    newPos2=newPos;
    if(!detectLastContact())
      continue;
    energy(newPos2,&DE2);
    mask(NULL,&DE2,&(newPos._HThetaLL));
    DEBUG_GRADIENT("DDE-LL",(newPos._HThetaLL*dx).norm(),(newPos._HThetaLL*dx-(DE2-DE)/DELTA).norm())

    //debug DDE-P
    _pos=pos;
    _lastPos=lastPos;
    setPD(P+dx*DELTA,D);
    if(!detectLastContact())
      continue;
    energy(newPos2,&DE2);
    mask(NULL,&DE2,&(newPos._HThetaPTarget));
    DEBUG_GRADIENT("DDE-P",(newPos._HThetaPTarget*dx).norm(),(newPos._HThetaPTarget*dx-(DE2-DE)/DELTA).norm())

    //debug DDE-D
    _pos=pos;
    _lastPos=lastPos;
    setPD(P,D+dx*DELTA);
    if(!detectLastContact())
      continue;
    energy(newPos2,&DE2);
    mask(NULL,&DE2,&(newPos._HThetaDTarget));
    DEBUG_GRADIENT("DDE-D",(newPos._HThetaDTarget*dx).norm(),(newPos._HThetaDTarget*dx-(DE2-DE)/DELTA).norm())
    setPD(P,D);

    //debug DDE-Design
    setDesign(design+dx*DELTA);
    _pos.reset(*_body,setKinematic(pos._info._xM,_dt));
    _lastPos.reset(*_body,setKinematic(lastPos._info._xM,_dt));
    newPos2.reset(*_body,setKinematic(newPos._info._xM,_dt));
    if(!detectLastContact())
      continue;
    energy(newPos2,&DE2);
    mask(NULL,&DE2,NULL);
    DEBUG_GRADIENT("DDE-Design",(newPos._HThetaDesign*dx).norm(),(newPos._HThetaDesign*dx-(DE2-DE)/DELTA).norm())
    setDesign(design);

    //debug DDE-XL
    _pos=pos;
    _lastPos=lastPos;
    if(!detectLastContact())
      continue;
    //move convex hull
    dc=Vec::Random(newPos._nrVex*3);
    updateConvexPoints(dc*DELTA);
    newPos2.reset(*_body,newPos._info._xM); //after updating the vertices, we have to reset
    energy(newPos2,&DE2);
    mask(NULL,&DE2,NULL,&(newPos._HThetaD));
    DEBUG_GRADIENT("DDE-XL",(newPos._HThetaD*dc).norm(),(newPos._HThetaD*dc-(DE2-DE)/DELTA).norm())
    debugBVHEnergy(newPos);
    break;
  }
}
void ConvHullPBDSimulator::debugBackward(T scale,const T* customDelta) {
  DEFINE_NUMERIC_DELTA_T(T)
  if(customDelta)
    DELTA=*customDelta;
  while(true) {
    //generate random pose
    GradInfo nextPos,pos,lastPos;
    setDesign(Vec::Random(_body->nrDOF())*scale);
    _pos.reset(*_body,setKinematic(Vec::Random(_body->nrDOF())*scale,_t));
    _lastPos.reset(*_body,setKinematic(Vec::Random(_body->nrDOF())*scale,_t-_dt));
    if(!detectLastContact())
      continue;
    pos=_pos;
    lastPos=_lastPos;
    int nrJ=_body->nrJ();

    //generate random joint
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
    }
    T e=energy(_pos,NULL);
    if(!isfinite(e))
      continue;
    setPD(P,D);

    step();
    nextPos=_pos;
    _pos=pos;
    _lastPos=lastPos;
    backward(nextPos);
    Vec dx=Vec::Random(_body->nrDOF()),dc,X,design=_design;

    //debug DTDL
    _pos.reset(*_body,setKinematic(pos._info._xM+dx*DELTA,_t));
    _lastPos=lastPos;
    _alpha=1;
    step();
    DEBUG_GRADIENT("DTDL",(nextPos._HThetaL*dx).norm(),(nextPos._HThetaL*dx-(_pos._info._xM-nextPos._info._xM)/DELTA).norm())

    //debug DTDLL
    _pos=pos;
    _lastPos.reset(*_body,setKinematic(lastPos._info._xM+dx*DELTA,_t));
    _alpha=1;
    step();
    DEBUG_GRADIENT("DTDLL",(nextPos._HThetaLL*dx).norm(),(nextPos._HThetaLL*dx-(_pos._info._xM-nextPos._info._xM)/DELTA).norm())

    //debug DTDP
    _pos=pos;
    _lastPos=lastPos;
    setPD(P+dx*DELTA,D);
    _alpha=1;
    step();
    DEBUG_GRADIENT("DTDP",(nextPos._HThetaPTarget*dx).norm(),(nextPos._HThetaPTarget*dx-(_pos._info._xM-nextPos._info._xM)/DELTA).norm())

    //debug DTDD
    _pos=pos;
    _lastPos=lastPos;
    setPD(P,D+dx*DELTA);
    _alpha=1;
    step();
    DEBUG_GRADIENT("DTDD",(nextPos._HThetaDTarget*dx).norm(),(nextPos._HThetaDTarget*dx-(_pos._info._xM-nextPos._info._xM)/DELTA).norm()) //recover D
    setPD(P,D);

    //debug Design
    setDesign(design+dx*DELTA);
    _pos.reset(*_body,setKinematic(pos._info._xM,_dt));
    _lastPos.reset(*_body,setKinematic(lastPos._info._xM,_dt));
    step();
    DEBUG_GRADIENT("DTDDesign",(nextPos._HThetaDesign*dx).norm(),(nextPos._HThetaDesign*dx-(_pos._info._xM-nextPos._info._xM)/DELTA).norm())
    setDesign(design);

    //debug DTDXLx
    _pos=pos;
    _lastPos=lastPos;
    dc=Vec::Random(_pos._nrVex*3);
    updateConvexPoints(dc*DELTA);
    _alpha=1;
    step();
    DEBUG_GRADIENT("DTDXL",(nextPos._HThetaD*dc).norm(),(nextPos._HThetaD*dc-(_pos._info._xM-nextPos._info._xM)/DELTA).norm())
    break;
  }
}
void ConvHullPBDSimulator::debugBVHEnergy(GradInfo& grad) {
  DEFINE_NUMERIC_DELTA_T(T)
  T E1=0,E2=0,val;
  //BVH-assisted
  detectContact(grad._info._TM);
  for(auto& m:_manifolds) {
    GJKPolytope<T>& mA=m._jidA<0?_obs[m._sidA]:grad._polytopes[m._jidA];
    GJKPolytope<T>& mB=m._jidB<0?_obs[m._sidB]:grad._polytopes[m._jidB];
    if(!mA.mesh() || !mB.mesh())
      continue;
    CCBarrierConvexEnergy<T,Barrier> cc(mA,mB,_barrier,0,&grad,_coefBarrier);
    ASSERT_MSG(cc.eval(&val,_body.get(),&grad,NULL,NULL,NULL),"Invalid debug configuration!")
    E1+=val;
  }
  //brute force
  int nrJ=_body->nrJ();
  int nrS=(int)_obs.size();
  for(int k1=0; k1<nrJ; k1++) {
    GJKPolytope<T>& mA=grad._polytopes[k1];
    if(!mA.mesh())
      continue;
    for(int k2=k1+1; k2<nrJ; k2++) {
      if(_contact->getExclude().find(Eigen::Matrix<int,2,1>(k1,k2))!=_contact->getExclude().end())
        continue;
      GJKPolytope<T>& mB=grad._polytopes[k2];
      if(!mB.mesh())
        continue;
      //class check
      if(_body->joint(_body->joint(k1)._parent)._class!=_body->joint(k2)._class && _body->joint(_body->joint(k2)._parent)._class!=_body->joint(k1)._class)
        if(_body->joint(k1)._class!=_body->joint(k2)._class && _body->joint(k2)._class!=_body->joint(k1)._class) {
          CCBarrierConvexEnergy<T,Barrier> cc(mA,mB,_barrier,0,&grad,_coefBarrier);
          ASSERT_MSGV(cc.eval(&val,_body.get(),&grad,NULL,NULL,NULL),"Invalid dynamic-dynamic configuration between (%d,%d)!",k1,k2)
          E2+=val;
        }
    }
    for(int k2=0; k2<nrS; k2++) {
      CCBarrierConvexEnergy<T,Barrier> cc(mA,_obs[k2],_barrier,0,&grad,_coefBarrier);
      ASSERT_MSGV(cc.eval(&val,_body.get(),&grad,NULL,NULL,NULL),"Invalid dynamic-static configuration between (%d,%d)!",k1,k2)
      E2+=val;
    }
  }
  DEBUG_GRADIENT("BVH-E",E1,E1-E2)
}
void ConvHullPBDSimulator::setIsDesign(const std::vector<char>& isDesign,int startJID) {
  for(int k=startJID; k<_body->nrJ(); k++) {
    //param
    auto& param=getJointPhysicsParameter(k);
    //set design
    Joint& J=_body->joint(k);
    int nrDJ=J.nrDOF();
    for(int c=0; c<nrDJ; c++)
      param._isDesign[c]=isDesign[J._offDOF+c];
  }
}
void ConvHullPBDSimulator::setPD(const Vec& P, const Vec& D,int st,int en,T* kp,T* kd) {
  _P=P;
  _D=D;
  for(int k=st; k<_body->nrJ()-en; k++) {
    //control
    Joint& J=_body->joint(k);
    J._control.setOnes(J.nrDOF());
    //param
    auto& param=getJointPhysicsParameter(k);
    if(kp) param._kp=*kp;
    if(kd) param._kd=*kd;
    param._tarP=[&](T,int)->Vec {
      return _P.segment(J._offDOF,J.nrDOF());
    };
    param._tarD=[&](T,int)->Vec {
      return _D.segment(J._offDOF,J.nrDOF());
    };
  }
}
void ConvHullPBDSimulator::updateConvexPoints(const Vec& X) {
  int nrJ=_body->nrJ();
  for(int k=0; k<nrJ; k++)
    if(_pos._polytopes[k].jid()>=0) {
      Eigen::Matrix<int,2,1> off=_pos._polytopes[k].getVertexId();
      Vec dc=X.segment(off[0]*3,(off[1]-off[0])*3);
      moveMesh(_body->joint(k),dc);
    }
}
//helper
void ConvHullPBDSimulator::detectContact(const Mat3XT& t) {
  _manifolds.clear();
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  GEOMETRY_SCALAR x0=10000;//2*(_barrier._x0+_d0)/(1-_barrier._x0);
  _contact->generateManifolds(x0,true,_manifolds,t.template cast<GEOMETRY_SCALAR>());
}
bool ConvHullPBDSimulator::detectLastContact() {
  detectContact(_pos._info._TM);
  Vec DE;
  T e=normalEnergy(_pos,&DE,false);
  if(!isfinite(e))
    return false;
  _manifoldsLast.clear();
  for(auto &m:_manifolds)
    _manifoldsLast.push_back(m);
  return !_manifoldsLast.empty();
}
void ConvHullPBDSimulator::update(const GradInfo& newPos,GradInfo& newPos2,Vec& D,const Vec& DE,const MatT& DDE,T alpha) const {
  MatT DDER=DDE;
  DDER.diagonal().array()+=_alpha;
  //compute
  _sol->compute(DDER);
  //update
  D=-_sol->solve(MatT(DE));
  newPos2.reset(*_body,newPos._info._xM+D);
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::energy(GradInfo& grad,Vec* DE) {
  T E=0;
  Vec3T P,PTotal;
  T MTotal=0;
  Mat3T PPT;
  Mat3XT GB,MRR,MRt,MtR,Mtt;
  detectContact(grad._info._TM);
  int nrJ=_body->nrJ();
  int nrD=_body->nrDOF();
  if(DE) {
    DE->setZero(nrD);
    grad._DTG.setZero(3,4*nrJ);
    _MRR.setZero(3,3*nrJ);
    _MRt.setZero(3,3*nrJ);
    _MtR.setZero(3,3*nrJ);
    _Mtt.setZero(3,3*nrJ);
    grad._HTheta.setZero(nrD,nrD);
  }
  grad._centre.resize(nrJ);
  //contact
  E+=normalEnergy(grad,DE);
  E+=tangentEnergy(grad,DE);
  for(int k=0; k<nrJ; k++) {
    const Joint& J=_body->joint(k);
    nrD=J.nrDOF();
    //update MC,MCCT
    PPT.setZero();
    P.setZero();
    PTotal.setZero();
    if(std::dynamic_pointer_cast<MeshExact>(J._mesh)) {
      std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
      T rho=J._M/(1.0*local->vss().size());
      for(int i=0; i<(int)local->vss().size(); i++) {
        P+=local->vss()[i].template cast<T>()*rho;
        PPT+=(local->vss()[i]*(local->vss()[i]).transpose()).template cast<T>()*rho;
        PTotal+=local->vss()[i].template cast<T>()/(1.0*local->vss().size());
      }
      grad._centre[k]=(ROT(TRANSI(grad._info._TM,k))*PTotal+CTR(TRANSI(grad._info._TM,k)));
      MTotal+=J._M;
    }
    //inertial
    E+=energyInertial(grad._info,_pos._info,_lastPos._info,
                      k,0,J._M,P,PPT,
                      DE?&(grad._DTG):NULL,_MRR,_MRt,_MtR,_Mtt);
    //PD controller
    E+=energyPDController(mapV2CV(grad._info._xM),mapV2CV(_pos._info._xM),k,J,nrD,DE,&(grad._HTheta));
    //joint limit
    if(_hardLimit) {
      if(!JointLimit::energy(mapV2CV(grad._info._xM),J,nrD,E,DE,&(grad._HTheta),_barrier))
        return std::numeric_limits<double>::infinity();
    } else E+=JointLimit::energy(mapV2CV(grad._info._xM),J,nrD,DE,&(grad._HTheta));
  }
  if(!isfinite(E))
    return E;
  //joints
  for(const auto& joint:_joints)
    E+=joint.energy(*_body,grad._info,DE,grad._HTheta,grad._DTG,_MRR,_MRt,_MtR,_Mtt,true);
  if(DE) {
    //gradient
    grad._info.DTG(*_body,mapM(GB=grad._DTG),mapV(*DE));
    //hessian
    grad._info.toolAB(*_body,mapM(MRR=_MRR),mapM(MRt=_MRt),mapM(MtR=_MtR),mapM(Mtt=_Mtt),mapM(GB=grad._DTG),[&](int r,int c,T val) {
      grad._HTheta(r,c)+=val;
    });
  }
  //custom PBD energy
  if(_custom)
    E+=_custom->energy(grad,_pos,_lastPos,_dt,DE);
  return E;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::normalEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifolds.size(); id++) {
    if(!isfinite(E))
      continue;
    ContactManifold& m=_manifolds[id];
    GJKPolytope<T>& mA=m._jidA<0?_obs[m._sidA]:grad._polytopes[m._jidA];
    GJKPolytope<T>& mB=m._jidB<0?_obs[m._sidB]:grad._polytopes[m._jidB];
    if(!mA.mesh() || !mB.mesh())
      continue;
    //compute energy/gradient/hessian
    T val=0;
    CCBarrierConvexEnergy<T,Barrier> cc(mA,mB,_barrier,_d0,&grad,_coefBarrier);
    if(!backward) {
      if(!cc.eval(&val,_body.get(),DE?&grad:NULL,&m._DNDX,NULL,NULL))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else {
        parallelAdd<T>(E,val);
        m._x=cc.getX();
      }
    } else {
      std::vector<MatX3T> HThetaD1,HThetaD2;
      cc.evalBackward(_body.get(),&grad,&HThetaD1,&HThetaD2);
    }
  }
  return E;
}
ConvHullPBDSimulator::T ConvHullPBDSimulator::tangentEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifoldsLast.size(); id++) {
    if(!isfinite(E))
      continue;
    ContactManifold& m=_manifoldsLast[id];
    GJKPolytope<T>& mA=m._jidA<0?_obs[m._sidA]:grad._polytopes[m._jidA];
    GJKPolytope<T>& mB=m._jidB<0?_obs[m._sidB]:grad._polytopes[m._jidB];
    GJKPolytope<T>& mALast=m._jidA<0?_obs[m._sidA]:_pos._polytopes[m._jidA];
    GJKPolytope<T>& mBLast=m._jidB<0?_obs[m._sidB]:_pos._polytopes[m._jidB];
    Vec4T x=m._x;
    if(!mA.mesh() || !mB.mesh())
      continue;
    //compute energy/gradient/hessian
    T val=0;
    CCBarrierFrictionEnergy<T,Barrier> cf(mA,mB,mALast,mBLast,x,_barrier,_d0,&grad,_coefBarrier,_dt);
    cf.fri()=std::max<T>(m._jidA>=0?_params[m._jidA]._friction:0,m._jidB>=0?_params[m._jidB]._friction:0);
    if(!backward) {
      if(!cf.eval(&val,_body.get(),DE?&grad:NULL,NULL,NULL))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else parallelAdd<T>(E,val);
    } else {
      std::vector<MatX3T> HThetaD1,HThetaD2;
      cf.evalBackward(_body.get(),&grad,&_pos,&m._DNDX,&HThetaD1,&HThetaD2);
    }
  }
  return E;
}
}
