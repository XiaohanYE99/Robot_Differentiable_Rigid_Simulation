#include <iostream>
#include "ThetaTrajectory.h"
#include "TrajectorySIPEnergy.h"
#include <Utils/SparseUtils.h>
#include <Utils/Interp.h>
#include <random>

namespace PHYSICSMOTION {
template <typename T>
bool ThetaTrajectory<T>::Segment::read(std::istream& is,IOData*) {
  readBinaryData(_variableToCP,is);
  readBinaryData(_variableToCPDim,is);
  readBinaryData(_index,is);
  return is.good();
}
template <typename T>
bool ThetaTrajectory<T>::Segment::write(std::ostream& os,IOData*) const {
  writeBinaryData(_variableToCP,os);
  writeBinaryData(_variableToCPDim,os);
  writeBinaryData(_index,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> ThetaTrajectory<T>::Segment::copy() const {
  return std::shared_ptr<SerializableBase>(new ThetaTrajectory::Segment);
}
template <typename T>
std::string ThetaTrajectory<T>::Segment::type() const {
  return typeid(ThetaTrajectory::Segment).name();
}
template <typename T>
ThetaTrajectory<T>::ThetaTrajectory(int dimension,int order,int numSegment) {
  if(numSegment==0) {
    _numSegment=0;
    _numNode=1;
    _order=0;
    _dimension=dimension;
    return;
  }
  if(order<=4)
    throw std::invalid_argument("Order must be greater than 4!");
  _order=order;
  _dimension=dimension;
  _numSegment=numSegment;
  _numNode=0;
  //The first/last segment introduces (order-2) independent nodes/control points
  //Every intermediary segment introduces (order-5) independent nodes/control points
  //Every joint has 3 nodes,providing 3 control points for2 curves
  std::vector<STrips> tripss;
  for(int i=0; i<numSegment; i++) {
    STrips trips;
    int k=0;
    if(i==0) {
      for(int j=0; j<order-2; j++)
        addNode(trips,k);
      addJointLeft(trips,k);
    } else if(i==numSegment-1) {
      addJointRight(trips,k);
      for(int j=0; j<order-2; j++)
        addNode(trips,k);
    } else {
      addJointRight(trips,k);
      for(int j=0; j<order-5; j++)
        addNode(trips,k);
      addJointLeft(trips,k);
    }
    tripss.push_back(trips);
  }
  //segment
  for(int i=0; i<numSegment; i++) {
    Segment seg;
    seg._variableToCP.resize(order+1,_numNode);
    seg._variableToCP.setFromTriplets(tripss[i].begin(),tripss[i].end());
    seg._variableToCPDim=kronecker(seg._variableToCP,dimension);
    seg._index=i;
    _segments.push_back(seg);
  }
  //subdivide stencil
  _curve=BezierCurve<T>(order,dimension);
}
template <typename T>
bool ThetaTrajectory<T>::read(std::istream& is,IOData*) {
  readBinaryData(_curve,is);
  readBinaryData(_segments,is);
  readBinaryData(_numSegment,is);
  readBinaryData(_numNode,is);
  readBinaryData(_order,is);
  readBinaryData(_dimension,is);
  return is.good();
}
template <typename T>
bool ThetaTrajectory<T>::write(std::ostream& os,IOData*) const {
  writeBinaryData(_curve,os);
  writeBinaryData(_segments,os);
  writeBinaryData(_numSegment,os);
  writeBinaryData(_numNode,os);
  writeBinaryData(_order,os);
  writeBinaryData(_dimension,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> ThetaTrajectory<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new ThetaTrajectory<T>);
}
template <typename T>
std::string ThetaTrajectory<T>::type() const {
  return typeid(ThetaTrajectory).name();
}
template <typename T>
const std::vector<typename ThetaTrajectory<T>::Segment>& ThetaTrajectory<T>::getSegments() const {
  return _segments;
}
template <typename T>
std::vector<typename ThetaTrajectory<T>::Segment>& ThetaTrajectory<T>::getSegments() {
  return _segments;
}
template <typename T>
const BezierCurve<T>& ThetaTrajectory<T>::getCurve() const {
  return _curve;
}
template <typename T>
int ThetaTrajectory<T>::getNumSegment() const {
  return _numSegment;
}
template <typename T>
int ThetaTrajectory<T>::getDimension() const {
  return _dimension;
}
template <typename T>
int ThetaTrajectory<T>::getNumNode() const {
  return _numNode;
}
template <typename T>
int ThetaTrajectory<T>::getNumDOF() const {
  return _numNode*_dimension;
}
template <typename T>
int ThetaTrajectory<T>::getOrder() const {
  return _order;
}
template <typename T>
T ThetaTrajectory<T>::getLength(const Vec& variable,T interval) const {
  T len=0;
  for(T t=0; t+interval<getNumSegment(); t+=interval)
    len+=(getPoint(variable,t)-getPoint(variable,t+interval)).norm();
  return len;
}
template <typename T>
typename ThetaTrajectory<T>::Vec ThetaTrajectory<T>::getPoint(const Vec& variable,T t) const {
  if(_segments.empty()) {
    ASSERT_MSG(t==0,"StaticPose can only accept t=0!")
    return variable;
  }
  for(const auto& seg:_segments) {
    if(t<=1) {
      Vec controlPointsFlatten=seg._variableToCPDim*variable;
      return _curve.getPoint(Eigen::Map<const MatT>(controlPointsFlatten.data(),_dimension,_order+1),t);
    } else t-=1;
  }
  ASSERT_MSG(false,"Time Out of Range in GetPoint")
}
template <typename T>
typename ThetaTrajectory<T>::Vec ThetaTrajectory<T>::getDerivative(const Vec& variable,T t,int d) const {
  if(_segments.empty()) {
    ASSERT_MSG(t==0,"StaticPose can only accept t=0!")
    return Vec::Zero(variable.size());
  }
  for(const auto&seg:_segments) {
    if(t<=1) {
      Vec controlPointsFlatten=seg._variableToCPDim*variable;
      return _curve.getDerivative(Eigen::Map<const MatT>(controlPointsFlatten.data(),_dimension,_order+1),t,d);
    } else t-=1;
  }
  ASSERT_MSG(false,"Time Out of Range in GetDerivative")
}
template <typename T>
void ThetaTrajectory<T>::getControlPointCoeff(STrips& trips,int rowOff,int segment,int dimension,int cpId) const {
  ASSERT(segment>=0 && segment<(int)_segments.size())
  ASSERT(dimension>=0 && dimension<_dimension)
  ASSERT(cpId>=0 && cpId<_order+1)

  SMatT unit;
  unit.resize(_segments[segment]._variableToCPDim.rows(),1);
  unit.coeffRef(cpId*_dimension+dimension,0)=1;
  unit=unit.transpose()*_segments[segment]._variableToCPDim;
  addBlock(trips,rowOff,0,unit);
}
template <typename T>
typename ThetaTrajectory<T>::SMatT ThetaTrajectory<T>::getControlPointCoeffAll(int dimension) const {
  if(_segments.empty()) {
    MatT ret=MatT::Zero(1,getNumDOF());
    ret(0,dimension)=1;
    return ret.sparseView();
  }
  STrips trips;
  int rowOff=0;
  for(int i=0; i<(int)_segments.size(); i++)
    for(int j=0; j<_order+1; j++) {
      if(i>0 && j==0)
        continue;
      getControlPointCoeff(trips,rowOff,i,dimension,j);
      rowOff++;
    }
  SMatT ret;
  ret.resize(rowOff,getNumDOF());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T>
void ThetaTrajectory<T>::debugControlPointAllDifference() const {
  T minDiff=std::numeric_limits<double>::max();
  SMatT CP=getControlPointCoeffAll(0);
  for(int i=0; i<CP.rows(); i++)
    for(int j=0; j<i; j++) {
      T diff=(CP.row(i)-CP.row(j)).toDense().norm();
      minDiff=std::min<T>(diff,minDiff);
    }
  std::cout << "MinDiff in control point rows: " << minDiff << std::endl;
}
template <typename T>
void ThetaTrajectory<T>::assembleEnergy(T t,GFunc* G,HFunc* H,const Vec* GTheta,const MatT* HTheta) const {
  if(_segments.empty()) {
    ASSERT_MSG(t==0,"StaticPose can only accept t=0!")
    if(G)
      (*G)(0,*GTheta);
    if(H)
      (*H)(0,0,*HTheta);
    return;
  }
  for(int segNum=0; segNum<(int)_segments.size(); ++segNum) {
    const auto seg=_segments[segNum];
    if(t<=1) {
      Vec timePolyBasis=_curve.getCoefficientMatrix()*_curve.getBasis(t);
      if(G) {
        ASSERT_MSG(GTheta,"G Must come with GTheta!")
        Vec tmpG=Vec::Zero(_dimension*(_order+1));
        for(int i=0; i<_order+1; ++i) {
          tmpG.block(i*_dimension,0,_dimension,1)=(*GTheta)*timePolyBasis(i);
        }
        (*G)(0,seg._variableToCPDim.transpose()*tmpG.sparseView());
      }
      if(H) {
        ASSERT_MSG(HTheta,"H Must come with HTheta!")
        MatT tmpH=MatT::Zero(_dimension*(_order+1),_dimension*(_order+1));
        for(int i=0; i<_order+1; ++i) {
          for(int j=0; j<_order+1; ++j) {
            tmpH.block(i*_dimension,j*_dimension,_dimension,_dimension)=(*HTheta)*timePolyBasis(i)*timePolyBasis(j);
          }
        }
        (*H)(0,0,seg._variableToCPDim.transpose()*tmpH.sparseView()*seg._variableToCPDim);
      }
      return;
    } else t-=1;
  }
  ASSERT_MSG(false,"Time Out of Range in AssembleEnergy")
}
template <typename T>
void ThetaTrajectory<T>::assembleEnergy(T t,Vec* G,MatT* H,const Vec* GTheta,const MatT* HTheta) const {
  GFunc GF=[&](int off,const Vec& val)->void{
    SIPEnergy<T>::template parallelAdd<-1>(*G,off,val);
  };
  HFunc HF=[&](int offr,int offc,const MatT& val)->void{
    SIPEnergy<T>::template parallelAdd<-1,-1>(*H,offr,offc,val);
  };
  assembleEnergy(t,G?&GF:NULL,H?&HF:NULL,GTheta,HTheta);
}
template <typename T>
void ThetaTrajectory<T>::smoothnessRegularizer(const Vec& variable,T* E,Vec* G,MatT* H,int d) const {
  if(_segments.empty())
    return;
  SMatT hessian;
  hessian.resize(variable.size(),variable.size());
  for(const auto& seg:_segments) {
    Vec controlPointsFlatten=seg._variableToCPDim*variable;
    hessian+=seg._variableToCPDim.transpose()*_curve.getBasisCrossIntegralFlatten(d)*seg._variableToCPDim;
  }
  if(E)
    (*E)+=variable.dot(hessian*variable)*0.5;
  if(G)
    (*G)+=hessian*variable;
  if(H)
    (*H)+=hessian;
}
template <typename T>
void ThetaTrajectory<T>::smoothnessBF(const Vec& variable,T* E,int segments,int d) const {
  T res=0;
  for(const auto& seg:_segments) {
    Vec controlPointsFlatten=seg._variableToCPDim*variable;
    Eigen::Map<const MatT> controlPoints=Eigen::Map<const MatT>(controlPointsFlatten.data(),_dimension,_order+1);
    for(int i=0; i<segments; ++i) {
      Vec theta=_curve.getDerivative(controlPoints,((T)i*2+1)/2./segments,d);
      res+=pow(theta.norm(),2)/segments;
    }
  }
  (*E)+=res*0.5;
}
template <typename T>
void ThetaTrajectory<T>::debugSmoothness(int segments,int d) const {
  DEFINE_NUMERIC_DELTA_T(T)
  T E=rand()/(T)RAND_MAX,E0=E,E2=E;
  Vec variable=Vec::Random(getNumDOF());
  Vec dVariable=Vec::Random(getNumDOF());
  Vec G=Vec::Random(getNumDOF()),G0=G,G2=G;
  MatT H=MatT::Random(getNumDOF(),getNumDOF()),H0=H;
  smoothnessRegularizer(variable,&E,&G,&H,d);
  smoothnessBF(variable,&E2,segments,d);
  DEBUG_GRADIENT("E",E-E0,E-E2)

  E2=E0;
  smoothnessRegularizer(variable+dVariable*DELTA,&E2,&G2,NULL,d);
  DEBUG_GRADIENT("G",(G-G0).dot(dVariable),(G-G0).dot(dVariable)-(E2-E)/DELTA)
  DEBUG_GRADIENT("H",((H-H0)*dVariable).norm(),((H-H0)*dVariable-(G2-G)/DELTA).norm())
}
template <typename T>
void ThetaTrajectory<T>::debugEnergy(T t) const {
  Vec GTheta=Vec::Random(_dimension),GThetaVariable,GThetaVariable2;
  MatT HTheta=MatT::Random(_dimension,_dimension);
  HTheta=(HTheta+HTheta.transpose()).eval();

  DEFINE_NUMERIC_DELTA_T(T)
  Vec G=Vec::Random(getNumDOF()),G0=G,G2=G;
  MatT H=MatT::Random(getNumDOF(),getNumDOF()),H0=H;
  Vec variable=Vec::Random(getNumDOF());
  Vec dVariable=Vec::Random(getNumDOF());
  Vec theta=getPoint(variable,t);
  Vec theta2=getPoint(variable+dVariable*DELTA,t);
  GThetaVariable=HTheta*theta+GTheta;
  assembleEnergy(t,&G,&H,&GThetaVariable,&HTheta);
  GThetaVariable2=HTheta*theta2+GTheta;
  assembleEnergy(t,&G2,NULL,&GThetaVariable2,NULL);

  T E=theta.dot(GTheta)+theta.dot(HTheta*theta)/2;
  T E2=theta2.dot(GTheta)+theta2.dot(HTheta*theta2)/2;
  DEBUG_GRADIENT("G",(G-G0).dot(dVariable),(G-G0).dot(dVariable)-(E2-E)/DELTA)
  DEBUG_GRADIENT("H",((H-H0)*dVariable).norm(),((H-H0)*dVariable-(G2-G)/DELTA).norm())
}
//helper
template <typename T>
void ThetaTrajectory<T>::addNode(STrips& trips,int& k) {
  trips.push_back(STrip(k,_numNode,1));
  k++;
  _numNode++;
}
template <typename T>
void ThetaTrajectory<T>::addJointLeft(STrips& trips,int& k) const {
  trips.push_back(STrip(k,_numNode,1));
  k++;
  trips.push_back(STrip(k,_numNode,0.5));
  trips.push_back(STrip(k,_numNode+1,0.5));
  k++;
  trips.push_back(STrip(k,_numNode,0.25));
  trips.push_back(STrip(k,_numNode+1,0.5));
  trips.push_back(STrip(k,_numNode+2,0.25));
  k++;
}
template <typename T>
void ThetaTrajectory<T>::addJointRight(STrips& trips,int& k) {
  trips.push_back(STrip(k,_numNode,0.25));
  trips.push_back(STrip(k,_numNode+1,0.5));
  trips.push_back(STrip(k,_numNode+2,0.25));
  k++;
  trips.push_back(STrip(k,_numNode+1,0.5));
  trips.push_back(STrip(k,_numNode+2,0.5));
  k++;
  trips.push_back(STrip(k,_numNode+2,1));
  k++;
  _numNode+=3;
}
template <typename T>
typename ThetaTrajectory<T>::Vec ThetaTrajectory<T>::getMaxGrad(const Vec& variable,T t0,T t1) const {
  if(_segments.empty()) {
    ASSERT_MSG(t0==0 && t1==0,"StaticPose can only accept t0=t1=0!")
    return Vec::Zero(variable.size());
  }
  if(t1<t0)
    std::swap(t0,t1);
  int timeOffset=0;
  Vec maxThetaGrad=Vec::Constant(_dimension,-std::numeric_limits<double>::max());
  for(const auto& seg:_segments) {
    T segT0=std::min<T>(std::max<T>(t0-timeOffset,0),1);
    T segT1=std::min<T>(std::max<T>(t1-timeOffset,0),1);
    if(segT0<segT1) {
      Vec controlPointsFlatten=seg._variableToCPDim*variable;
      MatT controlPoints=Eigen::Map<MatT>(controlPointsFlatten.data(),_dimension,_order+1);
      maxThetaGrad=maxThetaGrad.cwiseMax(_curve.getMaxGrad(controlPoints,segT0,segT1));
    }
    timeOffset+=1;
  }
  return maxThetaGrad;
}
template <typename T>
void ThetaTrajectory<T>::debugMaxGrad(int res) const {
  Vec variable=Vec::Random(getNumDOF());
  T t0=rand()/(T)RAND_MAX*getNumSegment();
  T t1=rand()/(T)RAND_MAX*getNumSegment();
  if(t1<t0)
    std::swap(t0,t1);
  std::cout << "Testing maxGrad for time segment: [" << t0 << "," << t1 << "]" << std::endl;
  Vec maxGrad=getMaxGrad(variable,t0,t1);
  Vec maxGradRef=Vec::Constant(_dimension,0);
  for(int i=0; i<res; i++) {
    T alpha=(i+0.5)/res;
    maxGradRef=maxGradRef.cwiseMax(getDerivative(variable,interp1D(t0,t1,alpha),1).cwiseAbs());
  }
  for(int i=0; i<_dimension; i++) {
    std::cout << "maxGradRef[" << i << "]=" << maxGradRef[i] << " maxGrad[" << i << "]=" << maxGrad[i] << std::endl;
    ASSERT_MSGV(maxGradRef[i]<maxGrad[i]+Epsilon<T>::defaultEps(),"maxGradRef[%d](%f)>=maxGrad[%d](%f)",i,(double)maxGradRef[i],i,(double)maxGrad[i])
  }
}
//instance
template class ThetaTrajectory<FLOAT>;
}
