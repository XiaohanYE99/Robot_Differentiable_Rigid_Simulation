#include "BezierCurve.h"
#include <Utils/Interp.h>

namespace PHYSICSMOTION {
template <typename T>
BezierCurve<T>::BezierCurve() {}
template <typename T>
BezierCurve<T>::BezierCurve(int order,int dimension) {
  _order=order;
  _dimension=dimension;
  _coefficientMatrix=getBezierCoefficients();
  _leftSubd=BezierCurve::getControlPoints(_order,0,0.5,false);
  _rightSubd=BezierCurve::getControlPoints(_order,0.5,1,false);
}
template <typename T>
typename BezierCurve<T>::MatT BezierCurve<T>::getBezierCoefficients() const {
  MatT coefficientMatrix(_order+1,_order+1);
  for(int j=0; j<_order+1; j++) {
    T coeff=getFactorial(_order)/getFactorial(_order-j);
    for(int i=0; i<_order+1; i++)
      if(i<=j)
        coefficientMatrix(i,j)=coeff*std::pow(-1,i+j)/(getFactorial(i)*getFactorial(j-i));
      else coefficientMatrix(i,j)=0;
  }
  return coefficientMatrix;
}
template <typename T>
bool BezierCurve<T>::read(std::istream& is,IOData*) {
  readBinaryData(_order,is);
  readBinaryData(_dimension,is);
  readBinaryData(_coefficientMatrix,is);
  readBinaryData(_leftSubd,is);
  readBinaryData(_rightSubd,is);
  return is.good();
}
template <typename T>
bool BezierCurve<T>::write(std::ostream& os,IOData*) const {
  writeBinaryData(_order,os);
  writeBinaryData(_dimension,os);
  writeBinaryData(_coefficientMatrix,os);
  writeBinaryData(_leftSubd,os);
  writeBinaryData(_rightSubd,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> BezierCurve<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new BezierCurve);
}
template <typename T>
std::string BezierCurve<T>::type() const {
  return typeid(BezierCurve).name();
}
template <typename T>
typename BezierCurve<T>::MatT BezierCurve<T>::getCoefficientMatrix() const {
  return _coefficientMatrix;
}
template <typename T>
int BezierCurve<T>::getOrder() const {
  return _order;
}
template <typename T>
int BezierCurve<T>::getDimension() const {
  return _dimension;
}
template <typename T>
int BezierCurve<T>::getFactorial(int n) const {
  ASSERT(n >= 0)
  if(n==0)
    return 1;
  return n*getFactorial(n-1);
}
template <typename T>
const typename BezierCurve<T>::SMatT& BezierCurve<T>::getLeftSubd() const {
  return _leftSubd;
}
template <typename T>
const typename BezierCurve<T>::SMatT& BezierCurve<T>::getRightSubd() const {
  return _rightSubd;
}
template <typename T>
void BezierCurve<T>::debugSubdivide() const {
  Vec controlPoints=Vec::Random((_order+1)*_dimension);
  Eigen::Map<const MatT> CPs(controlPoints.data(),_dimension,_order+1);
  MatT CPLs=CPs*getLeftSubd().transpose();
  MatT CPRs=CPs*getRightSubd().transpose();
  T t=rand()/(T)RAND_MAX;
  DEFINE_NUMERIC_DELTA_T(T)
  DEBUG_GRADIENT("subdivideLeft",getPoint(CPs,t).norm(),(getPoint(CPs,t)-getPoint(Eigen::Map<const MatT>(CPLs.data(),CPLs.rows(),CPLs.cols()),t*2)).norm())
  DEBUG_GRADIENT("subdivideRight",getPoint(CPs,t).norm(),(getPoint(CPs,t)-getPoint(Eigen::Map<const MatT>(CPRs.data(),CPRs.rows(),CPRs.cols()),t*2-1)).norm())
}
//monomial basis at time t
template <typename T>
typename BezierCurve<T>::Vec BezierCurve<T>::getBasis(T t) const {
  T tPow=1;
  Vec tvec(_order+1);
  for(int j=0; j<_order+1; j++) {
    tvec(j)=tPow;
    tPow*=t;
  }
  return tvec;
}
//flat variable at time t
template <typename T>
typename BezierCurve<T>::Vec BezierCurve<T>::getPoint(Eigen::Map<const MatT> controlPoints,T t) const {
  return controlPoints*_coefficientMatrix*getBasis(t);
}
//time derivative of flat variable at time t of order d
template <typename T>
std::pair<BezierCurve<T>,typename BezierCurve<T>::MatT> BezierCurve<T>::getDerivative(int derivativeOrder) const {
  int prefix=1;
  for(int i=0; i<derivativeOrder; ++i)
    prefix*=_order-i;
  MatT controlPointMapAll=MatT::Identity(_order+1,_order+1);
  for(int k=0; k<derivativeOrder; ++k) {
    for(int i=0; i<_order-k; ++i) {
      MatT controlPointMap=MatT::Identity(_order+1,_order+1);
      controlPointMap(i+1,i)=1;
      controlPointMap(i,i)=-1;
      controlPointMapAll*=controlPointMap;
    }
  }
  BezierCurve derivativedCurve(_order-derivativeOrder,_dimension);
  return std::make_pair(derivativedCurve,controlPointMapAll.block(0,0,_order+1,_order+1-derivativeOrder)*prefix);
}
template <typename T>
typename BezierCurve<T>::Vec BezierCurve<T>::getDerivative(Eigen::Map<const MatT> controlPoints,T t,int derivativeOrder) const {
  int prefix=1;
  for(int i=0; i<derivativeOrder; ++i)
    prefix*=_order-i;
  MatT tempControlPoints=controlPoints;
  for(int k=0; k<derivativeOrder; ++k)
    for(int i=0; i<_order-k; ++i)
      tempControlPoints.col(i)=tempControlPoints.col(i+1)-tempControlPoints.col(i);
  BezierCurve derivativedCurve(_order-derivativeOrder,_dimension);
  Vec derivative=derivativedCurve.getPoint(Eigen::Map<const MatT>(tempControlPoints.data(),tempControlPoints.rows(),tempControlPoints.cols()-derivativeOrder),t);
  return prefix*derivative;
}
template <typename T>
void BezierCurve<T>::debugDerivative(Eigen::Map<const MatT> controlPoints,T t,int d) const {
  DEFINE_NUMERIC_DELTA_T(T)
  Vec point=getDerivative(controlPoints,t,d);
  Vec point2=getDerivative(controlPoints,t+DELTA,d);
  Vec dPoint=getDerivative(controlPoints,t,d+1);
  if(d==0) {
    Vec pointRef=getPoint(controlPoints,t);
    DEBUG_GRADIENT("point",pointRef.norm(),(pointRef-point).norm())
  } else {
    std::pair<BezierCurve,MatT> curve=getDerivative(d);
    MatT controlPointsDerivative=controlPoints*curve.second;
    Vec pointRef=curve.first.getPoint(Eigen::Map<const MatT>(controlPointsDerivative.data(),controlPointsDerivative.rows(),controlPointsDerivative.cols()),t);
    DEBUG_GRADIENT("point",pointRef.norm(),(pointRef-point).norm())
  }
  DEBUG_GRADIENT("derivative",dPoint.norm(),(dPoint-(point2-point)/DELTA).norm())
}
template <typename T>
void BezierCurve<T>::debugDerivative() const {
  T t=rand()/(T)RAND_MAX;
  Vec controlPoints=Vec::Random((_order+1)*_dimension);
  for(int d=0; d<_order; d++)
    debugDerivative(Eigen::Map<const MatT>(controlPoints.data(),_dimension,_order+1),t,d);
}
//return the subdivision stencil mapping control points to the control points representing sub-trajectory [t0,t1]
template <typename T>
typename BezierCurve<T>::SMatT BezierCurve<T>::getControlPoints(int order,T t0,T t1,bool debugOutput) {
  //combination
  ASSERT_MSG((t0>=0&&t1<=1), "Wrong time slot!")
  std::vector<std::vector<long int>> combination;
  combination.resize(order+1);  //largest order
  combination[0].resize(1);
  combination[0][0]=1;
  for(int i=1; i<=order; ++i) {
    combination[i].resize(i+1);
    long long temp=1;
    for(int j=0; j<=i; ++j) {
      combination[i][j]=temp;
      //std::cout << temp << " ";
      temp=temp*(i-j)/(j+1);
    }
    //std::cout << std::endl;
  }
  if(debugOutput)
    std::cout << "Combination size is: " << combination.size() << std::endl;

  //blossom
  std::vector<T> list_t0;
  std::vector<T> list_t1;
  std::vector<T> list_1_t0;
  std::vector<T> list_1_t1;
  list_t0.resize(order+1);
  list_t1.resize(order+1);
  list_1_t0.resize(order+1);
  list_1_t1.resize(order+1);
  T _t0=1;
  T _t1=1;
  T _1_t0=1;
  T _1_t1=1;
  for(int i=0; i<=order; ++i) {
    list_t0[i]=_t0;
    _t0*=t0;
    list_1_t0[i]=_1_t0;
    _1_t0*=1-t0;
    list_t1[i]=_t1;
    _t1*=t1;
    list_1_t1[i]=_1_t1;
    _1_t1*=1-t1;
  }

  MatT coeff;
  coeff.setZero(order+1,order+1);
  //row is t0 t1 number change,first ORDER t0 end ORDER t1
  //col is j-th control point
  for(int i=0; i<=order; ++i)
    for(int j=0; j<=order; ++j) {
      int maxk;
      if(i+j<order) {
        maxk=std::min(i,j);
        for(int k=0; k<=maxk; ++k) {
          coeff(i,j)+=combination[order-i][j-k]*combination[i][k] *
                      list_1_t0[order-i-j+k]*list_1_t1[i-k]*list_t0[j-k]*list_t1[k];
        }
      } else {
        maxk=std::min(order-i,order-j);
        for(int k=0; k<=maxk; ++k) {
          coeff(i,j)+=combination[order-i][k]*combination[i][order-j-k] *
                      list_1_t0[k]*list_1_t1[order-j-k]*list_t0[order-i-k]*list_t1[i+j-order+k];
        }
      }
    }
  if(debugOutput)
    std::cout << "Coeff is: \n" << coeff.transpose() << std::endl;
  return coeff.sparseView();
}
//subdivide into line segments, return the ordered nodes
template <typename T>
std::vector<typename BezierCurve<T>::Vec> BezierCurve<T>::getSubdividedPoints(Eigen::Map<const MatT> controlPoints,T epsilon) const {
  T lenPoly=0,len=(controlPoints.col(0)-controlPoints.col(_order)).norm();
  for(int i=0; i<_order; i++)
    lenPoly+=(controlPoints.col(i)-controlPoints.col(i+1)).norm();
  if(lenPoly<len*(1+epsilon))
    return std::vector<Vec>(1,controlPoints.col(_order));
  else {
    MatT CL=controlPoints*_leftSubd.transpose();
    std::vector<Vec> left=getSubdividedPoints(Eigen::Map<const MatT>(CL.data(),_dimension,_order+1),epsilon);
    MatT CR=controlPoints*_rightSubd.transpose();
    std::vector<Vec> right=getSubdividedPoints(Eigen::Map<const MatT>(CR.data(),_dimension,_order+1),epsilon);
    left.insert(left.end(),right.begin(),right.end());
    return left;
  }
}
//get integral hessian
template <typename T>
typename BezierCurve<T>::MatT BezierCurve<T>::getBasisCrossIntegral(int d) const {
  /*
  *Calculate the integral of the following integrand from 0 to 1
  * _  _
  *| 1  |
  *| t  |
  *| t^2|
  *| t^3|*[1,t,t^2,t^3,...,t^n]
  *| .  |
  *| .  |
  *| .  |
  *| t^n|
  * -  -
   */
  if(d>0) {
    std::pair<BezierCurve,MatT> curve=getDerivative(d);
    return curve.second*curve.first.getBasisCrossIntegral(0)*curve.second.transpose();
  } else {
    MatT mat=MatT::Zero(_order+1,_order+1);
    for(int i=0; i<_order+1; ++i) {
      for(int j=0; j<_order+1; ++j) {
        mat(i,j)=1./(i+1+j);
      }
    }
    return _coefficientMatrix*mat*_coefficientMatrix.transpose();
  }
}
template <typename T>
typename BezierCurve<T>::SMatT BezierCurve<T>::getBasisCrossIntegralFlatten(int d) const {
  return kronecker<T,0,int>(getBasisCrossIntegral(d).sparseView(),_dimension);
}
//get maximum gradient
template <typename T>
typename BezierCurve<T>::Vec BezierCurve<T>::getMaxGrad(MatT controlPoints,T t0,T t1) const {
  if(t1<t0)
    std::swap(t0,t1);
  controlPoints*=getControlPoints(_order,t0,t1,false).transpose();
  std::pair<BezierCurve,MatT> curve=getDerivative(1);
  MatT derivative=controlPoints*curve.second/(t1-t0);
  //return row-wise maximum
  Vec ret=Vec::Zero(derivative.rows());
  for(int i=0; i<ret.size(); i++)
    ret[i]=derivative.row(i).cwiseAbs().maxCoeff();
  return ret;
}
template <typename T>
void BezierCurve<T>::debugMaxGrad(int res) const {
  MatT controlPoints=MatT::Random(_dimension,_order+1);
  Eigen::Map<const MatT> controlPointsMap(controlPoints.data(),controlPoints.rows(),controlPoints.cols());
  T t0=rand()/(T)RAND_MAX;
  T t1=rand()/(T)RAND_MAX;
  if(t1<t0)
    std::swap(t0,t1);
  std::cout << "Testing maxGrad for time segment: [" << t0 << "," << t1 << "]" << std::endl;
  Vec maxGrad=getMaxGrad(controlPoints,t0,t1);
  Vec maxGradRef=Vec::Constant(_dimension,0);
  for(int i=0; i<res; i++) {
    T alpha=(i+0.5)/res;
    maxGradRef=maxGradRef.cwiseMax(getDerivative(controlPointsMap,interp1D(t0,t1,alpha),1).cwiseAbs());
  }
  for(int i=0; i<_dimension; i++) {
    std::cout << "maxGradRef[" << i << "]=" << maxGradRef[i] << " maxGrad[" << i << "]=" << maxGrad[i] << std::endl;
    ASSERT_MSGV(maxGradRef[i]<maxGrad[i]+Epsilon<T>::defaultEps(),"maxGradRef[%d](%f)>=maxGrad[%d](%f)",i,(double)maxGradRef[i],i,(double)maxGrad[i])
  }
}
//instance
template class BezierCurve<FLOAT>;
}
