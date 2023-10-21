#ifndef BEZIER_CURVE_H
#define BEZIER_CURVE_H

#include <Utils/IO.h>
#include <Utils/Pragma.h>
#include <Utils/SparseUtils.h>
#include <Utils/ParallelVector.h>

namespace PHYSICSMOTION {
template <typename T>
class BezierCurve : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  BezierCurve();
  BezierCurve(int order,int dimension);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  MatT getBezierCoefficients() const;
  MatT getCoefficientMatrix() const;
  int getOrder() const;
  int getDimension() const;
  int getFactorial(int n) const;
  const SMatT& getLeftSubd() const;
  const SMatT& getRightSubd() const;
  void debugSubdivide() const;
  //monomial basis at time t
  Vec getBasis(T t) const;
  //flat variable at time t
  Vec getPoint(Eigen::Map<const MatT> controlPoints,T t) const;
  //time derivative of flat variable at time t of order d
  std::pair<BezierCurve,BezierCurve::MatT> getDerivative(int derivativeOrder) const;
  Vec getDerivative(Eigen::Map<const MatT> controlPoints,T t,int d) const;
  void debugDerivative(Eigen::Map<const MatT> controlPoints,T t,int d) const;
  void debugDerivative() const;
  //return the subdivision stencil mapping control points to the control points representing sub-trajectory [t0,t1]
  static SMatT getControlPoints(int order,T t0,T t1,bool debugOutput);
  //subdivide into line segments, return the ordered nodes
  std::vector<Vec> getSubdividedPoints(Eigen::Map<const MatT> controlPoints,T epsilon) const;
  //get integral hessian
  MatT getBasisCrossIntegral(int d) const;
  SMatT getBasisCrossIntegralFlatten(int d) const;
  //get maximum gradient
  Vec getMaxGrad(MatT controlPoints,T t0,T t1) const;
  void debugMaxGrad(int res) const;
 private:
  int _order;
  int _dimension;
  MatT _coefficientMatrix;
  SMatT _leftSubd,_rightSubd;
};
}

#endif
