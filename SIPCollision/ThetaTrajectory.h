#ifndef THETA_TRAJECTORY_H
#define THETA_TRAJECTORY_H

#include "BezierCurve.h"
#include <Utils/ParallelVector.h>

namespace PHYSICSMOTION {
template <typename T>
class ThetaTrajectory : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  typedef std::function<void(int,const Vec&)> GFunc;
  typedef std::function<void(int,int,const MatT&)> HFunc;
  struct Segment : public SerializableBase {
    virtual bool read(std::istream& is,IOData* dat) override;
    virtual bool write(std::ostream& os,IOData* dat) const override;
    virtual std::shared_ptr<SerializableBase> copy() const override;
    virtual std::string type() const override;
    SMatT _variableToCP;
    SMatT _variableToCPDim;
    int _index;
  };
  ThetaTrajectory()=default;
  ThetaTrajectory(int dimension,int order,int numSegment);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  const std::vector<Segment>& getSegments() const;
  std::vector<Segment>& getSegments();
  const BezierCurve<T>& getCurve() const;
  int getNumSegment() const;
  int getDimension() const;
  int getNumNode() const;
  int getNumDOF() const;
  int getOrder() const;
  T getLength(const Vec& variable,T interval=1e-3) const;
  Vec getPoint(const Vec& variable,T t) const;
  Vec getDerivative(const Vec& variable,T t,int d=1) const;
  void getControlPointCoeff(STrips& trips,int rowOff,int segment,int dimension,int cpId) const;
  SMatT getControlPointCoeffAll(int dimension) const;
  void debugControlPointAllDifference() const;
  void assembleEnergy(T time,GFunc* G,HFunc* H,const Vec* GTheta,const MatT* HTheta) const;
  void assembleEnergy(T time,Vec* G,MatT* H,const Vec* GTheta,const MatT* HTheta) const;
  void smoothnessRegularizer(const Vec& variable,T* E,Vec* G,MatT* H,int d=1) const;
  void smoothnessBF(const Vec& variable,T* E,int segments,int d=1) const;
  void debugSmoothness(int segments,int d) const;
  void debugEnergy(T t) const;
  //get maximum gradient
  Vec getMaxGrad(const Vec& controlPoints,T t0,T t1) const;
  void debugMaxGrad(int res) const;
 protected:
  void addNode(STrips& trips,int& k);
  void addJointLeft(STrips& trips,int& k) const;
  void addJointRight(STrips& trips,int& k);
  //data
  BezierCurve<T> _curve;
  std::vector<Segment> _segments;
  int _numSegment,_numNode;
  int _order,_dimension;
};
}

#endif
