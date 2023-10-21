#ifndef DEFORMED_ENVIRONMENT_H
#define DEFORMED_ENVIRONMENT_H

#include "Environment.h"
#include <Utils/SparseUtils.h>

namespace PHYSICSMOTION {
template <typename T>
class DeformedEnvironment : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  struct NodeFetcher {
    T& operator[](int xyz);
    const T& operator[](int xyz) const;
    const Eigen::Matrix<int,3,1>& _stride;
    int _offset,_DIM;
    Vec& _nodes;
  };
  DeformedEnvironment();
  DeformedEnvironment(const Eigen::Matrix<int,3,1>& nrCell,std::shared_ptr<Environment<T>> env);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  void setHistoryPath(const std::string& path);
  std::vector<Vec3T> getIsoline(int dir,T x,T y) const;
  void writeVTK(const std::string& path,int res=3,bool bottomOnly=false) const;
  void writeFrameVTK(const std::string& path,T len=0.05,int res=3) const;
  void optimizeLocalGlobal(int heightRes=8,int volumeRes=8,T wACAP=0.5,T tol=1e-3,bool writeHistory=true);
  void optimizeLocalInjective(int heightRes=8,int volumeRes=8,T wACAP=0.5,T mu=1e-3,T tol=1e-1,bool writeHistory=true);
  void initializeNodes();
  void debug(int N=10);
  Vec3T nrCell() const;
  std::shared_ptr<Environment<T>> getOptimizedEnvironment(int res=3) const;
  T localInjectiveEnergy(int heightRes,int volumeRes,T wACAP,T mu,Vec* g,SMatT* h=NULL);
  Mat3T rotation(const Vec3T& alpha,Mat3T* w,Vec3T* posOut=NULL,Mat3T* JposOut=NULL) const;
  Eigen::Matrix<int,3,1> aFloorSafe(const Vec3T& a) const;
  void jacobian(int row,const Vec3T& a,STrips& trips,int DIM=-1) const;
  void jacobianJ(int row,const Vec3T& a,STrips& trips,T wACAP=1) const;
  SMatT jacobian(const Vec3T& a) const;
  SMatT jacobianJ(const Vec3T& a,T wACAP=1) const;
  Vec3T alpha(const Vec3T& pos,const T tol=1e-3) const;
  Vec3T forward(const Vec3T& a,Mat3T* J,Mat3T* hess) const;
  void backward(const Vec3T& a,const Vec3T& dEde,const Mat3T& dEdJ);
 private:
  static T forward(const Vec3T& a,NodeFetcher n,Vec3T* grad,Mat3T* hess);
  static void backward(const Vec3T& a,NodeFetcher n,T dEde,const Vec3T& dEdJ);
  static MatT DSDJ(const Eigen::JacobiSVD<Eigen::Matrix<double,3,3>>& SVD);
  //data
  std::string _historyPath;
  std::shared_ptr<Environment<T>> _env;
  Eigen::Matrix<int,3,1> _nrCell,_stride;
  Vec _nodes,_Dnodes;
  SMatT _fixBasis;
};
}

#endif
