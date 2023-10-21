#ifndef PBD_MATRIX_SOLVER_H
#define PBD_MATRIX_SOLVER_H

#include <Articulated/ArticulatedBody.h>
#include <Articulated/NEArticulatedGradientInfo.h>

namespace PHYSICSMOTION {
class PBDMatrixSolver {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  PBDMatrixSolver(std::shared_ptr<ArticulatedBody> body);
  virtual ~PBDMatrixSolver();
  virtual void compute(const Eigen::MatrixBase<MatT>& h)=0;
  virtual MatT solve(const Eigen::MatrixBase<MatT>& b) const=0;
 protected:
  std::shared_ptr<ArticulatedBody> _body;
};
class PBDMatrixSolverEigen : public PBDMatrixSolver {
 public:
  PBDMatrixSolverEigen(std::shared_ptr<ArticulatedBody> body);
  void compute(const Eigen::MatrixBase<MatT>& h) override;
  MatT solve(const Eigen::MatrixBase<MatT>& b) const override;
 protected:
  Eigen::LDLT<Eigen::Matrix<double,-1,-1>> _ldlt;
};
class PBDMatrixSolverCRBA : public PBDMatrixSolver {
 public:
  typedef NEArticulatedGradientInfo<T> GradInfo;
  PBDMatrixSolverCRBA(std::shared_ptr<ArticulatedBody> body);
  void compute(const Eigen::MatrixBase<MatT>& h) override;
  MatT solve(const Eigen::MatrixBase<MatT>& b) const override;
 protected:
  GradInfo _info;
};
class PBDMatrixSolverABA : public PBDMatrixSolver {
 public:
  typedef NEArticulatedGradientInfo<T> GradInfo;
  PBDMatrixSolverABA(std::shared_ptr<ArticulatedBody> body);
  void compute(const Eigen::MatrixBase<Vec>& q,const Eigen::MatrixBase<Mat3XT>& MRR,const Eigen::MatrixBase<Mat3XT>& MRt,const Eigen::MatrixBase<Mat3XT>& MtR,const Eigen::MatrixBase<Mat3XT>& Mtt,const Eigen::MatrixBase<Vec>& d);
  void compute(const Eigen::MatrixBase<Vec>& q,const Eigen::MatrixBase<Mat6XT>& I,const Eigen::MatrixBase<Vec>& d);
  void compute(const Eigen::MatrixBase<MatT>& h) override;
  MatT solve(const Eigen::MatrixBase<MatT>& b) const override;
 protected:
  GradInfo _info;
  Mat6XT _I;
  Vec _d;
};
}

#endif
