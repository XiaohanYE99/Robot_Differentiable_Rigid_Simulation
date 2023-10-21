#include "PBDMatrixSolver.h"

namespace PHYSICSMOTION {
//PBDMatrixSolver
PBDMatrixSolver::PBDMatrixSolver(std::shared_ptr<ArticulatedBody> body):_body(body) {}
PBDMatrixSolver::~PBDMatrixSolver() {}
//PBDMatrixSolverEigen
PBDMatrixSolverEigen::PBDMatrixSolverEigen(std::shared_ptr<ArticulatedBody> body):PBDMatrixSolver(body) {}
void PBDMatrixSolverEigen::compute(const Eigen::MatrixBase<MatT>& h) {
  _ldlt.compute(h.template cast<double>());
}
PBDMatrixSolverEigen::MatT PBDMatrixSolverEigen::solve(const Eigen::MatrixBase<MatT>& b) const {
  return _ldlt.solve(b.template cast<double>()).template cast<T>();
}
//PBDMatrixSolverCRBA
PBDMatrixSolverCRBA::PBDMatrixSolverCRBA(std::shared_ptr<ArticulatedBody> body):PBDMatrixSolver(body) {
  _info.reset(*_body);
}
void PBDMatrixSolverCRBA::compute(const Eigen::MatrixBase<MatT>& h) {
  ASSERT_MSG(_info._invHM.size()==h.size(),"Body's input matrix size != input matrix size!")
  _info._invHM=h;
  _info.LTDL();
}
PBDMatrixSolverCRBA::MatT PBDMatrixSolverCRBA::solve(const Eigen::MatrixBase<MatT>& b) const {
  MatT invMb=b;
  _info.LTDLSolve<MatT&>(invMb);
  return invMb;
}
//PBDMatrixSolverABA
PBDMatrixSolverABA::PBDMatrixSolverABA(std::shared_ptr<ArticulatedBody> body):PBDMatrixSolver(body) {
  _info.reset(*_body);
}
void PBDMatrixSolverABA::compute(const Eigen::MatrixBase<Vec>& q,const Eigen::MatrixBase<Mat3XT>& MRR,const Eigen::MatrixBase<Mat3XT>& MRt,const Eigen::MatrixBase<Mat3XT>& MtR,const Eigen::MatrixBase<Mat3XT>& Mtt,const Eigen::MatrixBase<Vec>& d) {
  Mat3T R;
  ASSERT_MSG(MRR.size()==3*3*_body->nrJ(),"MRR's input matrix size != input matrix size!")
  ASSERT_MSG(MRt.size()==3*3*_body->nrJ(),"MRt's input matrix size != input matrix size!")
  ASSERT_MSG(MtR.size()==3*3*_body->nrJ(),"MtR's input matrix size != input matrix size!")
  ASSERT_MSG(Mtt.size()==3*3*_body->nrJ(),"Mtt's input matrix size != input matrix size!")
  ASSERT_MSG(_info._qM.size()==d.size(),"d's input matrix size != input matrix size!")
  _info.reset(*_body,q,Vec::Zero(_body->nrDOF()));
  _I.resize(6,_body->nrJ()*6);
  for(int k=0; k<_body->nrJ(); k++) {
    R=ROT(_info.NEArticulatedGradientInfoMap<T>::getTrans(k));
    _I.template block<3,3>(0+0,k*6+0)=R.transpose()*MRR.template block<3,3>(0,k*3)*R;
    _I.template block<3,3>(0+0,k*6+3)=R.transpose()*MRt.template block<3,3>(0,k*3)*R;
    _I.template block<3,3>(0+3,k*6+0)=R.transpose()*MtR.template block<3,3>(0,k*3)*R;
    _I.template block<3,3>(0+3,k*6+3)=R.transpose()*Mtt.template block<3,3>(0,k*3)*R;
  }
  _d=d;
}
void PBDMatrixSolverABA::compute(const Eigen::MatrixBase<Vec>& q,const Eigen::MatrixBase<Mat6XT>& I,const Eigen::MatrixBase<Vec>& d) {
  ASSERT_MSG(_info._IM.size()==I.size(),"I's input matrix size != input matrix size!")
  ASSERT_MSG(_info._qM.size()==d.size(),"d's input matrix size != input matrix size!")
  _info.reset(*_body,q,Vec::Zero(_body->nrDOF()));
  _I=I;
  _d=d;
}
void PBDMatrixSolverABA::compute(const Eigen::MatrixBase<MatT>&) {
  FUNCTION_NOT_IMPLEMENTED
}
PBDMatrixSolverABA::MatT PBDMatrixSolverABA::solve(const Eigen::MatrixBase<MatT>& b) const {
  Vec bcol;
  MatT ret=b;
  GradInfo info=_info;
  for(int c=0; c<b.cols(); c++) {
    bcol=b.col(c);
    info.ABAInner(*_body,
                  GradInfo::mapCV((const Vec6T*)NULL),
                  GradInfo::mapCM((const Mat6XT*)NULL),
                  GradInfo::mapCV(bcol),GradInfo::mapCM(_I),GradInfo::mapCV(_d));
    ret.col(c)=info._ddqM;
  }
  return ret;
}
}
