#ifndef SPARSE_UTILS_H
#define SPARSE_UTILS_H

#include "Pragma.h"
#include "Epsilon.h"
#include "DebugGradient.h"
#include "ParallelVector.h"
#include <Eigen/Sparse>
#include <iostream>

namespace PHYSICSMOTION {
//sparse matrix building
//dense
template <typename MAT,typename Derived>
void addBlock(MAT& H,int r,int c,const Eigen::MatrixBase<Derived>& coef) {
  H.block(r,c,coef.rows(),coef.cols())+=coef;
}
template <typename MAT>
void addBlock(MAT& H,int r,int c,typename MAT::Scalar coef) {
  H(r,c)+=coef;
}
template <typename MAT,typename T>
void addBlockId(MAT& H,int r,int c,int nr,T coefId) {
  H.block(r,c,nr,nr).diagonal().array()+=coefId;
}
//sparse
template <typename Derived>
void addBlock(ParallelVector<Eigen::Triplet<typename Derived::Scalar,int>>& H,int r,int c,const Eigen::MatrixBase<Derived>& coef) {
  int nrR=coef.rows();
  int nrC=coef.cols();
  for(int i=0; i<nrR; i++)
    for(int j=0; j<nrC; j++)
      H.push_back(Eigen::Triplet<typename Derived::Scalar,int>(r+i,c+j,coef(i,j)));
}
template <typename T>
void addBlock(ParallelVector<Eigen::Triplet<T,int>>& H,int r,int c,T coef) {
  H.push_back(Eigen::Triplet<T,int>(r,c,coef));
}
template <typename Derived>
void addBlockI(ParallelVector<Eigen::Triplet<typename Derived::Scalar,int>>& H,int r,int c,const Eigen::MatrixBase<Derived>& coefI) {
  int nr=coefI.size();
  for(int i=0; i<nr; i++)
    H.push_back(Eigen::Triplet<typename Derived::Scalar,int>(r+i,c+i,coefI[i]));
}
template <typename T>
void addBlockId(ParallelVector<Eigen::Triplet<T,int>>& H,int r,int c,int nr,T coefId) {
  for(int i=0; i<nr; i++)
    H.push_back(Eigen::Triplet<T,int>(r+i,c+i,coefId));
}
template <typename T,int O,typename I>
void addBlock(ParallelVector<Eigen::Triplet<T,int>>& H,int r,int c,const Eigen::SparseMatrix<T,O,I>& coef) {
  for(int k=0; k<coef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(coef,k); it; ++it)
      H.push_back(Eigen::Triplet<T,int>(r+(int)it.row(),c+(int)it.col(),it.value()));
}
//sparseVec
template <typename Derived>
void addBlock(Eigen::SparseVector<typename Derived::Scalar,0,int>& H,int r,const Eigen::MatrixBase<Derived>& coef) {
  for(int i=0; i<coef.rows(); i++)
    H.coeffRef(r+i)+=coef[i];
}
//sparseBlk
template <typename T,int O,typename I>
Eigen::SparseMatrix<T,O,I> sparseBlk(const Eigen::SparseMatrix<T,O,I>& h,int r,int c,int nr,int nc) {
  ParallelVector<Eigen::Triplet<T,int>> trips;
  for(int k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it)
      if(it.row()>=r && it.row()<r+nr && it.col()>=c && it.col()<c+nc)
        trips.push_back(Eigen::Triplet<T,int>(it.row()-r,it.col()-c,it.value()));

  Eigen::SparseMatrix<T,O,I> ret;
  ret.resize(nr,nc);
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
//maxAbs
template <typename T,int r,int c>
T absMax(const Eigen::Matrix<T,r,c>& h) {
  return h.cwiseAbs().maxCoeff();
}
template <typename T,int O,typename I>
T absMax(const Eigen::SparseMatrix<T,O,I>& h) {
  T ret=0;
  for(int k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it)
      ret=std::max<T>(ret,fabs(it.value()));
  return ret;
}
template <typename T,int O,typename I>
T absMaxRel(const Eigen::SparseMatrix<T,O,I>& h,const Eigen::SparseMatrix<T,O,I>& hRef,bool detail) {
  int row=-1,col=-1;
  T ret=0,num=0,denom=0;
  T EPS=Epsilon<T>::defaultEps();
  //check against h
  for(int k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it) {
      T err=fabs(it.value()-hRef.coeff(it.row(),it.col()));
      T ref=fabs(hRef.coeff(it.row(),it.col()));
      T rel=err/std::max<T>(ref,EPS);
      if(err<EPS)
        continue;
      if(rel>ret) {
        num=err;
        denom=ref;
        row=it.row();
        col=it.col();
      }
      ret=std::max(ret,rel);
    }
  //check against hRef
  for(int k=0; k<hRef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(hRef,k); it; ++it) {
      T err=fabs(h.coeff(it.row(),it.col())-it.value());
      T ref=fabs(it.value());
      T rel=err/std::max<T>(ref,EPS);
      if(err<EPS)
        continue;
      if(rel>ret) {
        num=err;
        denom=ref;
        row=it.row();
        col=it.col();
      }
      ret=std::max(ret,rel);
    }
  if(detail)
    std::cout << "(" << row << "," << col << ") Num: " << num << " Denom: " << denom << std::endl;
  return ret;
}
//build KKT matrix
template <typename MAT>
MAT buildKKT(const MAT& h,const MAT& a,typename MAT::Scalar shift) {
  MAT kkt;
  kkt.setZero(h.rows()+a.rows(),h.rows()+a.rows());
  kkt.block(0,0,h.rows(),h.rows())=h;
  kkt.block(0,0,h.rows(),h.rows()).diagonal().array()+=shift;
  kkt.block(h.rows(),0,a.rows(),a.cols())=a;
  kkt.block(0,h.rows(),a.cols(),a.rows())=a.transpose();
  return kkt;
}
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> buildKKT(const Eigen::SparseMatrix<T,O,I>& h,const Eigen::SparseMatrix<T,O,I>& a,T shift,bool at=false,const Eigen::Matrix<T,-1,1>* shiftOffDiag=NULL) {
  typedef Eigen::SparseMatrix<T,O,I> SMat;
  SMat kkt;
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(int k=0; k<h.outerSize(); ++k)
    for(typename SMat::InnerIterator it(h,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>((int)it.row(),(int)it.col(),it.value()));
  if(shift!=0)
    for(int k=0; k<h.rows(); ++k)
      trips.push_back(typename Eigen::Triplet<T,I>(k,k,shift));
  for(int k=0; k<a.outerSize(); ++k)
    for(typename SMat::InnerIterator it(a,k); it; ++it)
      if(at) {
        trips.push_back(typename Eigen::Triplet<T,I>(h.rows()+it.col(),it.row(),it.value()));
        trips.push_back(typename Eigen::Triplet<T,I>(it.row(),h.rows()+it.col(),it.value()));
      } else {
        trips.push_back(typename Eigen::Triplet<T,I>((int)h.rows()+(int)it.row(),(int)it.col(),it.value()));
        trips.push_back(typename Eigen::Triplet<T,I>((int)it.col(),(int)h.rows()+(int)it.row(),it.value()));
      }
  int dualDim=at?a.cols():a.rows();
  if(shiftOffDiag) {
    ASSERT_MSG(shiftOffDiag->size()==dualDim,"BuiltKKT has incorrect shiftOffDiag size")
    for(int i=0; i<dualDim; i++)
      trips.push_back(typename Eigen::Triplet<T,I>((int)h.rows()+i,(int)h.cols()+i,shiftOffDiag->coeff(i)));
  }
  kkt.resize(h.rows()+dualDim,h.rows()+dualDim);
  kkt.setFromTriplets(trips.begin(),trips.end());
  return kkt;
}
//kronecker-product
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> kronecker(const Eigen::SparseMatrix<T,O,I>& h,int n) {
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(int k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it)
      for(int d=0; d<n; d++)
        trips.push_back(typename Eigen::Triplet<T,I>(it.row()*n+d,it.col()*n+d,it.value()));

  Eigen::SparseMatrix<T,O,I> ret;
  ret.resize(h.rows()*n,h.cols()*n);
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
//concat-diag
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> concatDiag(const Eigen::SparseMatrix<T,O,I>& a,const Eigen::SparseMatrix<T,O,I>& b) {
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(int k=0; k<a.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(a,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>(it.row(),it.col(),it.value()));
  for(int k=0; k<b.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(b,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>(a.rows()+it.row(),a.cols()+it.col(),it.value()));

  typename Eigen::SparseMatrix<T,O,I> ret;
  ret.resize(a.rows()+b.rows(),a.cols()+b.cols());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
//identity
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> identity(int r,int c) {
  ParallelVector<Eigen::Triplet<T,int>> trips;
  for(int i=0; i<std::min(r,c); i++)
    trips.push_back(Eigen::Triplet<T,int>(i,i,1));
  Eigen::SparseMatrix<T,O,I> ID;
  ID.resize(r,c);
  ID.setFromTriplets(trips.begin(),trips.end());
  return ID;
}
template <typename T,int O,typename I,typename MT>
Eigen::SparseMatrix<T,O,I> toSparse(const MT& m,T eps=0) {
  Eigen::SparseMatrix<T,O,I> ret;
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(int r=0; r<m.rows(); r++)
    for(int c=0; c<m.cols(); c++)
      if(fabs(m(r,c))>eps)
        trips.push_back(Eigen::Triplet<T,I>(r,c,m(r,c)));
  ret.resize(m.rows(),m.cols());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T,int O,typename I>
Eigen::SparseMatrix<T,O,I> concat(const Eigen::SparseMatrix<T,O,I>& A,const Eigen::SparseMatrix<T,O,I>& B) {
  ASSERT(A.cols()==B.cols())
  Eigen::SparseMatrix<T,O,I> M(A.rows()+B.rows(),A.cols());
  M.reserve(A.nonZeros()+B.nonZeros());
  for(int c=0; c<A.cols(); ++c) {
    M.startVec(c);
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator itA(A,c); itA; ++itA)
      M.insertBack(itA.row(),c)=itA.value();
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator itB(B,c); itB; ++itB)
      M.insertBack(itB.row()+A.rows(),c)=itB.value();
  }
  M.finalize();
  return M;
}
template <typename VECA,typename VECB>
Eigen::Matrix<typename VECA::Scalar,-1,1> concat(const VECA& A,const VECB& B) {
  Eigen::Matrix<typename VECA::Scalar,-1,1> ret;
  ret.resize(A.size()+B.size());
  ret << A,B;
  return ret;
}
template <typename MATA,typename MATB>
Eigen::Matrix<typename MATA::Scalar,-1,-1> concatRow(const MATA& A,const MATB& B) {
  Eigen::Matrix<typename MATA::Scalar,-1,-1> ret;
  ret.resize(A.rows()+B.rows(),A.cols());
  ret << A,B;
  return ret;
}
template <typename MATA,typename MATB>
Eigen::Matrix<typename MATA::Scalar,-1,-1> concatCol(const MATA& A,const MATB& B) {
  Eigen::Matrix<typename MATA::Scalar,-1,-1> ret;
  ret.resize(A.rows(),A.cols()+B.cols());
  ret << A,B;
  return ret;
}
//3x3 multiply Left/Right flatten
template <typename T>
Eigen::SparseMatrix<T,0,int> mat3x3MulRight(const Eigen::Matrix<T,3,3>& m) {
  typedef ParallelVector<Eigen::Triplet<T,int> > TRIPS;
  Eigen::SparseMatrix<T,0,int> M;
  M.resize(9,9);
  TRIPS trips;
  for(int r=0; r<3; r++)
    for(int c=0; c<3; c++)
      for(int k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,int>(r+c*3,r+k*3,m(k,c)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,int> mat3x3MulTRight(const Eigen::Matrix<T,3,3>& m) {
  typedef ParallelVector<Eigen::Triplet<T,int> > TRIPS;
  Eigen::SparseMatrix<T,0,int> M;
  M.resize(9,9);
  TRIPS trips;
  for(int r=0; r<3; r++)
    for(int c=0; c<3; c++)
      for(int k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,int>(r+c*3,k+r*3,m(k,c)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,int> mat3x3MulLeft(const Eigen::Matrix<T,3,3>& m) {
  typedef ParallelVector<Eigen::Triplet<T,int> > TRIPS;
  Eigen::SparseMatrix<T,0,int> M;
  M.resize(9,9);
  TRIPS trips;
  for(int r=0; r<3; r++)
    for(int c=0; c<3; c++)
      for(int k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,int>(r+c*3,k+c*3,m(r,k)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,int> mat3x3MulTLeft(const Eigen::Matrix<T,3,3>& m) {
  typedef ParallelVector<Eigen::Triplet<T,int> > TRIPS;
  Eigen::SparseMatrix<T,0,int> M;
  M.resize(9,9);
  TRIPS trips;
  for(int r=0; r<3; r++)
    for(int c=0; c<3; c++)
      for(int k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,int>(r+c*3,k+c*3,m(k,r)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
}

#endif
