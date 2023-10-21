#include <Utils/SparseUtils.h>
#include <Utils/SparseVisualizer.h>
#include <fstream>
#include <random>

using namespace PHYSICSMOTION;

template <typename T>
void debugSparse(bool vis,int N=10,int M=3) {
  std::cout << "-------------------------------------------------------------DebugSparse" << std::endl;
  //fill-in entries
  ParallelVector<Eigen::Triplet<T,int>> trips;
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      if(std::abs(i-j)<M) {
        Eigen::Matrix<T,3,3> blk;
        blk.setRandom();
        addBlock(trips,i*3,j*3,blk);
      }
  //build matrix
  Eigen::SparseMatrix<T,0,int> m,a;
  m.resize(N*3,N*3);
  m.setFromTriplets(trips.begin(),trips.end());
  //extract block
  m=sparseBlk(m,M,M,N*M-M,N*M-M);
  //do kronecker
  m=kronecker(m,2);
  //build KKT
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1,1);
  a.resize(m.rows(),m.cols());
  for(int i=0; i<a.rows(); i++)
    for(int j=0; j<a.cols(); j++)
      if(dis(gen)<0)
        a.coeffRef(i,j)+=1;
  m=buildKKT<T,0,int>(m,a,0);
  //concat
  a=concatDiag(m,a);
  a=toSparse<T,0,int>(a.toDense());
  a=concat(a,a);
  //visualize
  if(vis)
    visualizeSparse(0,NULL,a);
}
int main(int argc,char** argv) {
  mpfr_float::default_precision(100);
  debugSparse<float>(false);
  debugSparse<double>(false);
  debugSparse<float128>(false);
  debugSparse<mpfr_float>(true);
  return 0;
}
