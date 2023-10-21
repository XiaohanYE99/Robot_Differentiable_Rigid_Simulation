#include <Utils/IO.h>

using namespace PHYSICSMOTION;

template <typename T>
void debugIO(T val) {
  T val2;
  std::cout << "-------------------------------------------------------------debugIO(" << typeid(T).name() << ")" << std::endl;
  std::cout << val << std::endl;
  {
    std::ofstream os("dat",std::ios::binary);
    writeBinaryData(val,os);
  }
  {
    std::ifstream is("dat",std::ios::binary);
    readBinaryData(val2,is);
  }
  std::cout << val2 << std::endl;
}
int main(int argc,char** argv) {
  debugIO<int>(1);
  debugIO<char>(1);
  debugIO<bool>(1);
  debugIO<float>(1);
  debugIO<double>(1);
  debugIO<rational>(rational(2,5));
  debugIO<float128>(1);
  debugIO<mpfr_float>(1);
  debugIO<Eigen::Matrix<float128,3,3>>(Eigen::Matrix<float128,3,3>::Random());
  debugIO<Eigen::Matrix<float128,3,-1>>(Eigen::Matrix<float128,3,-1>::Random(3,3));
  debugIO<Eigen::Matrix<float128,-1,3>>(Eigen::Matrix<float128,-1,3>::Random(3,3));
  debugIO<Eigen::Matrix<float128,-1,-1>>(Eigen::Matrix<float128,-1,-1>::Random(3,3));
  return 0;
}
