#include <Utils/Utils.h>
#include <Eigen/Dense>
#include <random>

using namespace PHYSICSMOTION;

template <typename T>
void debugUtilsScalar() {
  T val=1;
  std::cout << "-------------------------------------------------------------DebugUtilsScalar(" << typeid(T).name() << ")" << std::endl;
  tinyxml2::XMLDocument doc;
  doc.InsertEndChild(doc.NewElement("root"));
  put(doc,"test",val);
  std::cout << get<T>(doc,"test") << std::endl;
  std::cout << get<T>(doc,"testDefault",0) << std::endl;
}
void debugUtilsString() {
  std::cout << "-------------------------------------------------------------DebugUtilsString" << std::endl;
  tinyxml2::XMLDocument doc;
  doc.InsertEndChild(doc.NewElement("root"));
  put<std::string>(doc,"test","test");
  std::cout << get<std::string>(doc,"test") << std::endl;
  std::cout << get<std::string>(doc,"testDefault","testDefault") << std::endl;
}
template <typename T,int r,int c>
void debugUtilsEigen() {
  std::cout << "-------------------------------------------------------------debugUtilsEigen" << std::endl;
  tinyxml2::XMLDocument doc;
  doc.InsertEndChild(doc.NewElement("root"));
  Eigen::Matrix<T,r,c> val;
  val.setRandom();
  putPtree<Eigen::Matrix<T,r,c>>(doc,"test",val);
  std::cout << parsePtree<Eigen::Matrix<T,r,c>>(doc,"test") << std::endl;
  std::cout << parsePtree<Eigen::Matrix<T,-1,-1>>(doc,"test") << std::endl;
}
int main(int argc,char** argv) {
  debugUtilsScalar<int>();
  debugUtilsScalar<bool>();
  debugUtilsScalar<char>();
  debugUtilsScalar<float>();
  debugUtilsScalar<double>();
  debugUtilsScalar<float128>();
  debugUtilsScalar<mpfr_float>();
  debugUtilsString();
  debugUtilsEigen<float,3,4>();
  debugUtilsEigen<double,3,4>();
  debugUtilsEigen<float128,3,4>();
  debugUtilsEigen<mpfr_float,3,4>();
  return 0;
}
