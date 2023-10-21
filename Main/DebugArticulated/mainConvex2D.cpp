#include <Environment/EnvironmentUtils.h>
#include <iostream>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  std::vector<Eigen::Matrix<double,3,1>,Eigen::aligned_allocator<Eigen::Matrix<double,3,1>>> vss;
  vss.push_back(Eigen::Matrix<double,3,1>(0,0,0));
  vss.push_back(Eigen::Matrix<double,3,1>(1,0,0));
  vss.push_back(Eigen::Matrix<double,3,1>(0,1,0));
  vss.push_back(Eigen::Matrix<double,3,1>(1,1,0));
  vss.push_back(Eigen::Matrix<double,3,1>(1,1,0));
  vss.push_back(Eigen::Matrix<double,3,1>(2,1,0));
  vss.push_back(Eigen::Matrix<double,3,1>(1,2,0));
  vss.push_back(Eigen::Matrix<double,3,1>(2,2,0));
  makeConvexProject(vss);
  for(auto v:vss)
    std::cout << v.transpose() << std::endl;
  return 0;
}
