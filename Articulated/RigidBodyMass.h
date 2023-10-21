#ifndef RIGID_BODY_MASS_H
#define RIGID_BODY_MASS_H

#include <Eigen/Dense>
#include <Utils/Pragma.h>
#include <vector>

namespace PHYSICSMOTION {
template <typename T>
struct RigidBodyMass {
  DECL_MAT_VEC_MAP_TYPES_T
  RigidBodyMass(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                const std::vector<Eigen::Matrix<int,3,1>>& iss);
  RigidBodyMass(const std::vector<Eigen::Matrix<double,2,1>>& vss,
                const std::vector<Eigen::Matrix<int,3,1>>& iss);
  Vec3T getCtr() const;
  Mat6T getMass() const;
  Mat6T getMassCOM() const;
  T getM() const;
  Vec3T getMC() const;
  Mat3T getMCCT() const;
 private:
  Mat6T _mat,_matCOM;
  Vec3T _ctr;
  Mat3T _MCCT;
};
}

#endif
