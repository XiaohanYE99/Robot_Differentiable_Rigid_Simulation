#ifndef PB_CENTROID_BODY_DYNAMICS_GRADIENT_INFO_H
#define PB_CENTROID_BODY_DYNAMICS_GRADIENT_INFO_H

#include "ArticulatedBody.h"

namespace PHYSICSMOTION {
template <typename T>
class DeformedEnvironment;
template <typename T>
struct PBCentroidBodyDynamicsGradientInfoMap {
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  PBCentroidBodyDynamicsGradientInfoMap();
  PBCentroidBodyDynamicsGradientInfoMap(const PBCentroidBodyDynamicsGradientInfoMap<T>& other);
  Vec6T C(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
          const PBCentroidBodyDynamicsGradientInfoMap<T>& IN,
          std::function<void(int,int,Vec3T)> DCDINN,
          std::function<void(int,int,Vec3T)> DCDIN,
          std::function<void(int,int,Vec3T)> DCDI,
          std::function<void(int,int,Mat3T)> DCDX,
          std::function<void(int,int,Mat3T)> DCDF) const;
  Vec6T C(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
          const PBCentroidBodyDynamicsGradientInfoMap<T>& IN,
          Mat6TM DCDINN,Mat6TM DCDIN,Mat6TM DCDI,
          Mat6XTM DCDX,Mat6XTM DCDF) const;
  Vec6T C(const PBCentroidBodyDynamicsGradientInfoMap<T>& INN,
          const PBCentroidBodyDynamicsGradientInfoMap<T>& IN) const;
  void setContactForce(int index,const Vec3T& pos,const Vec3T& force);
  void DTG(Mat3X4T GK,std::function<void(int,T)> DTG) const;
  Mat3X4T getTrans() const;
  Vec6T getDOF() const;
 protected:
  Vec3T IXX(const Mat3T& Ra,const Mat3T& Rb,Mat3TM dIdwa,Mat3TM dIdwb) const;
  Mat3X4TM _RTM;
  Mat3XTM _wM,_DTM;
  Mat3XTM _xM,_fM;
  Mat3T _XTX;
  Vec3T _g;
  T _dt,_M;
  bool _isEXP;
};
template <typename T>
struct PBCentroidBodyDynamicsGradientInfo : public PBCentroidBodyDynamicsGradientInfoMap<T> {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_RTM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_wM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_DTM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_xM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_fM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_XTX;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_g;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_dt;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_M;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::_isEXP;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::mapV;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::mapCV;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::mapM;
  using PBCentroidBodyDynamicsGradientInfoMap<T>::mapCM;
  PBCentroidBodyDynamicsGradientInfo();
  PBCentroidBodyDynamicsGradientInfo(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
                                     std::shared_ptr<DeformedEnvironment<T>> DEnv=NULL);
  PBCentroidBodyDynamicsGradientInfo& operator=(const PBCentroidBodyDynamicsGradientInfo& other);
  void setDeformedEnvironment(std::shared_ptr<DeformedEnvironment<T>> DEnv=NULL);
  Vec invReset(const ArticulatedBody& body,const Vec& DOF) const;
  void reset(const ArticulatedBody& body,const Vec& DOF);
  void DEnvReset(const ArticulatedBody& body,const Vec& DOF);
  void init(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
            std::shared_ptr<DeformedEnvironment<T>> DEnv=NULL);
  static void debug(const ArticulatedBody& body,int NContact,const Vec3T& g,T dt,
                    std::shared_ptr<DeformedEnvironment<T>> DEnv=NULL);
  static void debugContinuous(const ArticulatedBody& body,int NContact);
 private:
  std::shared_ptr<DeformedEnvironment<T>> _DEnv;
  void resetPtr();
  Mat3X4T _RT;
  Mat3XT _w,_DT;
  Mat3XT _x,_f;
};
}

#endif
