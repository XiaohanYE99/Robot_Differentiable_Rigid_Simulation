#include "EntityId.h"
#include "ConvexHullDistanceEnergy.h"
#include <Environment/MeshExact.h>

namespace PHYSICSMOTION {
//CollisionGradInfo
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo() {}
template <typename T>
CollisionGradInfo<T>::CollisionGradInfo(const ArticulatedBody& body,const Vec& theta) {
  reset(body,theta);
}
template <typename T>
void CollisionGradInfo<T>::reset(const ArticulatedBody& body,const Vec& theta) {
  _HTheta.setZero(body.nrDOF(),body.nrDOF());
  _DTG.setZero(3,4*body.nrJ());
  _info.reset(body,theta);
  //data
  _globalVss.resize(body.nrJ());
  _polytopes.resize(body.nrJ());
  for(int i=0; i<body.nrJ(); i++) {
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    std::shared_ptr<MeshExact> m;
    if(body.joint(i)._mesh){
      body.joint(i)._mesh->getMesh(vss,iss);
      m.reset(new MeshExact(vss,iss,false));
    }
    //std::shared_ptr<MeshExact> m=std::dynamic_pointer_cast<MeshExact>(body.joint(i)._mesh);
    if(!m) {
      _globalVss[i].resize(3,0);
    } else {
      _globalVss[i].resize(3,m->vss().size());
      for(int j=0; j<(int)m->vss().size(); j++)
        _globalVss[i].col(j)=ROTI(_info._TM,i)*m->vss()[j].template cast<T>()+CTRI(_info._TM,i);
      _polytopes[i].reset(new GJKPolytope<T>(i,body,*this));
    }
  }
}

//instance
template struct CollisionGradInfo<FLOAT>;

}
