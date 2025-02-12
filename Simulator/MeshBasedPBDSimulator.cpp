#include "MeshBasedPBDSimulator.h"
#include <ConvexHull/ConvexHullMeshDistanceEnergy.h>
#include <Environment/ConvexHullExact.h>
#include <Articulated/ArticulatedUtils.h>

namespace PHYSICSMOTION {
MeshBasedPBDSimulator::MeshBasedPBDSimulator(T dt):ConvHullPBDSimulator(dt) {}
void MeshBasedPBDSimulator::step() {
  _termMap.clear();
  ConvHullPBDSimulator::step();
}
void MeshBasedPBDSimulator::detectContact(const Mat3XT& t) {
  _manifolds.clear();
  if(!_contact)
    _contact.reset(new ContactGenerator(_body,_shapes));
  GEOMETRY_SCALAR x0=10000;//_barrier._x0+_d0;
  _contact->generateManifolds(x0,true,_manifolds,t.template cast<GEOMETRY_SCALAR>());
}
MeshBasedPBDSimulator::T MeshBasedPBDSimulator::normalEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifolds.size(); id++) {
    if(!isfinite(E))
      continue;
    ContactManifold& m=_manifolds[id];
    GJKPolytope<T>& mA=m._jidA<0?_obs[m._sidA]:grad._polytopes[m._jidA];
    GJKPolytope<T>& mB=m._jidB<0?_obs[m._sidB]:grad._polytopes[m._jidB];
    if(!mA.mesh() || !mB.mesh())
      continue;
    //compute energy/gradient/hessian
    T val=0;
    CCBarrierMeshEnergy<T,Barrier> cc(mA,mB,_barrier,_d0,&grad,_coefBarrier);
    if(!backward) {
      if(!cc.eval(&val,_body.get(),DE?&grad:NULL,&m._DNDX,NULL,NULL))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else {
        parallelAdd<T>(E,val);
        m._x=cc.getX();
      }
    } else cc.evalbackward(&val,_body.get(),&grad);
  }
  return E;
}
MeshBasedPBDSimulator::T MeshBasedPBDSimulator::tangentEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifoldsLast.size(); id++) {
    if(!isfinite(E))
      continue;
    ContactManifold& m=_manifoldsLast[id];
    GJKPolytope<T>& mA=m._jidA<0?_obs[m._sidA]:grad._polytopes[m._jidA];
    GJKPolytope<T>& mB=m._jidB<0?_obs[m._sidB]:grad._polytopes[m._jidB];
    GJKPolytope<T>& mALast=m._jidA<0?_obs[m._sidA]:_pos._polytopes[m._jidA];
    GJKPolytope<T>& mBLast=m._jidB<0?_obs[m._sidB]:_pos._polytopes[m._jidB];
    if(!mA.mesh() || !mB.mesh())
      continue;
    //find term
    std::vector<CCBarrierMeshFrictionEnergy<T,Barrier>::FrictionTerm>* term=NULL;
    OMP_CRITICAL_ {
      auto it=_termMap.find(m);
      if(it!=_termMap.end())
        term=&(it->second);
    }
    //compute energy/gradient/hessian
    T val=0;
    CCBarrierMeshFrictionEnergy<T,Barrier> cf
    (mA,mB,mALast,mBLast,_barrier,_d0,&grad,_coefBarrier,_dt,term);
    //store term
    OMP_CRITICAL_ {
      auto it=_termMap.find(m);
      if(it==_termMap.end())
        _termMap[m]=cf.terms();
    }
    //evaluate
    //ASSERT(!backward)
    if(!backward) {
      if(!cf.eval(&val,_body.get(),DE?&grad:NULL,NULL,NULL))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else parallelAdd<T>(E,val);
    } else cf.evalbackward(NULL,_body.get(),&grad);
  }
  return E;
}
}
