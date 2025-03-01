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
  GEOMETRY_SCALAR x0=_barrier._x0+_d0;
  _contact->generateManifolds(x0,true,_manifolds,t.template cast<GEOMETRY_SCALAR>());
}
bool MeshBasedPBDSimulator::detectLastContact() {
  //detectContact(_pos._info._TM);
  Vec DE;
  _manifoldsLast.clear();
  T e=normalLastEnergy(_pos,&DE,false);
  if(!isfinite(e))
    return false;
  return !_manifoldsLast.empty();
}
MeshBasedPBDSimulator::T MeshBasedPBDSimulator::normalEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  //std::cout<<_manifolds.size()<<std::endl;
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
      if(!cc.eval(&val,_body.get(),DE?&grad:NULL,&m._DNDX,NULL,NULL,NULL))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else {
        parallelAdd<T>(E,val);
        m._x=cc.getX();
      }
    } else cc.evalbackward(&val,_body.get(),&grad);
  }
  return E;
}
MeshBasedPBDSimulator::T MeshBasedPBDSimulator::normalLastEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  //std::cout<<_manifolds.size()<<std::endl;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifolds.size(); id++) {
    std::vector<ContactManifold> manifoldsLast;
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
      if(!cc.eval(&val,_body.get(),DE?&grad:NULL,&m._DNDX,NULL,NULL,&manifoldsLast))
        parallelAdd<T>(E,std::numeric_limits<T>::infinity());
      else {
        parallelAdd<T>(E,val);
        m._x=cc.getX();
      }
    } else cc.evalbackward(&val,_body.get(),&grad);
    #pragma omp critical
    {
        _manifoldsLast.insert(_manifoldsLast.end(), manifoldsLast.begin(), manifoldsLast.end());
    }
  }
  //std::cout<<_manifoldsLast.size()<<std::endl;
  return E;
}
MeshBasedPBDSimulator::T MeshBasedPBDSimulator::tangentEnergy(GradInfo& grad,Vec* DE,bool backward) {
  T E=0;
  OMP_PARALLEL_FOR_
  for(int id=0; id<(int)_manifoldsLast.size(); id++) {
    if(!isfinite(E))
      continue;
    ContactManifold& m=_manifoldsLast[id];
    GJKPolytope<T> mA=GJKPolytope<T>(m._jidA,std::dynamic_pointer_cast<MeshExact>(m._sA),grad);
    GJKPolytope<T> mB=GJKPolytope<T>(m._jidB,std::dynamic_pointer_cast<MeshExact>(m._sB),grad);
    GJKPolytope<T> mALast=GJKPolytope<T>(m._jidA,std::dynamic_pointer_cast<MeshExact>(m._sA),_pos);
    GJKPolytope<T> mBLast=GJKPolytope<T>(m._jidB,std::dynamic_pointer_cast<MeshExact>(m._sB),_pos);
    if(!mA.mesh() || !mB.mesh())
      continue;
    //find term
    std::vector<CCBarrierMeshFrictionEnergy<T,Barrier>::FrictionTerm>* term=NULL;
    //compute energy/gradient/hessian
    T val=0;
    CCBarrierMeshFrictionEnergy<T,Barrier> cf
    (mA,mB,mALast,mBLast,_barrier,_d0,&grad,_coefBarrier,_dt,term);
    cf.setX(m._x);
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
