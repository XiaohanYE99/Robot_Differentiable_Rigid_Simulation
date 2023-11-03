#include "ContactGenerator.h"
#include "CompositeShapeExact.h"
#include "SphericalBBoxExact.h"
#include "GJKK.h"
#include "SAT.h"
#include <stack>

namespace PHYSICSMOTION {
//ContactManifold
ContactGenerator::T ContactGenerator::ContactPoint::depth() const {
  return (_ptA-_ptB).dot(_nA2B);
}
void ContactGenerator::ContactPoint::swap() {
  std::swap(_ptA,_ptB);
  _nA2B*=-1;
}
void ContactGenerator::ContactPoint::transform(const Mat3X4T& t) {
  _ptA=ROT(t)*_ptA+CTR(t);
  _ptB=ROT(t)*_ptB+CTR(t);
  _nA2B=ROT(t)*_nA2B;
}
ContactGenerator::ContactManifold::ContactManifold():_sidA(-1),_sidB(-1),_jidA(-1),_jidB(-1) {}
void ContactGenerator::ContactManifold::swap() {
  for(auto& p:_points)
    p.swap();
  std::swap(_sA,_sB);
  std::swap(_pA,_pB);
  std::swap(_sidA,_sidB);
  std::swap(_jidA,_jidB);
  std::swap(_tA,_tB);
}
//ContactGenerator
ContactGenerator::ContactGenerator(std::shared_ptr<ArticulatedBody> body,std::vector<std::shared_ptr<ShapeExact>> shapes):_body(body),_staticShapes(shapes) {
  Node<int,BBoxExact> n;
  //static BVH
  for(int i=0; i<(int)_staticShapes.size(); i++) {
    n._nrCell=1;
    n._bb=shapes[i]->getBB();
    n._cell=i;
    _staticBVH.push_back(n);
  }
  Node<int,BBoxExact>::buildBVHBottomUpAll(_staticBVH);

  //dynamic BVH
  if(!_body)
    return;
  ArticulatedBody::Mat3XT t=body->getT(ArticulatedBody::Vec::Zero(body->nrDOF()));
  std::unordered_map<int,int> jid2nid;
  for(int i=0; i<_body->nrJ(); i++)
    if(_body->joint(i)._mesh) {
      n._nrCell=1;
      n._cell=i;
      n._bb=_body->joint(i).getBB(TRANSI(t,n._cell));
      jid2nid[i]=(int)_dynamicBVH.size();
      _dynamicBVH.push_back(n);
    }
  //edge map
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  for(int i=0; i<_body->nrJ(); i++)
    if(_body->joint(i)._mesh) {
      int p=_body->joint(i)._parent;
      while(p>=0 && !_body->joint(p)._mesh)
        p=_body->joint(p)._parent;
      if(p>=0) {
        //ith joint is linked to jth joint, both of which have associated meshes
        edgeMap[Eigen::Matrix<int,2,1>(i,p)]=std::make_pair(jid2nid[i],jid2nid[p]);
        edgeMap[Eigen::Matrix<int,2,1>(p,i)]=std::make_pair(jid2nid[p],jid2nid[i]);
        _exclude.insert(Eigen::Matrix<int,2,1>(i,p));
        _exclude.insert(Eigen::Matrix<int,2,1>(p,i));
      }
    }
  Node<int,BBoxExact>::buildBVHBottomUp(_dynamicBVH,edgeMap,true);
}
void ContactGenerator::generateManifolds(std::vector<ContactManifold>& manifolds,Mat3XT t,T x0,int status, bool usingCCD) {
  updateBVH(t,x0);
  std::stack<std::pair<int,int>> ss;
  //static-static
  if((status&STATIC_STATIC) && !_staticBVH.empty()) {
    ss.push(std::make_pair((int)_staticBVH.size()-1,(int)_staticBVH.size()-1));
    while(!ss.empty()) {
      int ia=ss.top().first;
      int ib=ss.top().second;
      const auto& na=_staticBVH[ia];
      const auto& nb=_staticBVH[ib];
      ss.pop();
      if(!na._bb.intersect(nb._bb))
        continue;
      else if(na._cell>=0 && nb._cell>=0) {
        if(na._cell>=nb._cell)
          continue;
        ContactManifold m;
        m._sidA=na._cell;
        m._sidB=nb._cell;
        m._sA=_staticShapes[m._sidA];
        m._sB=_staticShapes[m._sidB];
        m._tA.setIdentity();
        m._tB.setIdentity();
        if(! usingCCD) generateManifold(manifolds,m);
        else manifolds.push_back(m);
      } else if(na._cell>=0) {
        ss.push(std::make_pair(ia,nb._l));
        ss.push(std::make_pair(ia,nb._r));
      } else if(nb._cell>=0) {
        ss.push(std::make_pair(na._l,ib));
        ss.push(std::make_pair(na._r,ib));
      } else {
        ss.push(std::make_pair(na._l,nb._l));
        ss.push(std::make_pair(na._l,nb._r));
        ss.push(std::make_pair(na._r,nb._l));
        ss.push(std::make_pair(na._r,nb._r));
      }
    }
  }
  //static-dynamic
  if((status&STATIC_DYNAMIC) && !_staticBVH.empty() && !_dynamicBVH.empty()) {
    ss.push(std::make_pair((int)_staticBVH.size()-1,(int)_dynamicBVH.size()-1));
    while(!ss.empty()) {
      int ia=ss.top().first;
      int ib=ss.top().second;
      const auto& na=_staticBVH[ia];
      const auto& nb=_dynamicBVH[ib];
      ss.pop();
      if(!na._bb.intersect(nb._bb))
        continue;
      else if(na._cell>=0 && nb._cell>=0) {
        ContactManifold m;
        m._sidA=na._cell;
        m._jidB=nb._cell;
        m._sA=_staticShapes[m._sidA];
        m._sB=_body->joint(m._jidB)._mesh;
        m._tA.setIdentity();
        m._tB=TRANSI(t,m._jidB);
        if(! usingCCD) generateManifold(manifolds,m);
        else manifolds.push_back(m);
      } else if(na._cell>=0) {
        ss.push(std::make_pair(ia,nb._l));
        ss.push(std::make_pair(ia,nb._r));
      } else if(nb._cell>=0) {
        ss.push(std::make_pair(na._l,ib));
        ss.push(std::make_pair(na._r,ib));
      } else {
        ss.push(std::make_pair(na._l,nb._l));
        ss.push(std::make_pair(na._l,nb._r));
        ss.push(std::make_pair(na._r,nb._l));
        ss.push(std::make_pair(na._r,nb._r));
      }
    }
  }
  //dynamic-dynamic
  if((status&DYNAMIC_DYNAMIC) && !_dynamicBVH.empty()) {
    ss.push(std::make_pair((int)_dynamicBVH.size()-1,(int)_dynamicBVH.size()-1));
    while(!ss.empty()) {
      int ia=ss.top().first;
      int ib=ss.top().second;
      const auto& na=_dynamicBVH[ia];
      const auto& nb=_dynamicBVH[ib];
      ss.pop();
      if(!na._bb.intersect(nb._bb))
        continue;
      else if(na._cell>=0 && nb._cell>=0) {
        if(na._cell>=nb._cell)
          continue;
        if(_exclude.find(Eigen::Matrix<int,2,1>(na._cell,nb._cell))!=_exclude.end())
          continue;
        ContactManifold m;
        m._jidA=na._cell;
        m._jidB=nb._cell;
        m._sA=_body->joint(m._jidA)._mesh;
        m._sB=_body->joint(m._jidB)._mesh;
        m._tA=TRANSI(t,m._jidA);
        m._tB=TRANSI(t,m._jidB);
        if(! usingCCD) generateManifold(manifolds,m);
        else manifolds.push_back(m);
      } else if(na._cell>=0) {
        ss.push(std::make_pair(ia,nb._l));
        ss.push(std::make_pair(ia,nb._r));
      } else if(nb._cell>=0) {
        ss.push(std::make_pair(na._l,ib));
        ss.push(std::make_pair(na._r,ib));
      } else {
        ss.push(std::make_pair(na._l,nb._l));
        ss.push(std::make_pair(na._l,nb._r));
        ss.push(std::make_pair(na._r,nb._l));
        ss.push(std::make_pair(na._r,nb._r));
      }
    }
  }
}
void ContactGenerator::updateBVH(Mat3XT& t, T x0) {
  //static BVH
  for(int i=0; i<(int)_staticBVH.size(); i++) {
    auto& n=_staticBVH[i];
    if(n._cell>=0)
      n._bb=_staticShapes[n._cell]->getBB();
    else {
      n._bb=BBoxExact();
      n._bb.setUnion(_staticBVH[n._l]._bb);
      n._bb.setUnion(_staticBVH[n._r]._bb);
    }
    n._bb.extendUnion(x0);
  }
  //dynamic BVH
  for(int i=0; i<(int)_dynamicBVH.size(); i++) {
    auto& n=_dynamicBVH[i];
    if(n._cell>=0) {
      ArticulatedBody::Mat3X4T tGOut;
      n._bb=_body->joint(n._cell).getBB(TRANSI(t,n._cell).template cast<ArticulatedBody::T>(),&tGOut);
      TRANSI(t,n._cell)=tGOut.template cast<T>();
    } else {
      n._bb=BBoxExact();
      n._bb.setUnion(_dynamicBVH[n._l]._bb);
      n._bb.setUnion(_dynamicBVH[n._r]._bb);
    }
    n._bb.extendUnion(x0);
  }
}
BBoxExact ContactGenerator::getBB() const {
  if(_staticBVH.empty() && _dynamicBVH.empty())
    return BBoxExact();
  else if(_staticBVH.empty())
    return _dynamicBVH.back()._bb;
  else if(_dynamicBVH.empty())
    return _staticBVH.back()._bb;
  else {
    BBoxExact bb=_staticBVH.back()._bb;
    bb.setUnion(_dynamicBVH.back()._bb);
    return bb;
  }
}
ContactGenerator::T ContactGenerator::epsDist() {
  return _epsDist;
}
ContactGenerator::T ContactGenerator::epsDir() {
  return _epsDir;
}
//helper
void ContactGenerator::generateManifold(std::vector<ContactManifold>& manifolds,ContactManifold m) {
  auto cA=std::dynamic_pointer_cast<CompositeShapeExact>(m._sA);
  auto cB=std::dynamic_pointer_cast<CompositeShapeExact>(m._sB);
  if(cA && cB) {
    for(int i=0; i<(int)cA->getGeoms().size(); i++)
      for(int j=0; j<(int)cB->getGeoms().size(); j++) {
        ContactManifold m2=m;
        m2._sA=cA->getGeoms()[i];
        m2._sB=cB->getGeoms()[j];
        APPLY_TRANS(m2._tA,cA->getTrans()[i].template cast<T>(),m._tA);
        APPLY_TRANS(m2._tB,cB->getTrans()[j].template cast<T>(),m._tB);
        generateManifold(manifolds,m2);
      }
  } else if(cA) {
    for(int i=0; i<(int)cA->getGeoms().size(); i++) {
      ContactManifold m2=m;
      m2._sA=cA->getGeoms()[i];
      APPLY_TRANS(m2._tA,cA->getTrans()[i].template cast<T>(),m._tA);
      generateManifold(manifolds,m2);
    }
  } else if(cB) {
    for(int j=0; j<(int)cB->getGeoms().size(); j++) {
      ContactManifold m2=m;
      m2._sB=cB->getGeoms()[j];
      APPLY_TRANS(m2._tB,cB->getTrans()[j].template cast<T>(),m._tB);
      generateManifold(manifolds,m2);
    }
  } else {
    //case-by-case generation
    //sphere-sphere
    if(generateManifoldSphereSphere(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    //sphere-capsule
    if(generateManifoldSphereCapsule(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    m.swap();
    //capsule-sphere
    if(generateManifoldSphereCapsule(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    //capsule-capsule
    if(generateManifoldCapsuleCapsule(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    //sphere-box
    if(generateManifoldSphereBox(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    m.swap();
    //box-sphere
    if(generateManifoldSphereBox(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    //capsule-box
    if(generateManifoldCapsuleBox(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    m.swap();
    //box-capsule
    if(generateManifoldCapsuleBox(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    //box-box
    if(generateManifoldBoxBox(manifolds,m)) {
      if(!m._points.empty())
        manifolds.push_back(m);
      return;
    }
    ASSERT_MSG(false,"Not supported contact case!")
  }
}
void ContactGenerator::generateManifoldSphereSphereInternal(std::vector<ContactManifold>& manifolds,ContactManifold& m,const Vec3T& cA1,const Vec3T& cA2,const Vec3T& cB1,const Vec3T& cB2) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sB);
  Vec4T distSqrs((cA1-cB1).squaredNorm(),(cA2-cB1).squaredNorm(),(cA1-cB2).squaredNorm(),(cA2-cB2).squaredNorm());
  Vec4T::Index id;
  T distSqr=distSqrs.minCoeff(&id),dist=0;
  T sumRad=sA->radius()+sB->radius(),sumRadSqr=sumRad*sumRad;
  const Vec3T& cA=(id%2==0)?cA1:cA2;
  const Vec3T& cB=(id<2)?cB1:cB2;
  //not in contact
  if(distSqr>sumRadSqr)
    return;
  //in contact
  if(distSqr>_epsDist*_epsDist) {
    //a single contact point
    p._nA2B=cB-cA;
    p._nA2B/=dist=sqrt((double)distSqr);
  } else {
    //overlapping degenerate case
    distSqrs.maxCoeff(&id);
    const Vec3T& cAn=(id%2==0)?cA1:cA2;
    const Vec3T& cBn=(id<2)?cB1:cB2;
    p._nA2B=cBn-cAn;
    p._nA2B/=p._nA2B.template cast<double>().norm();
  }
  p._ptA=cA+p._nA2B*sA->radius();
  p._ptB=cB-p._nA2B*sB->radius();
  m._points.push_back(p);
}
bool ContactGenerator::generateManifoldSphereSphere(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sB);
  if(!sA || !sB)
    return false;
  else if(sA->isSphere() && sB->isSphere()) {
    Vec3T cA=ROT(m._tA)*sA->minCorner().template cast<T>()+CTR(m._tA);
    Vec3T cB=ROT(m._tB)*sB->minCorner().template cast<T>()+CTR(m._tB);
    T distSqr=(cA-cB).squaredNorm(),dist=0;
    T sumRad=sA->radius()+sB->radius(),sumRadSqr=sumRad*sumRad;
    //not in contact
    if(distSqr>sumRadSqr)
      return true;
    //in contact
    if(distSqr>_epsDist*_epsDist) {
      //a single contact point
      p._nA2B=cB-cA;
      p._nA2B/=dist=sqrt((double)distSqr);
    } else {
      //overlapping degenerate case
      p._nA2B=Vec3T(0,0,1);
    }
    p._ptA=cA+p._nA2B*sA->radius();
    p._ptB=cB-p._nA2B*sB->radius();
    m._points.push_back(p);
    return true;
  } else return false;
}
void ContactGenerator::generateManifoldSphereCapsuleInternal(std::vector<ContactManifold>& manifolds,ContactManifold& m,const Vec3T& cA,const Vec3T& cB1,const Vec3T& cB2) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sB);
  Vec3T n=cB2-cB1;
  T nLenSqr=n.squaredNorm(),nLen=sqrt((double)nLenSqr);
  n/=nLen;
  T d=(cA-cB1).dot(n);
  T sumRad=sA->radius()+sB->radius(),sumRadSqr=sumRad*sumRad;
  //three cases
  if(d<=0) {
    T distSqr=(cA-cB1).squaredNorm(),dist=0;
    //not in contact
    if(distSqr>sumRadSqr)
      return;
    //in contact
    if(distSqr>_epsDist*_epsDist) {
      //a single contact point
      p._nA2B=cB1-cA;
      p._nA2B/=dist=sqrt((double)distSqr);
    } else {
      p._nA2B=n;
    }
    p._ptA=cA+p._nA2B*sA->radius();
    p._ptB=cB1-p._nA2B*sB->radius();
  } else if(d>=nLen) {
    T distSqr=(cA-cB2).squaredNorm(),dist=0;
    //not in contact
    if(distSqr>sumRadSqr)
      return;
    //in contact
    if(distSqr>_epsDist*_epsDist) {
      //a single contact point
      p._nA2B=cB2-cA;
      p._nA2B/=dist=sqrt((double)distSqr);
    } else {
      p._nA2B=-n;
    }
    p._ptA=cA+p._nA2B*sA->radius();
    p._ptB=cB2-p._nA2B*sB->radius();
  } else if(d>0 && d<nLen) {
    Vec3T dir=cA-cB1-n*d;
    T distSqr=dir.squaredNorm(),dist=0;
    //not in contact
    if(distSqr>sumRadSqr)
      return;
    //in contact
    if(distSqr>_epsDist*_epsDist) {
      p._nA2B=-dir;
      p._nA2B/=dist=sqrt((double)distSqr);
    } else {
      Vec3T::Index id;
      n.cwiseAbs().minCoeff(&id);
      p._nA2B=n.cross(Vec3T::Unit(id));
      p._nA2B/=p._nA2B.template cast<double>().norm();
    }
    p._ptA=cA+p._nA2B*sA->radius();
    p._ptB=cA-dir-p._nA2B*sB->radius();
  }
  m._points.push_back(p);
}
bool ContactGenerator::generateManifoldSphereCapsule(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sB);
  if(!sA || !sB)
    return false;
  else if(sA->isSphere() && sB->isCapsule()) {
    Vec3T cA=ROT(m._tA)*sA->minCorner().template cast<T>()+CTR(m._tA);
    Vec3T cB1=ROT(m._tB)*sB->minCorner().template cast<T>()+CTR(m._tB);
    Vec3T cB2=ROT(m._tB)*sB->maxCorner().template cast<T>()+CTR(m._tB);
    generateManifoldSphereCapsuleInternal(manifolds,m,cA,cB1,cB2);
    return true;
  } else return false;
}
bool ContactGenerator::generateManifoldCapsuleCapsule(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sB);
  if(!sA || !sB)
    return false;
  else if(sA->isCapsule() && sB->isCapsule()) {
    Vec3T cA1=ROT(m._tA)*sA->minCorner().template cast<T>()+CTR(m._tA);
    Vec3T cA2=ROT(m._tA)*sA->maxCorner().template cast<T>()+CTR(m._tA);
    Vec3T cB1=ROT(m._tB)*sB->minCorner().template cast<T>()+CTR(m._tB);
    Vec3T cB2=ROT(m._tB)*sB->maxCorner().template cast<T>()+CTR(m._tB);
    Vec3T nA=cA2-cA1,nB=cB2-cB1;
    T nLenASqr=nA.squaredNorm(),nLenA=sqrt((double)nLenASqr);
    T nLenBSqr=nB.squaredNorm(),nLenB=sqrt((double)nLenBSqr);
    nA/=nLenA;
    nB/=nLenB;
    if(abs(nA.dot(nB))>1-_epsDir) {
      //nearly parallel
      T dB1=(cB1-cA1).dot(nA);
      T dB2=(cB2-cA1).dot(nA);
      if(dB1<=0 && dB2<=0) {
        //sphere-sphere
        generateManifoldSphereSphereInternal(manifolds,m,cA1,cA2,cB1,cB2);
      } else if(dB1>=nLenA && dB2>=nLenA) {
        //sphere-sphere
        generateManifoldSphereSphereInternal(manifolds,m,cA1,cA2,cB1,cB2);
      } else {
        //range
        Vec3T dir=cB1-cA1-dB1*nA;
        T distSqr=dir.squaredNorm(),dist=0;
        T sumRad=sA->radius()+sB->radius(),sumRadSqr=sumRad*sumRad;
        //not in contact
        if(distSqr>sumRadSqr)
          return true;
        //in contact
        if(distSqr>_epsDist*_epsDist) {
          p._nA2B=dir;
          p._nA2B/=dist=sqrt((double)distSqr);
        } else {
          Vec3T::Index id;
          nA.cwiseAbs().minCoeff(&id);
          p._nA2B=nA.cross(Vec3T::Unit(id));
          p._nA2B/=p._nA2B.template cast<double>().norm();
        }
        //two contacts
        Vec2T range(std::max<T>(0,std::min(dB1,dB2)),std::min(nLenA,std::max(dB1,dB2)));
        for(const T& r:range) {
          p._ptA=cA1+nA*r+p._nA2B*sA->radius();
          p._ptB=cA1+nA*r+dir-p._nA2B*sB->radius();
          m._points.push_back(p);
        }
      }
    } else {
      //not parallel
      Mat3X2T LHS;
      LHS.col(0)=-(cA2-cA1);
      LHS.col(1)= (cB2-cB1);
      Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(cA1-cB1));
      if((bary.array()>=0).all() && (bary.array()<=1).all()) {
        Vec3T cA=cA1*(1-bary[0])+cA2*bary[0];
        Vec3T cB=cB1*(1-bary[1])+cB2*bary[1];
        T distSqr=(cA-cB).squaredNorm();
        T sumRad=sA->radius()+sB->radius(),sumRadSqr=sumRad*sumRad;
        if(distSqr>sumRadSqr)
          return true;
        p._nA2B=nA.cross(nB);
        p._nA2B/=p._nA2B.template cast<double>().norm();
        if((cB-cA).dot(p._nA2B)<0)
          p._nA2B*=-1;
        p._ptA=cA+p._nA2B*sA->radius();
        p._ptB=cB-p._nA2B*sB->radius();
        m._points.push_back(p);
      } else {
        int nC=(int)m._points.size();
        generateManifoldSphereCapsuleInternal(manifolds,m,cA1,cB1,cB2);
        generateManifoldSphereCapsuleInternal(manifolds,m,cA2,cB1,cB2);
        m.swap();
        generateManifoldSphereCapsuleInternal(manifolds,m,cB1,cA1,cA2);
        generateManifoldSphereCapsuleInternal(manifolds,m,cB2,cA1,cA2);
        if((int)m._points.size()==nC)
          return true;  //no new points added
        for(int i=nC+1; i<(int)m._points.size(); i++)
          if(m._points[i].depth()>m._points[nC].depth())
            m._points[nC]=m._points[i];
        m._points.resize(nC+1);
      }
    }
    return true;
  } else return false;
}
bool ContactGenerator::generateManifoldSphereBox(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<BBoxExact>(m._sB);
  if(sA && sA->isSphere() && sB) {
    BBoxExact sBL=sB->enlarged(Vec3T::Constant(_epsDist));  //for numerical safety
    Vec3T cA=ROT(m._tA)*sA->minCorner().template cast<T>()+CTR(m._tA);
    Vec3T cAB=ROT(m._tB).transpose()*(cA-CTR(m._tB));
    if(sBL.contain(cAB)) {
      T minDist=std::numeric_limits<double>::max(),dist=0;
      for(int d=0; d<3; d++) {
        dist=cAB[d]-(T)sBL.minCorner()[d];
        if(dist<=minDist) {
          minDist=dist;
          p._nA2B=Vec3T::Unit(d);
        }
        dist=(T)sBL.maxCorner()[d]-cAB[d];
        if(dist<=minDist) {
          minDist=dist;
          p._nA2B=-Vec3T::Unit(d);
        }
      }
      p._ptA=cAB+p._nA2B*(T)sA->radius();
      p._ptB=cAB-p._nA2B*minDist;
      p.transform(m._tB);
      m._points.push_back(p);
    } else {
      Vec3T dist=Vec3T::Zero();
      for(int i=0; i<3; i++) {
        if (cAB[i]<sB->minCorner()[i])
          dist[i]=cAB[i]-(T)sB->minCorner()[i];
        else if(cAB[i]>sB->maxCorner()[i])
          dist[i]=cAB[i]-(T)sB->maxCorner()[i];
      }
      T distSqr=dist.squaredNorm();
      if(distSqr>=sA->radius()*sA->radius())
        return true;
      p._nA2B=-dist/sqrt((double)distSqr);
      p._ptA=cAB+p._nA2B*sA->radius();
      p._ptB=cAB-dist;
      p.transform(m._tB);
      m._points.push_back(p);
    }
    return true;
  } else return false;
}
bool ContactGenerator::generateManifoldCapsuleBox(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  ContactPoint p;
  auto sA=std::dynamic_pointer_cast<SphericalBBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<BBoxExact>(m._sB);
  if(sA && sA->isCapsule() && sB) {
    Vec3T pAL,pBL;
    bool intersect;
    OMP_CRITICAL_ {
      //facets
      if(_facetCache.find(sA)==_facetCache.end())
        _facetCache[sA]=sA->facets();
      if(_facetCache.find(sB)==_facetCache.end())
        _facetCache[sB]=sB->facets();
      //edges
      if(_edgeCache.find(sA)==_edgeCache.end())
        _edgeCache[sA]=sA->edges();
      if(_edgeCache.find(sB)==_edgeCache.end())
        _edgeCache[sB]=sB->edges();
    }
    auto FA=_facetCache.find(sA)->second;
    auto FB=_facetCache.find(sB)->second;
    auto EA=_edgeCache.find(sA)->second;
    auto EB=_edgeCache.find(sB)->second;
    T distSqr=GJKK::runGJK(sA,sB,m._tA,m._tB,pAL,pBL,&intersect);
    if(intersect || distSqr<_epsDist*_epsDist) {
      SAT::generateManifold(sA,sB,FA,FB,EA,EB,m._tA,m._tB,m);
      for(auto& p:m._points) {
        p._ptA+=p._nA2B*sA->radius();
      }
    } else if(distSqr<sA->radius()*sA->radius()) {
      Facet fA;
      Vec3T cA1=ROT(m._tA)*sA->minCorner().template cast<T>()+CTR(m._tA);
      Vec3T cA2=ROT(m._tA)*sA->maxCorner().template cast<T>()+CTR(m._tA);
      fA._boundary.push_back(ROT(m._tB).transpose()*(cA1-CTR(m._tB)));
      fA._boundary.push_back(ROT(m._tB).transpose()*(cA2-CTR(m._tB)));
      for(const auto& fB:FB)
        if(abs(fB._n.dot(fA._boundary[0]-fA._boundary[1]))<_epsDir && (fA._boundary[0]-fB._boundary[0]).dot(fB._n)>0) {
          //we can return multiple contacts
          SAT::clip(fA,fB);
          for(const auto& pA:fA._boundary) {
            p._ptA=ROT(m._tB)*(pA-fB._n*sA->radius())+CTR(m._tB);
            p._ptB=ROT(m._tB)*(pA-fB._n*(pA-fB._boundary[0]).dot(fB._n))+CTR(m._tB);
            p._nA2B=-ROT(m._tB)*fB._n;
            m._points.push_back(p);
          }
          return true;
        }
      //just return one closest point
      if(distSqr>_epsDist*_epsDist) {
        p._ptA=ROT(m._tA)*pAL+CTR(m._tA);
        p._ptB=ROT(m._tB)*pBL+CTR(m._tB);
        p._nA2B=(p._ptB-p._ptA)/sqrt((double)distSqr);
        p._ptA+=p._nA2B*sA->radius();
        m._points.push_back(p);
      }
      return true;
    }
    return true;
  } else return false;
}
bool ContactGenerator::generateManifoldBoxBox(std::vector<ContactManifold>& manifolds,ContactManifold& m) {
  auto sA=std::dynamic_pointer_cast<BBoxExact>(m._sA);
  auto sB=std::dynamic_pointer_cast<BBoxExact>(m._sB);
  if(sA && sB) {
    OMP_CRITICAL_ {
      //facets
      if(_facetCache.find(sA)==_facetCache.end())
        _facetCache[sA]=sA->facets();
      if(_facetCache.find(sB)==_facetCache.end())
        _facetCache[sB]=sB->facets();
      //edges
      if(_edgeCache.find(sA)==_edgeCache.end())
        _edgeCache[sA]=sA->edges();
      if(_edgeCache.find(sB)==_edgeCache.end())
        _edgeCache[sB]=sB->edges();
    }
    auto FA=_facetCache.find(sA)->second;
    auto FB=_facetCache.find(sB)->second;
    auto EA=_edgeCache.find(sA)->second;
    auto EB=_edgeCache.find(sB)->second;
    SAT::generateManifold(sA,sB,FA,FB,EA,EB,m._tA,m._tB,m);
    return true;
  } else return false;
}
ContactGenerator::T ContactGenerator::_epsDir=1e-3;
ContactGenerator::T ContactGenerator::_epsDist=1e-3;
}
