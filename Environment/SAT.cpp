#include "SAT.h"
#include <Utils/Utils.h>
#include <Utils/Interp.h>

namespace PHYSICSMOTION {
SAT::T SAT::depth(const ProjRange& rng) {
  return std::min(rng._rngA[1]-rng._rngB[0],rng._rngB[1]-rng._rngA[0]);
}
bool SAT::matchWitness(const ProjRange& rng,T witnessA,T*) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if(abs(witnessA-rng._rngA[1])<ContactGenerator::epsDist())
      return true;
  } else {
    //check matching witness
    if(abs(witnessA-rng._rngA[0])<ContactGenerator::epsDist())
      return true;
  }
  //witness does not match
  return false;
}
bool SAT::matchWitness(const ProjRange& rng,T*,T witnessB) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if(abs(witnessB-rng._rngB[0])<ContactGenerator::epsDist())
      return true;
  } else {
    //check matching witness
    if(abs(witnessB-rng._rngB[1])<ContactGenerator::epsDist())
      return true;
  }
  //witness does not match
  return false;
}
bool SAT::matchWitness(const ProjRange& rng,T witnessA,T witnessB) {
  if(rng._rngA[1]-rng._rngB[0]<rng._rngB[1]-rng._rngA[0]) {
    //check matching witness
    if( abs(witnessA-rng._rngA[1])<ContactGenerator::epsDist() &&
        abs(witnessB-rng._rngB[0])<ContactGenerator::epsDist())
      return true;
  } else {
    //check matching witness
    if( abs(witnessA-rng._rngA[0])<ContactGenerator::epsDist() &&
        abs(witnessB-rng._rngB[1])<ContactGenerator::epsDist())
      return true;
  }
  //witness does not match
  return false;
}
SAT::Vec2T SAT::project(std::shared_ptr<ShapeExact> S,
                        const Mat3X4T& trans,
                        const Vec3T& n) {
  Vec2T ret=S->project(ROT(trans).transpose()*n);
  ret=(ret.array()+n.dot(CTR(trans))).matrix();
  sort2(ret[0],ret[1]);
  return ret;
}
SAT::ProjRange SAT::runSAT(std::shared_ptr<ShapeExact> A,
                           std::shared_ptr<ShapeExact> B,
                           const std::vector<Facet>& FA,
                           const std::vector<Facet>& FB,
                           const std::vector<Edge>& EA,
                           const std::vector<Edge>& EB,
                           const Mat3X4T& transA,
                           const Mat3X4T& transB,bool* intersect) {
  ProjRange minRng,rng;
  minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
  minRng._depth=std::numeric_limits<double>::max();
  if(intersect)
    *intersect=false;
  //faceA
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(int i=0; i<(int)FA.size(); i++) {
    rng._fidA=i;
    rng._n=ROT(transA)*FA[i]._n;
    rng._n=rng._n.template cast<double>().normalized().template cast<T>();
    rng._rngA=project(A,transA,rng._n);
    rng._rngB=project(B,transB,rng._n);
    rng._depth=depth(rng);
    if(rng._depth<=0) {
      if(intersect)
        *intersect=false;
      minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
      return minRng;
    } else if(rng._depth<minRng._depth) {
      Vec3T cA=ROT(transA)*FA[i]._boundary[0]+CTR(transA);
      if(!matchWitness(rng,cA.dot(rng._n),(T*)NULL))
        continue;
      if(intersect)
        *intersect=true;
      minRng=rng;
    }
  }
  //faceB
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(int i=0; i<(int)FB.size(); i++) {
    rng._fidB=i;
    rng._n=ROT(transB)*FB[i]._n;
    rng._n=rng._n.template cast<double>().normalized().template cast<T>();
    rng._rngA=project(A,transA,rng._n);
    rng._rngB=project(B,transB,rng._n);
    rng._depth=depth(rng);
    if(rng._depth<=0) {
      if(intersect)
        *intersect=false;
      minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
      return minRng;
    } else if(rng._depth<minRng._depth) {
      Vec3T cB=ROT(transB)*FB[i]._boundary[0]+CTR(transB);
      if(!matchWitness(rng,(T*)NULL,cB.dot(rng._n)))
        continue;
      if(intersect)
        *intersect=true;
      minRng=rng;
    }
  }
  //edgeAB
  rng._fidA=rng._eidA=rng._fidB=rng._eidB=-1;
  for(rng._eidA=0; rng._eidA<(int)EA.size(); rng._eidA++) {
    Vec3T cA1=ROT(transA)*EA[rng._eidA]._a+CTR(transA);
    Vec3T cA2=ROT(transA)*EA[rng._eidA]._b+CTR(transA);
    Vec3T dA=ROT(transA)*(EA[rng._eidA]._b-EA[rng._eidA]._a);
    for(rng._eidB=0; rng._eidB<(int)EB.size(); rng._eidB++) {
      rng._n=dA.cross(ROT(transB)*(EB[rng._eidB]._b-EB[rng._eidB]._a));
      T nLenSqr=rng._n.squaredNorm();
      if(nLenSqr<Epsilon<double>::defaultEps())
        continue;
      //project range with witness, to ensure the closest point pair lies on the line segment
      rng._n/=sqrt((double)nLenSqr);
      rng._rngA=project(A,transA,rng._n);
      rng._rngB=project(B,transB,rng._n);
      rng._depth=depth(rng);
      if(rng._depth<=0) {
        if(intersect)
          *intersect=false;
        minRng._fidA=minRng._eidA=minRng._fidB=minRng._eidB=-1;
        return minRng;
      } else if(rng._depth<minRng._depth-ContactGenerator::epsDist()) { //we prefer face intersection
        //edge-edge case, only a single contact is generated
        Vec3T cB1=ROT(transB)*EB[rng._eidB]._a+CTR(transB);
        Vec3T cB2=ROT(transB)*EB[rng._eidB]._b+CTR(transB);
        //find intersection coordinates
        Mat3X2T LHS;
        LHS.col(0)=-(cA2-cA1);
        LHS.col(1)= (cB2-cB1);
        Vec2T bary=(LHS.transpose()*LHS).inverse()*(LHS.transpose()*(cA1-cB1));
        if((bary.array()<0).any() || (bary.array()>1).any())
          continue;
        if(!matchWitness(rng,cA1.dot(rng._n),cB1.dot(rng._n)))
          continue;
        if(intersect)
          *intersect=true;
        //only (projected) intersecting line segments form potential closest point pair
        rng._ptA=cA1*(1-bary[0])+cA2*bary[0];
        rng._ptB=cB1*(1-bary[1])+cB2*bary[1];
        minRng=rng;
      }
    }
  }
  return minRng;
}
SAT::ProjRange SAT::runSAT(std::shared_ptr<ShapeExact> A,
                           std::shared_ptr<ShapeExact> B,
                           const Mat3X4T& transA,
                           const Mat3X4T& transB,bool* intersect) {
  return runSAT(A,B,A->facets(),B->facets(),A->edges(),B->edges(),transA,transB,intersect);
}
//not only run contact, but also generate manifold
SAT::ProjRange SAT::generateManifold(std::shared_ptr<ShapeExact> A,
                                     std::shared_ptr<ShapeExact> B,
                                     const std::vector<Facet>& FA,
                                     const std::vector<Facet>& FB,
                                     const std::vector<Edge>& EA,
                                     const std::vector<Edge>& EB,
                                     const Mat3X4T& transA,
                                     const Mat3X4T& transB,
                                     ContactManifold& m,bool* Intersect) {
  bool intersect;
  ProjRange rng=runSAT(A,B,FA,FB,EA,EB,transA,transB,&intersect);
  if(Intersect)
    *Intersect=intersect;
  if(!intersect)
    return rng;
  if(rng._eidA>=0 && rng._eidB>=0) {
    //edge-edge case, only a single contact is generated
    ContactPoint p;
    p._ptA=rng._ptA;
    p._ptB=rng._ptB;
    p._nA2B=rng._n;
    if(p.depth()<0)
      p._nA2B*=-1;
    m._points.push_back(p);
  } else if(rng._fidA>=0) {
    generateManifoldFace(A,B,FA,FB,EA,EB,transA,transB,rng._fidA,m._points);
  } else {
    int nrP=(int)m._points.size();
    generateManifoldFace(B,A,FB,FA,EB,EA,transB,transA,rng._fidB,m._points);
    for(int i=nrP; i<(int)m._points.size(); i++)
      m._points[i].swap();
  }
  return rng;
}
SAT::ProjRange SAT::generateManifold(std::shared_ptr<ShapeExact> A,
                                     std::shared_ptr<ShapeExact> B,
                                     const Mat3X4T& transA,
                                     const Mat3X4T& transB,
                                     ContactManifold& m,bool* Intersect) {
  return generateManifold(A,B,A->facets(),B->facets(),A->edges(),B->edges(),transA,transB,m,Intersect);
}
void SAT::generateManifoldFace(std::shared_ptr<ShapeExact> A,
                               std::shared_ptr<ShapeExact> B,
                               const std::vector<Facet>& FA,
                               const std::vector<Facet>& FB,
                               const std::vector<Edge>& EA,
                               const std::vector<Edge>& EB,
                               const Mat3X4T& transA,
                               const Mat3X4T& transB,
                               int fidA,std::vector<ContactPoint>& points) {
  ASSERT(fidA>=0)
  const Vec3T& n=FA[fidA]._n;
  //clip B against A
  ContactPoint p;
  Facet b;
  if(FB.empty()) {
    //this is a capsule
    ASSERT((int)EB.size()==1)
    b._boundary.push_back(EB[0]._a);
    b._boundary.push_back(EB[0]._b);
  } else {
    //find most negative facet
    int minId=-1;
    T minVal=1,val=0;
    for(int i=0; i<(int)FB.size(); i++) {
      val=(ROT(transB)*FB[i]._n).dot(ROT(transA)*n);
      if(val<minVal) {
        minVal=val;
        minId=i;
      }
    }
    //the facet must be negative
    if(minVal>=0)
      return;
    b=FB[minId];
  }

  //transform B->A and clip
  for(auto& v:b._boundary) {
    v=ROT(transB)*v+CTR(transB);
    v=ROT(transA).transpose()*(v-CTR(transA));
  }
  clip(b,FA[fidA]);

  //retain negative vertices
  for(const auto& v:b._boundary) {
    T d=(v-FA[fidA]._boundary[0]).dot(n);
    if(d<0) {
      p._ptA=ROT(transA)*(v-d*n)+CTR(transA);
      p._ptB=ROT(transA)*v+CTR(transA);
      p._nA2B=ROT(transA)*n;
      points.push_back(p);
    }
  }
}
void SAT::clip(Facet& f,const Vec3T& pos,const Vec3T& inward) {
  int nr=(int)f._boundary.size();
  bool hasIn=false,hasOut=false;
  std::vector<T> inside(nr);
  for(int i=0; i<nr; i++) {
    inside[i]=(f._boundary[i]-pos).dot(inward);
    hasIn |=inside[i]>0;
    hasOut|=inside[i]<=0;
  }
  if(!hasIn) {
    f._boundary.clear();
    return;
  }
  if(!hasOut)
    return;
  //Sutherland–Hodgman algorithm
  std::vector<Vec3T> boundaryNew;
  for(int curr=0, next=1; curr<nr; curr++, next=(curr+1)%nr)
    if(inside[curr]>0 && inside[next]>0)            //inside->inside
      boundaryNew.push_back(f._boundary[next]);
    else if(inside[curr]>0 && inside[next]<=0) {    //inside->outside
      if(abs(inside[curr]-inside[next])<Epsilon<double>::defaultEps())   //safety check
        continue;
      T alpha=inside[curr]/(inside[curr]-inside[next]);
      boundaryNew.push_back(interp1D<Vec3T,T>(f._boundary[curr],f._boundary[next],alpha));
    } else if(inside[curr]<=0 && inside[next]>0) {  //outside->inside
      if(abs(inside[curr]-inside[next])<Epsilon<double>::defaultEps())   //safety check
        continue;
      T alpha=inside[curr]/(inside[curr]-inside[next]);
      boundaryNew.push_back(interp1D<Vec3T,T>(f._boundary[curr],f._boundary[next],alpha));
      boundaryNew.push_back(f._boundary[next]);
    } //outside->outside
  boundaryNew.swap(f._boundary);
}
void SAT::clip(Facet& f,const Facet& ref) {
  int nr=(int)ref._boundary.size();
  for(int i=0; i<nr; i++) {
    Vec3T dir=ref._boundary[(i+1)%nr]-ref._boundary[i];
    clip(f,ref._boundary[i],ref._n.cross(dir));
    if(f._boundary.empty())
      return;
  }
}
}
