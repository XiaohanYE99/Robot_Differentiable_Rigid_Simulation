#ifndef SHAPE_EXACT_H
#define SHAPE_EXACT_H

#include "BVHNode.h"
#include <Utils/Pragma.h>

namespace PHYSICSMOTION {
template <typename T>
struct VTKWriter;
struct BBoxExact;
struct SphericalBBoxExact;
//#define USE_RATIONAL
#ifdef USE_RATIONAL
#define GEOMETRY_SCALAR rational
#else
#define  GEOMETRY_SCALAR FLOAT
#endif
struct ShapeExact : public SerializableBase {
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  struct Facet {
    std::vector<Vec3T> _boundary;
    Vec3T _n;
  };
  struct Edge {
    Vec3T _a,_b;
  };
  virtual const BBoxExact& getBB() const=0;
  virtual bool empty() const=0;
  virtual void getMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                       std::vector<Eigen::Matrix<int,3,1>>& iss) const=0;
  virtual bool closestInner(const Vec3T& pt,Vec3T& n,Vec3T& normal,Mat3T& hessian,
                            T& rad,Eigen::Matrix<int,2,1>& feat,bool cache=false,
                            std::vector<Vec3T>* history=NULL) const=0;
  virtual void scale(T coef)=0;
  template <typename T2>
  T2 closest(const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& n,Eigen::Matrix<T2,3,1>& normal,Eigen::Matrix<T2,3,3>& hessian,
             Eigen::Matrix<int,2,1>& feat,bool cache=false,
             std::vector<Eigen::Matrix<T2,3,1>>* history=NULL) const {
    bool inside;
    T rad=0;
    Mat3T hessianR;
    Vec3T ptR=pt.template cast<T>(),nR,normalR;
    std::vector<Vec3T> historyT;
    inside=closestInner(ptR,nR,normalR,hessianR,rad,feat,cache,history?&historyT:NULL);
    if(history) {
      history->resize(historyT.size());
      for(int i=0; i<(int)history->size(); i++)
        history->at(i)=historyT[i].template cast<T2>();
    }
    //cast
    n=nR.template cast<T2>();
    normal=normalR.template cast<T2>();
    hessian=hessianR.template cast<T2>();
    //post process
    T2 nLen=n.norm();
    hessian/=std::max(std::numeric_limits<T2>::epsilon(),nLen);
    if(nR[0]==0 && nR[1]==0 && nR[2]==0) {
      nLen=0;
      normal/=std::max(std::numeric_limits<T2>::epsilon(),normal.norm());
      //handle additional radius
      if(rad>0) {
        n+=normal*(T2)rad;
        nLen=(T2)-rad;
        inside=true;
      }
    } else {
      normal=n/std::max(std::numeric_limits<T2>::epsilon(),nLen);
      //handle additional radius
      if(rad>0) {
        if(inside) {
          n+=normal*(T2)rad;
        } else {
          inside=nLen<(T2)rad;
          n-=normal*(T2)rad;
        }
        nLen=n.norm();
      }
      if(inside)
        nLen*=-1;
      else {
        normal*=-1;
        hessian*=-1;
      }
    }
    return nLen;
  }
  //for GJK
  virtual Vec3T support(const Vec3T& D,int& id) const;
  //for SAT
  virtual Vec2T project(const Vec3T& d) const;
  virtual std::vector<Facet> facets() const;
  virtual std::vector<Edge> edges() const;
  void checkPositiveFacets() const;
  virtual void writeVTK(VTKWriter<double>& os,const Mat3X4T& trans) const=0;
  virtual void writeFacetsVTK(VTKWriter<double>& os,const Mat3X4T& trans) const;
  virtual void writeEdgesVTK(VTKWriter<double>& os,const Mat3X4T& trans) const;
};
}

#endif
