#ifndef GJK_H
#define GJK_H

#include "ShapeExact.h"
#include "PointCloudExact.h"
#include <stack>

namespace PHYSICSMOTION {
struct GJK {
  typedef GEOMETRY_SCALAR T;
  DECL_MAT_VEC_MAP_TYPES_T
  struct GJKPoint {
    void calculate(const Mat3X4T& transA,const Mat3X4T& transB);
    Vec3T _ptAL,_ptBL,_ptAB;
    int _idA,_idB;
  };
  template <typename T2>
  struct PDEntry {
    Eigen::Matrix<T2,3,1> _pAL,_pBL;
    int _idBVH;
    T2 _PD;
  };
  static Vec3T computeD(const GJKPoint v[4],int nrP,T* bary,
                        const Mat3X4T& transA,
                        const Mat3X4T& transB,
                        Vec3T& pAL,Vec3T& pBL);
  static T runGJK(std::shared_ptr<ShapeExact> A,
                  std::shared_ptr<ShapeExact> B,
                  const Mat3X4T& transA,
                  const Mat3X4T& transB,
                  Vec3T& pAL,Vec3T& pBL,bool* intersect=NULL);
  template <typename T2>
  static T2 closest(std::shared_ptr<MeshExact> A,const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& n) {
    Eigen::Matrix<T2,3,1> normal;
    Eigen::Matrix<T2,3,3> hessian;
    Eigen::Matrix<int,2,1> feat;
    return A->template closest<T2>(pt,n,normal,hessian,feat);
  }
  template <typename T2>
  static T2 closest(std::shared_ptr<PointCloudExact> A,const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& n) {
    Eigen::Matrix<T2,3,1> normal;
    n.setZero();
    return A->template closestSDF<T2>(pt,normal);
  }
  template <typename T2,typename TYPE_ROBOT=MeshExact>
  static T2 runPD(std::shared_ptr<TYPE_ROBOT> A,
                  std::shared_ptr<PointCloudExact> B,
                  const Eigen::Matrix<T2,3,4>& transA,
                  const Eigen::Matrix<T2,3,4>& transB,
                  Eigen::Matrix<T2,3,1>& pAL,Eigen::Matrix<T2,3,1>& pBL) {
    Eigen::Matrix<T2,3,4> invTA,invTATB;
    INV(invTA,transA)
    APPLY_TRANS(invTATB,invTA,transB)
    const std::vector<Node<int,BBoxExact>>& bvh=B->getBVH();
    //query function
    auto queryPD=[&](int idBVH)->PDEntry<T2> {
      Eigen::Matrix<T2,3,1> n;
      PDEntry<T2> entry;
      entry._idBVH=idBVH;
      if(bvh[idBVH]._cell>=0)
        entry._pBL=B->vss()[bvh[idBVH]._cell].template cast<T2>();
      else entry._pBL=(bvh[idBVH]._bb.minCorner()+bvh[idBVH]._bb.maxCorner()).template cast<T2>()/2;
      entry._pAL=ROT(invTATB)*entry._pBL+CTR(invTATB);
      entry._PD=-closest<T2>(A,entry._pAL,n);
      entry._pAL+=n;
      return entry;
    };
    //main loop
    T2 rad;
    std::stack<PDEntry<T2>> ss;
    PDEntry<T2> maxPD,entry,entryL,entryR;
    ss.push(queryPD((int)bvh.size()-1));
    maxPD._PD=-std::numeric_limits<double>::max();
    while(!ss.empty()) {
      entry=ss.top();
      ss.pop();
      if(bvh[entry._idBVH]._cell>=0) {
        if(entry._PD>maxPD._PD)
          maxPD=entry;
      } else {
        const Node<int,BBoxExact>& node=bvh[entry._idBVH];
        rad=(node._bb.maxCorner()-node._bb.minCorner()).template cast<T2>().norm()/2;
        if(entry._PD+(T2)rad<maxPD._PD)
          continue;
        else {
          //first explore node with larger PD
          entryL=queryPD(node._l);
          entryR=queryPD(node._r);
          if(entryL._PD>entryR._PD)
            std::swap(entryL,entryR);
          ss.push(entryL);
          ss.push(entryR);
        }
      }
    }
    pAL=maxPD._pAL;
    pBL=maxPD._pBL;
    return maxPD._PD;
  }
  static void writeVTK(const std::string& path,
                       std::shared_ptr<MeshExact> A,
                       std::shared_ptr<MeshExact> B,
                       const Mat3X4T& transA,
                       const Mat3X4T& transB);
  static void writeVTK(const std::string& path,
                       std::shared_ptr<MeshExact> A,
                       std::shared_ptr<PointCloudExact> B,
                       const Mat3X4T& transA,
                       const Mat3X4T& transB);
};
}
#endif
