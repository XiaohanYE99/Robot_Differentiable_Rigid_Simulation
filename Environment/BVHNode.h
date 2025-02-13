#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "EnvironmentUtils.h"
#include <Utils/Serializable.h>

namespace PHYSICSMOTION {
template <typename T,typename BBOX>
struct Node : public SerializableBase {
  typedef BBOX BoxType;
  using SerializableBase::read;
  using SerializableBase::write;
  using SerializableBase::operator<;
  using SerializableBase::operator=;
  Node();
  bool read(std::istream& is,IOData *dat);
  bool write(std::ostream& os,IOData *dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  virtual std::string type() const override;
  Node<T,BBOX>& operator=(const Node<T,BBOX>& other);
  static void buildBVHTriangleBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge=false);
  static void buildBVHVertexBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge=false);
  static void buildBVHBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap,bool forceMerge=false);
  static void buildBVHBottomUpAll(std::vector<Node<T,BBOX>>& bvh);

  static void buildBSHTriangleBottomUp(std::vector<Node<T,BBOX>>& bsh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge=false);
  static void buildBSHVertexBottomUp(std::vector<Node<T,BBOX>>& bsh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge=false);
  static void buildBSHBottomUp(std::vector<Node<T,BBOX>>& bsh,const std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap,bool forceMerge=false);
  static void buildBSHBottomUpAll(std::vector<Node<T,BBOX>>& bsh);

  static int nrEmpty(const std::vector<Node<T,BBOX>>& bvh);
  static bool isEmpty(const std::vector<Node<T,BBOX>>& bvh,int id);
  static bool hasEmpty(const std::vector<Node<T,BBOX>>& bvh,int nr0);
  static void addEmpty(std::vector<Node<T,BBOX>>& bvh,int id);
  static void reserveEmpty(std::vector<Node<T,BBOX>>& bvh,int nr);
  static int getEmpty(std::vector<Node<T,BBOX>>& bvh);
  static void recomputeNrCell(std::vector<Node<T,BBOX>>& bvh,int pid);
  static void updateBVH(std::vector<Node<T,BBOX>>& bvh,int i);
  static void insertLeaf(std::vector<Node<T,BBOX>>& bvh,T val,int& i,std::function<void(Node<T,BBOX>&)> updateNode,int reserveBatch=256);
  static void removeLeaf(std::vector<Node<T,BBOX>>& bvh,int i);
  static void compact(std::vector<Node<T,BBOX>>& bvh);
  static void reconnect(std::vector<Node<T,BBOX>>& bvh,int id);
  static typename BBOX::Vec3T::Scalar surfaceArea(const std::vector<Node<T,BBOX>>& bvh,int bid);
  static typename BBOX::Vec3T::Scalar surfaceArea(const std::vector<Node<T,BBOX>>& bvh,int l,int r);
  static bool tryRotate(std::vector<Node<T,BBOX>>& bvh,int pid,int nlid);
  static bool tryRecombine(std::vector<Node<T,BBOX>>& bvh,int pid);
  static bool optimizeBVH(std::vector<Node<T,BBOX>>& bvh,typename BBOX::T expand=0,bool dynamic=true);
  static int findLeaf(std::vector<Node<T,BBOX>>& bvh,const BBOX& bb);
  static int parityCheck(const std::vector<Node<T,BBOX>>& bvh);
  static int parityCheck(const std::vector<Node<T,BBOX>>& bvh,int i);
  static int depth(const std::vector<Node<T,BBOX>>& bvh,int i=-1);
  bool operator==(Node<T,BBOX>& node);
  void print(const std::string& prefix="") const;
  int _l,_r,_parent,_nrCell;
  int _num;//number of triangles
  BBOX _bb;
  T _cell;
};
}

#endif
