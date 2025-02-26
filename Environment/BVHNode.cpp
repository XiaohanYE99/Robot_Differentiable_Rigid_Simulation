#include "BVHNode.h"
#include "BBoxExact.h"
#include <Utils/Heap.h>
#include <Utils/IO.h>
#include <stack>

namespace PHYSICSMOTION {
template <int D>
struct SurfaceArea;
template <>
struct SurfaceArea<3> {
  template <typename BoxType>
  static typename BoxType::Vec3T::Scalar area(const BoxType& bb) {
    typedef typename BoxType::Vec3T Vec3T;
    Vec3T ext=Vec3T::Zero().cwiseMax(bb.maxCorner()-bb.minCorner());
    return (ext[0]*ext[1]+ext[0]*ext[2]+ext[1]*ext[2])*2.0f;
  }
};
template <>
struct SurfaceArea<2> {
  template <typename BoxType>
  static typename BoxType::Vec3T::Scalar area(const BoxType& bb) {
    typedef typename BoxType::Vec3T Vec3T;
    Vec3T ext=Vec3T::Zero().cwiseMax(bb.maxCorner()-bb.minCorner());
    return (ext[0]+ext[1])*2.0f;
  }
};
//Node
template <typename T,typename BBOX>
Node<T,BBOX>::Node():_l(-1),_r(-1),_parent(-1),_nrCell(-1) {}
template <typename T,typename BBOX>
bool Node<T,BBOX>::read(std::istream& is,IOData* dat) {
  readBinaryData(_bb,is);
  readBinaryData(_cell,is,dat);
  readBinaryData(_l,is);
  readBinaryData(_r,is);
  readBinaryData(_parent,is);
  readBinaryData(_nrCell,is);
  return is.good();
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_bb,os);
  writeBinaryData(_cell,os,dat);
  writeBinaryData(_l,os);
  writeBinaryData(_r,os);
  writeBinaryData(_parent,os);
  writeBinaryData(_nrCell,os);
  return os.good();
}
template <typename T,typename BBOX>
std::shared_ptr<SerializableBase> Node<T,BBOX>::copy() const {
  return std::shared_ptr<SerializableBase>(new Node<T,BBOX>);
}
template <typename T,typename BBOX>
std::string Node<T,BBOX>::type() const {
  return typeid(Node<T,BBOX>).name();
}
template <typename T,typename BBOX>
Node<T,BBOX>& Node<T,BBOX>::operator=(const Node<T,BBOX>& other) {
  _bb=other._bb;
  _cell=other._cell;
  _l=other._l;
  _r=other._r;
  _parent=other._parent;
  _nrCell=other._nrCell;
  return* this;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::buildBVHTriangleBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge) {
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  buildEdge(iss,edgeMap);
  buildBVHBottomUp(bvh,edgeMap,forceMerge);
}
template <typename T,typename BBOX>
void Node<T,BBOX>::buildBSHTriangleBottomUp(std::vector<Node<T,BBOX>>& bsh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge) {
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  buildEdge(iss,edgeMap);
  buildBSHBottomUp(bsh,edgeMap,forceMerge);
}
template <typename T,typename BBOX>
void Node<T,BBOX>::buildBVHVertexBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::vector<Eigen::Matrix<int,3,1>>& iss,bool forceMerge) {
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  for(const Eigen::Matrix<int,3,1>& t:iss)
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,2,1> e(t[(d+1)%3],t[(d+2)%3]);
      if(e[0]>e[1])
        std::swap(e[0],e[1]);
      edgeMap[e]=std::make_pair(e[0],e[1]);
    }
  buildBVHBottomUp(bvh,edgeMap,forceMerge);
}
template <typename T,typename BBOX>
void Node<T,BBOX>::buildBVHBottomUp(std::vector<Node<T,BBOX>>& bvh,const std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap,bool forceMerge) {
  int nrLeaf=0;
  for(const Node<T,BBOX>& n:bvh)
    if(n._l==-1)
      nrLeaf++;
  //initialize hash
  std::vector<int> heap;
  std::vector<int> heapOffsets;
  std::vector<std::pair<int,int>> ess;
  std::vector<typename BBOX::Vec3T::Scalar> cost;
  for(const auto& e:edgeMap) {
    heapOffsets.push_back(-1);
    ess.push_back(e.second);
    BBOX bb=bvh[e.second.first]._bb;
    if(e.second.second>=0)
      bb.setUnion(bvh[e.second.second]._bb);
    typename BBOX::Vec3T::Scalar c=SurfaceArea<3>::area(bb);
    cost.push_back(c);
  }
  for(int i=0; i<(int)ess.size(); i++)
    pushHeapDef(cost,heapOffsets,heap,i);
  //merge BVH
  int err;
  while(!heap.empty()) {
    int i=popHeapDef(cost,heapOffsets,heap,err);
    int t0=ess[i].first,t1=ess[i].second;
    //boundary edge
    if(t1==-1)
      continue;
    //find parent
    while(bvh[t0]._parent>=0)
      t0=bvh[t0]._parent;
    while(bvh[t1]._parent>=0)
      t1=bvh[t1]._parent;
    //check already merged
    if(t0==t1)
      continue;
    //merge
    BBOX bb=bvh[t0]._bb;
    bb.setUnion(bvh[t1]._bb);
    typename BBOX::Vec3T::Scalar c=SurfaceArea<3>::area(bb);
    if(c>cost[i]) {
      cost[i]=c;
      pushHeapDef(cost,heapOffsets,heap,i);
    } else {
      Node<T,BBOX> n;
      n._l=t0;
      n._r=t1;
      n._parent=-1;
      n._cell=-1;
      n._bb=bb;
      n._nrCell=bvh[n._l]._nrCell+bvh[n._r]._nrCell;
      bvh[t0]._parent=(int)bvh.size();
      bvh[t1]._parent=(int)bvh.size();
      bvh.push_back(n);
    }
  }
  if(forceMerge)
    buildBVHBottomUpAll(bvh);
  else {
    ASSERT_MSG((int)bvh.size()==nrLeaf*2-1,"Multi-component mesh detected!")
  }
}

template <typename T,typename BBOX>
void Node<T,BBOX>::buildBSHBottomUp(std::vector<Node<T,BBOX>>& bsh,const std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap,bool forceMerge) {
  int nrLeaf=0;
  for(const Node<T,BBOX>& n:bsh)
    if(n._l==-1)
      nrLeaf++;
  //initialize hash
  std::vector<int> heap;
  std::vector<int> heapOffsets;
  std::vector<std::pair<int,int>> ess;
  std::vector<typename BBOX::Vec3T::Scalar> cost;
  for(const auto& e:edgeMap) {
    heapOffsets.push_back(-1);
    ess.push_back(e.second);
    BBOX bb=bsh[e.second.first]._bb;
    if(e.second.second>=0)
      bb.setUnion(bsh[e.second.second]._bb);
    typename BBOX::Vec3T::Scalar c=SurfaceArea<3>::area(bb);
    cost.push_back(c);
  }
  for(int i=0; i<(int)ess.size(); i++)
    pushHeapDef(cost,heapOffsets,heap,i);
  //merge BVH
  int err;
  while(!heap.empty()) {
    int i=popHeapDef(cost,heapOffsets,heap,err);
    int t0=ess[i].first,t1=ess[i].second;
    //boundary edge
    if(t1==-1)
      continue;
    //find parent
    while(bsh[t0]._parent>=0)
      t0=bsh[t0]._parent;
    while(bsh[t1]._parent>=0)
      t1=bsh[t1]._parent;
    //check already merged
    if(t0==t1)
      continue;
    //merge
    BBOX bb=bsh[t0]._bb;
    bb.setUnion(bsh[t1]._bb);
    typename BBOX::Vec3T::Scalar c=SurfaceArea<3>::area(bb);
    if(c>cost[i]) {
      cost[i]=c;
      pushHeapDef(cost,heapOffsets,heap,i);
    } else {
      Node<T,BBOX> n;
      n._l=t0;
      n._r=t1;
      n._parent=-1;
      n._cell=-1;
      n._bb=bb;
      n._nrCell=bsh[n._l]._nrCell+bsh[n._r]._nrCell;
      bsh[t0]._parent=(int)bsh.size();
      bsh[t1]._parent=(int)bsh.size();
      bsh.push_back(n);
    }
  }
  if(forceMerge)
    buildBSHBottomUpAll(bsh);
  else {
    ASSERT_MSG((int)bsh.size()==nrLeaf*2-1,"Multi-component mesh detected!")
  }
}

template <typename T,typename BBOX>
void Node<T,BBOX>::buildBVHBottomUpAll(std::vector<Node<T,BBOX>>& bvh) {
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  std::vector<int> roots;
  for(int i=0; i<(int)bvh.size(); i++)
    if(bvh[i]._parent==-1)
      roots.push_back(i);
  if((int)roots.size()>1) {
    std::cout << "Merging " << roots.size() << " isolated BVH!" << std::endl;
    for(int i=0; i<(int)roots.size(); i++)
      for(int j=i+1; j<(int)roots.size(); j++)
        edgeMap[Eigen::Matrix<int,2,1>(roots[i],roots[j])]=std::make_pair(roots[i],roots[j]);
    buildBVHBottomUp(bvh,edgeMap,false);
  }
}
template <typename T,typename BBOX>
void Node<T,BBOX>::buildBSHBottomUpAll(std::vector<Node<T,BBOX>>& bvh) {
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  std::vector<int> roots;
  for(int i=0; i<(int)bvh.size(); i++)
    if(bvh[i]._parent==-1)
      roots.push_back(i);
  if((int)roots.size()>1) {
    std::cout << "Merging " << roots.size() << " isolated BVH!" << std::endl;
    for(int i=0; i<(int)roots.size(); i++)
      for(int j=i+1; j<(int)roots.size(); j++)
        edgeMap[Eigen::Matrix<int,2,1>(roots[i],roots[j])]=std::make_pair(roots[i],roots[j]);
    buildBSHBottomUp(bvh,edgeMap,false);
  }
}
template <typename T,typename BBOX>
int Node<T,BBOX>::nrEmpty(const std::vector<Node<T,BBOX>>& bvh) {
  if(bvh.empty() || !isEmpty(bvh,0))
    return 0;
  int curr=0,nr=1;
  while(bvh[curr]._parent>=0) {
    curr=bvh[curr]._parent;
    nr++;
  }
  return nr;
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::isEmpty(const std::vector<Node<T,BBOX>>& bvh,int id) {
  return bvh[id]._cell==-1 && bvh[id]._l==-1 && bvh[id]._r==-1;
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::hasEmpty(const std::vector<Node<T,BBOX>>& bvh,int nr0) {
  if(bvh.empty() || !isEmpty(bvh,0))
    return 0>=nr0;
  int curr=0,nr=1;
  if(nr>=nr0)
    return true;
  while(bvh[curr]._parent>=0) {
    curr=bvh[curr]._parent;
    nr++;
    if(nr>=nr0)
      return true;
  }
  return false;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::addEmpty(std::vector<Node<T,BBOX>>& bvh,int id) {
  int tmp=bvh[0]._parent;
  bvh[0]._parent=id;
  bvh[id]._parent=tmp;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::reserveEmpty(std::vector<Node<T,BBOX>>& bvh,int nr) {
  int sz=(int)bvh.size();
  bvh.resize(sz+nr);
  // reserve
  std::copy_backward(bvh.begin(),bvh.begin()+sz,bvh.end());
  for(int i=nr; i<nr+sz; i++) {
    if(bvh[i]._l>=0)
      bvh[i]._l+=nr;
    if(bvh[i]._r>=0)
      bvh[i]._r+=nr;
    if(bvh[i]._parent>=0)
      bvh[i]._parent+=nr;
  }
  // reconnect
  for(int i=0; i<nr; i++) {
    bvh[i]._l=bvh[i]._r=bvh[i]._parent=-1;
    bvh[i]._cell=-1;
    bvh[i]._parent=i+1;
  }
  if((int)bvh.size()>nr && !isEmpty(bvh,nr))
    bvh[nr-1]._parent=-1;
}
template <typename T,typename BBOX>
int Node<T,BBOX>::getEmpty(std::vector<Node<T,BBOX>>& bvh) {
  ASSERT_MSG(isEmpty(bvh,0),"No empty nodes to use!")
  if(isEmpty(bvh,(int)bvh.size()-1)) {
    bvh[bvh.size()-2]._parent=-1;
    return (int)bvh.size()-1;
  } else if(bvh[0]._parent>=0) {
    int ret=bvh[0]._parent;
    bvh[0]._parent=bvh[ret]._parent;
    return ret;
  } else return 0;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::recomputeNrCell(std::vector<Node<T,BBOX>>& bvh,int pid) {
  while((pid=bvh[pid]._parent)>=0)
    bvh[pid]._nrCell=bvh[bvh[pid]._l]._nrCell+bvh[bvh[pid]._r]._nrCell;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::updateBVH(std::vector<Node<T,BBOX>>& bvh,int i) {
  while(i>=0) {
    bvh[i]._bb=BBOX();
    bvh[i]._bb.setUnion(bvh[bvh[i]._l]._bb);
    bvh[i]._bb.setUnion(bvh[bvh[i]._r]._bb);
    i=bvh[i]._parent;
  }
}
template <typename T,typename BBOX>
void Node<T,BBOX>::insertLeaf(std::vector<Node<T,BBOX>>& bvh,T val,int& i,std::function<void(Node<T,BBOX>&)> updateNode,int reserveBatch) {
  if(bvh.empty() || !hasEmpty(bvh,2)) {
    ASSERT_MSG(reserveBatch>0,"reserveBatch<=0, but we need more empty nodes!")
    reserveEmpty(bvh,reserveBatch);
    i+=reserveBatch;
  }

  if(i==(int)bvh.size()) {
    Node<T,BBOX>& c=bvh[getEmpty(bvh)];
    c._cell=val;
    c._nrCell=1;
    c._parent=-1;
    updateNode(c);
    return;
  }

  Node<T,BBOX>& c=bvh[i];
  if(c._cell==-1)
    std::cout << c._cell << std::endl;
  ASSERT_MSG(c._cell!=-1 && val!=-1,"Wrong cell!")
  //find parent
  c._l=getEmpty(bvh);
  c._r=getEmpty(bvh);
  Node<T,BBOX>& cl=bvh[c._l];
  cl._l=cl._r=-1;
  cl._parent=i;
  cl._cell=val;
  cl._nrCell=1;
  updateNode(cl);
  Node<T,BBOX>& cr=bvh[c._r];
  cr._l=cr._r=-1;
  cr._parent=i;
  cr._cell=c._cell;
  cr._nrCell=c._nrCell;
  updateNode(cr);
  //recompute
  //c._cell=-1; //we do not clear the cell of parent node
  c._nrCell=cl._nrCell+cr._nrCell;
  recomputeNrCell(bvh,i);
  updateBVH(bvh,i);
}
template <typename T,typename BBOX>
void Node<T,BBOX>::removeLeaf(std::vector<Node<T,BBOX>>& bvh,int i) {
  //make sure that first node is empty
  if(!hasEmpty(bvh,1)) {
    reserveEmpty(bvh,1);
    i+=1;
  }
  Node<T,BBOX>& c=bvh[i];
  ASSERT(c._cell!=-1)
  //find parent
  if(c._parent==-1)
    bvh.clear();
  else {
    //use parent as sibling
    Node<T,BBOX>& p=bvh[c._parent];
    int sid=p._l==i?p._r:p._l;
    Node<T,BBOX>& s=bvh[sid];
    p._bb=s._bb;
    p._cell=s._cell;
    p._nrCell=s._nrCell;
    recomputeNrCell(bvh,c._parent);
    p._l=s._l;
    if(p._l>=0)
      bvh[p._l]._parent=c._parent;
    p._r=s._r;
    if(p._r>=0)
      bvh[p._r]._parent=c._parent;
    //remove sibling
    s._l=s._r=s._parent=-1;
    s._cell=-1;
    addEmpty(bvh,sid);
    //remove current
    c._l=c._r=c._parent=-1;
    c._cell=-1;
    addEmpty(bvh,i);
  }
}
template <typename T,typename BBOX>
void Node<T,BBOX>::compact(std::vector<Node<T,BBOX>>& bvh) {
  std::vector<Node<T,BBOX>> tmp(bvh.size()-nrEmpty(bvh));
  if(tmp.empty()) {
    bvh.clear();
    return;
  }
  //assign
  std::stack<std::pair<int,int>> ss;
  tmp[(int)tmp.size()-1]=bvh.back();
  ss.push(std::make_pair((int)bvh.size()-1,(int)tmp.size()-1));
  int curr=ss.top().second;
  while(!ss.empty())
    if(bvh[ss.top().first]._cell==-1) {
      int l=bvh[ss.top().first]._l;
      int r=bvh[ss.top().first]._r;
      tmp[curr-1]=bvh[l];
      tmp[curr-2]=bvh[r];
      tmp[ss.top().second]._l=curr-1;
      tmp[ss.top().second]._r=curr-2;
      tmp[curr-1]._parent=ss.top().second;
      tmp[curr-2]._parent=ss.top().second;
      ss.pop();
      ss.push(std::make_pair(l,curr-1));
      ss.push(std::make_pair(r,curr-2));
      curr-=2;
    } else ss.pop();
  const_cast<std::vector<Node<T,BBOX>>&>(bvh)=tmp;
}
template <typename T,typename BBOX>
void Node<T,BBOX>::reconnect(std::vector<Node<T,BBOX>>& bvh,int id) {
  bvh[bvh[id]._l]._parent=id;
  bvh[bvh[id]._r]._parent=id;
  bvh[id]._bb=bvh[bvh[id]._l]._bb;
  bvh[id]._bb.setUnion(bvh[bvh[id]._r]._bb);
  bvh[id]._nrCell=bvh[bvh[id]._l]._nrCell+bvh[bvh[id]._r]._nrCell;
}
template <typename T,typename BBOX>
typename BBOX::Vec3T::Scalar Node<T,BBOX>::surfaceArea(const std::vector<Node<T,BBOX>>& bvh,int bid) {
  return SurfaceArea<3>::area(bvh[bid]._bb);
}
template <typename T,typename BBOX>
typename BBOX::Vec3T::Scalar Node<T,BBOX>::surfaceArea(const std::vector<Node<T,BBOX>>& bvh,int l,int r) {
  BBOX bb=BBOX();
  bb.setUnion(bvh[l]._bb);
  bb.setUnion(bvh[r]._bb);
  return SurfaceArea<3>::area(bb);
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::tryRotate(std::vector<Node<T,BBOX>>& bvh,int pid,int nlid) {
  int& otherId=bvh[pid]._l==nlid?bvh[pid]._r:bvh[pid]._l;
  //compute cost
  auto cost=surfaceArea(bvh,nlid)+surfaceArea(bvh,otherId);
  auto cost1=surfaceArea(bvh,bvh[nlid]._r)+surfaceArea(bvh,bvh[nlid]._l,otherId);
  auto cost2=surfaceArea(bvh,bvh[nlid]._l)+surfaceArea(bvh,bvh[nlid]._r,otherId);
  //rotate
  if(cost1<cost && cost1<cost2)
    std::swap(bvh[nlid]._r,otherId);
  else if(cost2<cost && cost2<cost1)
    std::swap(bvh[nlid]._l,otherId);
  else return false;
  //reconnect
  reconnect(bvh,nlid);
  reconnect(bvh,pid);
  //auto costFinal=surfaceArea(bvh[pid]._l)+surfaceArea(bvh[pid]._r);
  return true;
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::tryRecombine(std::vector<Node<T,BBOX>>& bvh,int pid) {
  int lid=bvh[pid]._l;
  int rid=bvh[pid]._r;
  //compute cost
  auto cost=surfaceArea(bvh,lid)+surfaceArea(bvh,rid);
  auto cost1=surfaceArea(bvh,bvh[lid]._l,bvh[rid]._l)+surfaceArea(bvh,bvh[lid]._r,bvh[rid]._r);
  auto cost2=surfaceArea(bvh,bvh[lid]._l,bvh[rid]._r)+surfaceArea(bvh,bvh[lid]._r,bvh[rid]._l);
  if(cost1<cost && cost1<cost2)
    std::swap(bvh[lid]._r,bvh[rid]._l);
  else if(cost2<cost && cost2<cost1)
    std::swap(bvh[lid]._r,bvh[rid]._r);
  else return false;
  //reconnect
  reconnect(bvh,lid);
  reconnect(bvh,rid);
  //SCALAR costFinal=surfaceArea(bvh[pid]._l)+surfaceArea(bvh[pid]._r);
  return true;
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::optimizeBVH(std::vector<Node<T,BBOX>>& bvh,typename BBOX::T expand,bool dynamic) {
  if(hasEmpty(bvh,1))
    compact(bvh);
  int nrN=(int)bvh.size();
  bool updated=false;
  for(int i=0; i<nrN; i++) {
    Node<T,BBOX>& n=const_cast<Node<T,BBOX>&>(bvh[i]);
    if(n._cell!=-1)
      n._bb.enlarged(BBOX::Vec3T::Constant(expand));
    else if(n._l>=0 && n._r>=0) {
      n._bb=bvh[n._l]._bb;
      n._bb.setUnion(bvh[n._r]._bb);
      //dynamic SAH-based adjustment
      if(!dynamic)
        continue;
      const Node<T,BBOX>& l=bvh[n._l];
      const Node<T,BBOX>& r=bvh[n._r];
      if(l._cell!=-1 && r._cell==-1)
        updated=tryRotate(bvh,i,n._r)||updated;
      else if(l._cell==-1 && r._cell!=-1)
        updated=tryRotate(bvh,i,n._l)||updated;
      else if(l._cell==-1 && r._cell==-1)
        updated=tryRecombine(bvh,i)||updated;
    }
  }
  return updated;
}
template <typename T,typename BBOX>
int Node<T,BBOX>::findLeaf(std::vector<Node<T,BBOX>>& bvh,const BBOX& bb) {
  if(bvh.empty())
    return 0;
  int curr=(int)bvh.size()-1;
  while(bvh[curr]._cell==-1) {
    BBOX bbl=bvh[bvh[curr]._l]._bb;
    BBOX bbr=bvh[bvh[curr]._r]._bb;
    auto costL=SurfaceArea<3>::area(bbr);
    auto costR=SurfaceArea<3>::area(bbl);
    bbl.setUnion(bb);
    bbr.setUnion(bb);
    costL+=SurfaceArea<3>::area(bbl);
    costR+=SurfaceArea<3>::area(bbr);
    if(costL < costR)
      curr=bvh[curr]._l;
    else curr=bvh[curr]._r;
  }
  return curr;
}
//check
template <typename T,typename BBOX>
int Node<T,BBOX>::parityCheck(const std::vector<Node<T,BBOX>>& bvh) {
  if(bvh.empty())
    return 0;
  ASSERT_MSG(bvh.back()._parent==-1,"root node failed!")
  int nrLeaf=parityCheck(bvh,(int)bvh.size()-1);
  ASSERT_MSG(bvh.back()._nrCell*2-1+nrEmpty(bvh)==(int)bvh.size(),"bvh node sum inconsistent!")
  return nrLeaf;
}
template <typename T,typename BBOX>
int Node<T,BBOX>::parityCheck(const std::vector<Node<T,BBOX>>& bvh,int i) {
  int ret=(bvh[i]._cell!=-1)?1:0;
  if(bvh[i]._l>=0) {
    ASSERT_MSG(bvh[bvh[i]._l]._parent==i && bvh[i]._cell==-1,"internal node left child failed!")
    ASSERT_MSGV(bvh[i]._nrCell==bvh[bvh[i]._l]._nrCell+bvh[bvh[i]._r]._nrCell,"%dth bvh node #cell inconsistent!",i)
    ret+=parityCheck(bvh,bvh[i]._l);
  } else {
    ASSERT_MSG(bvh[i]._cell!=-1,"leaf node failed!")
  }
  if(bvh[i]._r >= 0) {
    ASSERT_MSG(bvh[bvh[i]._r]._parent==i && bvh[i]._cell==-1,"internal node right child failed!")
    ASSERT_MSGV(bvh[i]._nrCell==bvh[bvh[i]._l]._nrCell+bvh[bvh[i]._r]._nrCell,"%dth bvh node #cell inconsistent!",i)
    ret+=parityCheck(bvh,bvh[i]._r);
  } else {
    ASSERT_MSG(bvh[i]._cell!=-1,"leaf node failed!")
  }
  return ret;
}
template <typename T,typename BBOX>
int Node<T,BBOX>::depth(const std::vector<Node<T,BBOX>>& bvh,int i) {
  if(i==-1)
    return depth(bvh,(int)bvh.size()-1);
  else if(bvh[i]._l==-1)
    return 1;
  else {
    int dl=depth(bvh,bvh[i]._l);
    int dr=depth(bvh,bvh[i]._r);
    return std::max(dl,dr)+1;
  }
}
template <typename T,typename BBOX>
bool Node<T,BBOX>::operator==(Node<T,BBOX>& node) {
  return (_l==node._l) && (_r==node._r) && (_parent==node._parent) && (_cell==node._cell) && (_nrCell==node._nrCell) && (_bb==node._bb);
}
template <typename T,typename BBOX>
void Node<T,BBOX>::print(const std::string& prefix) const {
  std::cout << prefix << "Left=" << _l << " Right=" << _r << " Parent=" << _parent << " nrCell=" << _nrCell <<   std::endl;
}
//instance
template struct Node<int,BBoxExact>;
}
