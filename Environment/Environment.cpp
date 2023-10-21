#include "Environment.h"
#include <Utils/DebugGradient.h>
#include <Utils/Interp.h>
#include <Utils/IO.h>
#include <random>

namespace PHYSICSMOTION {
//Environment
template <typename T>
T Environment<T>::phi(const Vec3T& x) const {
  return phi(x,NULL);
}
template <typename T>
typename Environment<T>::Vec3T Environment<T>::phiGrad(const Vec3T& x) const {
  return phiGrad(x,NULL);
}
template <typename T>
void Environment<T>::createHills(double xMin,double xMax,double yMin,double yMax,EnvironmentCallback& h,double res) {
  createHills(xMin,xMax,yMin,yMax,[&](double x,double y)->double {
    return h.height(x,y);
  },res);
}
template <typename T>
void Environment<T>::createStair(double x,double y,double x0,double z0,double slope,int n) {
#define RES 8
  std::vector<Eigen::Matrix<double,3,1>> vss;
  vss.push_back(Eigen::Matrix<double,3,1>(0,0,0));
  vss.push_back(Eigen::Matrix<double,3,1>(x,0,0));
  for(int i=1; i<=n; i++) {
    vss.push_back(Eigen::Matrix<double,3,1>(x+x0*(i-1)+slope*x0,0,i*z0));
    vss.push_back(Eigen::Matrix<double,3,1>(x+x0*i,0,i*z0));
  }
  double xMin=vss.front()[0];
  double xMax=vss.back()[0];
  createHills(xMin,xMax,-y,y,[&](double x,double)->double {
    if(x<vss.front().x())
      return vss.front().z();
    else if(x>=vss.back().x())
      return vss.back().z();
    else
      for(int i=0; i<(int)vss.size()-1; i++)
        if(x>=vss[i].x() && x<vss[i+1].x()) {
          double alpha=(x-vss[i].x())/(vss[i+1].x()-vss[i].x());
          return interp1D(vss[i].z(),vss[i+1].z(),alpha);
        }
    ASSERT(false)
    return (double)0;
  },std::min((xMax-xMin)/(n*RES),2*y/(n*RES)));
#undef RES
}
template <typename T>
void Environment<T>::createFloor(const Eigen::Matrix<double,4,1>& plane) {
  double sz=plane.template segment<3>(0).norm();
  createHills(-sz,sz,-sz,sz,[&](double x,double y) {
    //(x,y,z).dot(n)+d=0
    return -(Eigen::Matrix<double,2,1>(x,y).dot(plane.segment<2>(0))+plane[3])/plane[2];
  },sz/4);
}
//EnvironmentExact
template <typename T>
EnvironmentExact<T>::EnvironmentExact() {}
template <typename T>
EnvironmentExact<T>::EnvironmentExact(const MeshExact& obj):_obj(obj) {}
template <typename T>
bool EnvironmentExact<T>::read(std::istream& is,IOData* dat) {
  _obj.read(is,dat);
  return is.good();
}
template <typename T>
bool EnvironmentExact<T>::write(std::ostream& os,IOData* dat) const {
  _obj.write(os,dat);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> EnvironmentExact<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new EnvironmentExact(_obj));
}
template <typename T>
std::string EnvironmentExact<T>::type() const {
  return typeid(EnvironmentExact).name();
}
template <typename T>
T EnvironmentExact<T>::phi(const Vec3T& x,Vec3T* g) const {
  Mat3T hessian;
  Vec3T n,normal;
  Eigen::Matrix<int,2,1> feat;
  T ret=_obj.closest<T>(x,n,normal,hessian,feat);
  if(g)
    *g=normal;
  return ret;
}
template <typename T>
typename EnvironmentExact<T>::Vec3T EnvironmentExact<T>::phiGrad(const Vec3T& x,Mat3T* h) const {
  Mat3T hessian;
  Vec3T n,normal;
  Eigen::Matrix<int,2,1> feat;
  _obj.closest<T>(x,n,normal,hessian,feat);
  if(h)
    *h=hessian;
  return normal;
}
template <typename T>
const std::vector<Node<int,BBoxExact>>& EnvironmentExact<T>::getBVH() const {
  return _obj.getBVH();
}
template <typename T>
const BBoxExact& EnvironmentExact<T>::getBB() const {
  return _obj.getBB();
}
template <typename T>
bool EnvironmentExact<T>::empty() const {
  return _obj.empty();
}
template <typename T>
void EnvironmentExact<T>::createHills(double xMin,double xMax,double yMin,double yMax,std::function<double(double,double)> h,double res) {
#define GI(X,Y) (X)*(resY+1)+(Y)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  int resX=ceil((xMax-xMin)/res);
  int resY=ceil((yMax-yMin)/res);
  vss.resize((resX+1)*(resY+1));
  for(int r=0; r<=resX; r++)
    for(int c=0; c<=resY; c++) {
      double xx=xMin+(xMax-xMin)*r/resX;
      double yy=yMin+(yMax-yMin)*c/resY;
      vss.at(GI(r,c))=Eigen::Matrix<double,3,1>(xx,yy,h(xx,yy));
      if(r<resX && c<resY) {
        iss.push_back(Eigen::Matrix<int,3,1>(GI(r,c),GI(r+1,c),GI(r+1,c+1)));
        iss.push_back(Eigen::Matrix<int,3,1>(GI(r,c),GI(r+1,c+1),GI(r,c+1)));
      }
    }
  _obj=MeshExact(vss,iss);
#undef GI
}
template <typename T>
const MeshExact& EnvironmentExact<T>::getMesh() const {
  return _obj;
}
//EnvironmentHeight
template <typename T>
EnvironmentHeight<T>::EnvironmentHeight() {}
template <typename T>
bool EnvironmentHeight<T>::read(std::istream& is,IOData* dat) {
  readBinaryData(_height,is);
  readBinaryData(_bvh,is);
  readBinaryData(_invCellSzX,is);
  readBinaryData(_invCellSzY,is);
  readBinaryData(_resX,is);
  readBinaryData(_resY,is);
  _bb.read(is,dat);
  return is.good();
}
template <typename T>
bool EnvironmentHeight<T>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_height,os);
  writeBinaryData(_bvh,os);
  writeBinaryData(_invCellSzX,os);
  writeBinaryData(_invCellSzY,os);
  writeBinaryData(_resX,os);
  writeBinaryData(_resY,os);
  _bb.write(os,dat);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> EnvironmentHeight<T>::copy() const {
  std::shared_ptr<EnvironmentHeight<T>> ret(new EnvironmentHeight);
  *ret=*this;
  return ret;
}
template <typename T>
std::string EnvironmentHeight<T>::type() const {
  return typeid(EnvironmentHeight).name();
}
template <typename T>
T EnvironmentHeight<T>::phi(const Vec3T& x,Vec3T* g) const {
  T val[4][4];
  Vec2T pos;
  Eigen::Matrix<int,2,1> posId=getPosId(x,pos);
  for(int x=0; x<4; x++)
    for(int y=0; y<4; y++)
      val[x][y]=getHeight(posId+Eigen::Matrix<int,2,1>(x-1,y-1));
  //gradient
  if(g) {
    g->coeffRef(0)=-interp1DCubic<T>(pos[1],[&](int yId) {
      return interp1DCubicDiff<T>(pos[0],[&](int xId) {
        return val[xId][yId];
      });
    })*_invCellSzX;
    g->coeffRef(1)=-interp1DCubicDiff<T>(pos[1],[&](int yId) {
      return interp1DCubic<T>(pos[0],[&](int xId) {
        return val[xId][yId];
      });
    })*_invCellSzY;
    g->coeffRef(2)=1;
  }
  //value
  return x[2]-interp1DCubic<T>(pos[1],[&](int yId) {
    return interp1DCubic<T>(pos[0],[&](int xId) {
      return val[xId][yId];
    });
  });
}
template <typename T>
typename EnvironmentHeight<T>::Vec3T EnvironmentHeight<T>::phiGrad(const Vec3T& x,Mat3T* h) const {
  T val[4][4];
  Vec2T pos;
  Eigen::Matrix<int,2,1> posId=getPosId(x,pos);
  for(int x=0; x<4; x++)
    for(int y=0; y<4; y++)
      val[x][y]=getHeight(posId+Eigen::Matrix<int,2,1>(x-1,y-1));
  //gradient
  Vec3T g;
  g[0]=interp1DCubic<T>(pos[1],[&](int yId) {
    return interp1DCubicDiff<T>(pos[0],[&](int xId) {
      return val[xId][yId];
    });
  })*_invCellSzX;
  g[1]=interp1DCubicDiff<T>(pos[1],[&](int yId) {
    return interp1DCubic<T>(pos[0],[&](int xId) {
      return val[xId][yId];
    });
  })*_invCellSzY;
  g[2]=1;
  T len=g.norm(),gradLen=0;
  //hessian
  if(h) {
    Mat2T hessian;
    hessian(0,0)=interp1DCubic<T>(pos[1],[&](int yId) {
      return interp1DCubicDDiff<T>(pos[0],[&](int xId) {
        return val[xId][yId];
      });
    })*_invCellSzX*_invCellSzX;
    hessian(1,1)=interp1DCubicDDiff<T>(pos[1],[&](int yId) {
      return interp1DCubic<T>(pos[0],[&](int xId) {
        return val[xId][yId];
      });
    })*_invCellSzY*_invCellSzY;
    hessian(0,1)=hessian(1,0)=interp1DCubicDiff<T>(pos[1],[&](int yId) {
      return interp1DCubicDiff<T>(pos[0],[&](int xId) {
        return val[xId][yId];
      });
    })*_invCellSzX*_invCellSzY;
    //fill
    h->setZero();
    for(int d=0; d<2; d++) {
      gradLen=2*(g[0]*hessian(0,d)+g[1]*hessian(1,d));
      h->col(d)=Vec3T(-hessian(0,d),-hessian(1,d),0)/len;
      h->col(d)-=Vec3T(-g[0],-g[1],1)*(0.5*gradLen/len/len/len);
    }
  }
  return Vec3T(-g[0],-g[1],1)/len;
}
template <typename T>
const std::vector<Node<int,BBoxExact>>& EnvironmentHeight<T>::getBVH() const {
  return _bvh;
}
template <typename T>
const BBoxExact& EnvironmentHeight<T>::getBB() const {
  return _bb;
}
template <typename T>
bool EnvironmentHeight<T>::empty() const {
  return _height.empty();
}
template <typename T>
void EnvironmentHeight<T>::createHills(double xMin,double xMax,double yMin,double yMax,std::function<double(double,double)> h,double res) {
#define GI(X,Y) (X)+(Y)*(_resX+1)
  _bb=BBoxExact();
  _resX=ceil((xMax-xMin)/res);
  _resY=ceil((yMax-yMin)/res);
  _height.resize((_resX+1)*(_resY+1));
  for(int c=0; c<=_resY; c++)
    for(int r=0; r<=_resX; r++) {
      double xx=xMin+(xMax-xMin)*r/_resX;
      double yy=yMin+(yMax-yMin)*c/_resY;
      _height[GI(r,c)]=h(xx,yy);
      _bb.setUnion(Vec3T(xx,yy,_height[GI(r,c)]).template cast<typename BBoxExact::T>());
    }
  computeInvCellSize();
  buildBVH();
#undef GI
}
template <typename T>
typename EnvironmentHeight<T>::MatT EnvironmentHeight<T>::getHeightMatrix() const {
  return MatTCM(&_height[0],_resX+1,_resY+1,Eigen::OuterStride<>(_resX+1));
}
//helper
template <typename T>
void EnvironmentHeight<T>::computeInvCellSize() {
  Vec2T ext=(getBB().maxCorner()-getBB().minCorner()).template segment<2>(0).template cast<T>();
  _invCellSzX=_resX/ext[0];
  _invCellSzY=_resY/ext[1];
}
template <typename T>
T EnvironmentHeight<T>::getHeight(const Eigen::Matrix<int,2,1>& id) const {
#define GI(X,Y) (X)+(Y)*(_resX+1)
  int x=std::min(std::max(0,id[0]),_resX);
  int y=std::min(std::max(0,id[1]),_resY);
  return _height[GI(x,y)];
#undef GI
}
template <typename T>
Eigen::Matrix<int,2,1> EnvironmentHeight<T>::getPosId(const Vec3T& pt,Vec2T& pos) const {
  pos[0]=(pt[0]-(T)_bb.minCorner()[0])*_invCellSzX;
  pos[1]=(pt[1]-(T)_bb.minCorner()[1])*_invCellSzY;
  Eigen::Matrix<int,2,1> posId(floor(pos[0]),floor(pos[1]));
  posId=posId.cwiseMax(Eigen::Matrix<int,2,1>(1,1)).cwiseMin(Eigen::Matrix<int,2,1>(_resX-2,_resY-2));
  pos-=posId.template cast<T>();
  return posId;
}
template <typename T>
int EnvironmentHeight<T>::buildBVH(int xFrom,int xTo,int yFrom,int yTo) {
#define GI(X,Y) (X)+(Y)*_resX
  Node<int,BBoxExact> n;
  if(xTo-xFrom==1 && yTo-yFrom==1)
    return GI(xFrom,yFrom);
  else if(xTo-xFrom>yTo-yFrom) {
    n._l=buildBVH(xFrom,(xFrom+xTo)/2,yFrom,yTo);
    n._r=buildBVH((xFrom+xTo)/2,xTo,yFrom,yTo);
  } else {
    n._l=buildBVH(xFrom,xTo,yFrom,(yFrom+yTo)/2);
    n._r=buildBVH(xFrom,xTo,(yFrom+yTo)/2,yTo);
  }
  //node
  _bvh[n._l]._parent=(int)_bvh.size();
  _bvh[n._r]._parent=(int)_bvh.size();
  n._nrCell=_bvh[n._l]._nrCell+_bvh[n._r]._nrCell;
  n._bb=_bvh[n._l]._bb;
  n._bb.setUnion(_bvh[n._r]._bb);
  _bvh.push_back(n);
  return (int)_bvh.size()-1;
#undef GI
}
template <typename T>
BBoxExact EnvironmentHeight<T>::BBCell(int X,int Y) const {
  T cellSzX=1/_invCellSzX;
  T cellSzY=1/_invCellSzY;
  //XY bound of cell
  BBoxExact bb;
  bb.minCorner()[0]=_bb.minCorner()[0]+BBoxExact::T(cellSzX*X);
  bb.maxCorner()[0]=_bb.minCorner()[0]+BBoxExact::T(cellSzX*(X+1));
  bb.minCorner()[1]=_bb.minCorner()[1]+BBoxExact::T(cellSzY*Y);
  bb.maxCorner()[1]=_bb.minCorner()[1]+BBoxExact::T(cellSzY*(Y+1));
  //find maximum polynomial value as Z bound of cell
  T val[4][4];
  Eigen::Matrix<int,2,1> posId(X,Y);
  for(int x=0; x<4; x++)
    for(int y=0; y<4; y++)
      val[x][y]=getHeight(posId+Eigen::Matrix<int,2,1>(x-1,y-1));
  toBezier(val);
  bb.minCorner()[2]=_bb.minCorner()[2]-1e6;
  bb.maxCorner()[2]=_bb.minCorner()[2];
  for(int x=0; x<4; x++)
    for(int y=0; y<4; y++)
      bb.maxCorner()[2]=std::max<BBoxExact::T>(bb.maxCorner()[2],BBoxExact::T(val[x][y]));
  return bb;
}
template <typename T>
void EnvironmentHeight<T>::toBezier(T val[4][4]) const {
  /*T sx=T(rand())/T(RAND_MAX);
  T sy=T(rand())/T(RAND_MAX);
  T ref=interp1DCubic<T>(sy,[&](int yId) {
    return interp1DCubic<T>(sx,[&](int xId) {
      return val[xId][yId];
    });
  });*/
  for(int x=0; x<4; x++) {
    T a=val[x][0];
    T b=val[x][1];
    T c=val[x][2];
    T d=val[x][3];
    val[x][0]=b;
    val[x][1]=-a/6+b+c/6;
    val[x][2]= b/6+c-d/6;
    val[x][3]=c;
  }
  for(int y=0; y<4; y++) {
    T a=val[0][y];
    T b=val[1][y];
    T c=val[2][y];
    T d=val[3][y];
    val[0][y]=b;
    val[1][y]=-a/6+b+c/6;
    val[2][y]= b/6+c-d/6;
    val[3][y]=c;
  }
  /*T bezier=interp1DCubicBezier<T>(sy,[&](int yId) {
    return interp1DCubicBezier<T>(sx,[&](int xId) {
      return val[xId][yId];
    });
  });
  std::cout << "ref=" << ref << " bezier=" <<bezier << std::endl;*/
}
template <typename T>
void EnvironmentHeight<T>::buildBVH() {
#define GI(X,Y) (X)+(Y)*_resX
  Node<int,BBoxExact> n;
  n._nrCell=1;
  _bvh.resize(_resX*_resY);
  for(int x=0; x<_resX; x++)
    for(int y=0; y<_resY; y++) {
      n._cell=GI(x,y);
      n._bb=BBCell(x,y);
      _bvh[GI(x,y)]=n;
    }
  buildBVH(0,_resX,0,_resY);
#undef GI
}
template class EnvironmentExact<FLOAT>;
template class EnvironmentHeight<FLOAT>;
#ifdef FORCE_ADD_DOUBLE_PRECISION
template class EnvironmentExact<double>;
template class EnvironmentHeight<double>;
#endif
}
