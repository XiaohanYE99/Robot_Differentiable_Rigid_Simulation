#include "DeformedEnvironment.h"
#include <Utils/CrossSpatialUtils.h>
#include <Utils/VTKWriter.h>
#include <Utils/Interp.h>
#include <Utils/IO.h>

namespace PHYSICSMOTION {
template <typename T>
T& DeformedEnvironment<T>::NodeFetcher::operator[](int xyz) {
  ASSERT(xyz>=0 && xyz<64)
  int xy=xyz%16;
  Eigen::Matrix<int,3,1> off(xy%4,xy/4,xyz/16);
  return _nodes[(_offset+_stride.dot(off))*3+_DIM];
}
template <typename T>
const T& DeformedEnvironment<T>::NodeFetcher::operator[](int xyz) const {
  return const_cast<NodeFetcher&>(*this)[xyz];
}
template <typename T>
DeformedEnvironment<T>::DeformedEnvironment():_historyPath("history") {}
template <typename T>
DeformedEnvironment<T>::DeformedEnvironment(const Eigen::Matrix<int,3,1>& nrCell,std::shared_ptr<Environment<T>> env):_historyPath("history") {
  _env=env;
  Eigen::Matrix<int,3,1> sz=nrCell+Eigen::Matrix<int,3,1>::Constant(3);
  _nodes.resize(sz.prod()*3);
  _Dnodes.resize(sz.prod()*3);
  _nrCell=nrCell;
  _stride=Eigen::Matrix<int,3,1>(1,sz.x(),sz.x()*sz.y());
  //trips
  STrips trips;
  int nrBoundarySample=1,row=0;
  for(int x=0; x<=_nrCell[0]*nrBoundarySample; x++) {
    jacobian(row++,Vec3T(T(x)/nrBoundarySample,0,0),trips,1);
    jacobian(row++,Vec3T(T(x)/nrBoundarySample,_nrCell[1],0),trips,1);
  }
  for(int y=0; y<=_nrCell[1]*nrBoundarySample; y++) {
    jacobian(row++,Vec3T(0,T(y)/nrBoundarySample,0),trips,0);
    jacobian(row++,Vec3T(_nrCell[0],T(y)/nrBoundarySample,0),trips,0);
  }
  _fixBasis.resize(row,_nodes.size());
  _fixBasis.setFromTriplets(trips.begin(),trips.end());
}
template <typename T>
bool DeformedEnvironment<T>::read(std::istream& is,IOData* dat) {
  registerType<EnvironmentExact<T>>(dat);
  registerType<EnvironmentHeight<T>>(dat);
  readBinaryData(_env,is,dat);
  readBinaryData(_nrCell,is);
  readBinaryData(_stride,is);
  readBinaryData(_nodes,is);
  readBinaryData(_Dnodes,is);
  return is.good();
}
template <typename T>
bool DeformedEnvironment<T>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_env,os,dat);
  writeBinaryData(_nrCell,os);
  writeBinaryData(_stride,os);
  writeBinaryData(_nodes,os);
  writeBinaryData(_Dnodes,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> DeformedEnvironment<T>::copy() const {
  return std::shared_ptr<SerializableBase>(new DeformedEnvironment<T>());
}
template <typename T>
std::string DeformedEnvironment<T>::type() const {
  return typeid(DeformedEnvironment<T>).name();
}
template <typename T>
void DeformedEnvironment<T>::setHistoryPath(const std::string& path) {
  _historyPath=path+"/history";
}
template <typename T>
typename std::vector<typename DeformedEnvironment<T>::Vec3T>
DeformedEnvironment<T>::getIsoline(int dir,T x,T y) const {
  Mat3T J;
  Vec3T alpha;
  std::vector<Vec3T> ret;
  alpha[(dir+1)%3]=x;
  alpha[(dir+2)%3]=y;
  for(int c=0; c<=_nrCell[dir]; c++) {
    alpha[dir]=c;
    ret.push_back(forward(alpha,&J,NULL));
    ret.push_back(J.col(dir));
  }
  return ret;
}
template <typename T>
void DeformedEnvironment<T>::writeVTK(const std::string& path,int res,bool bottomOnly) const {
  std::vector<Eigen::Matrix<int,4,1>> iss;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  //x
  if(!bottomOnly)
    for(int x=0; x<=_nrCell[0]; x++) {
      int stride=_nrCell[2]*res+1;
      int offset=(int)vss.size();
      for(int y=0; y<=_nrCell[1]*res; y++)
        for(int z=0; z<=_nrCell[2]*res; z++) {
          Vec3T a(x,T(y)/res,T(z)/res);
          vss.push_back(Vec3T(forward(a,NULL,NULL)).template cast<double>());
          if(y<_nrCell[1]*res && z<_nrCell[2]*res) {
            int off=offset+y*stride+z;
            iss.push_back(Eigen::Matrix<int,4,1>(off,off+1,off+stride+1,off+stride));
          }
        }
    }
  //y
  if(!bottomOnly)
    for(int y=0; y<=_nrCell[1]; y++) {
      int stride=_nrCell[2]*res+1;
      int offset=(int)vss.size();
      for(int x=0; x<=_nrCell[0]*res; x++)
        for(int z=0; z<=_nrCell[2]*res; z++) {
          Vec3T a(T(x)/res,y,T(z)/res);
          vss.push_back(Vec3T(forward(a,NULL,NULL)).template cast<double>());
          if(x<_nrCell[0]*res && z<_nrCell[2]*res) {
            int off=offset+x*stride+z;
            iss.push_back(Eigen::Matrix<int,4,1>(off,off+1,off+stride+1,off+stride));
          }
        }
    }
  //z
  for(int z=0; z<=(bottomOnly?0:_nrCell[2]); z++) {
    int stride=_nrCell[1]*res+1;
    int offset=(int)vss.size();
    for(int x=0; x<=_nrCell[0]*res; x++)
      for(int y=0; y<=_nrCell[1]*res; y++) {
        Vec3T a(T(x)/res,T(y)/res,z);
        vss.push_back(Vec3T(forward(a,NULL,NULL)).template cast<double>());
        if(x<_nrCell[0]*res && y<_nrCell[1]*res) {
          int off=offset+x*stride+y;
          iss.push_back(Eigen::Matrix<int,4,1>(off,off+1,off+stride+1,off+stride));
        }
      }
  }
  //output
  VTKWriter<double> os("env",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::QUAD);
}
template <typename T>
void DeformedEnvironment<T>::writeFrameVTK(const std::string& path,T len,int res) const {
  std::vector<Eigen::Matrix<double,4,1>> css;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  for(int x=0; x<=_nrCell[0]*res; x++)
    for(int y=0; y<=_nrCell[1]*res; y++)
      for(int z=0; z<=_nrCell[2]*res; z++) {
        Vec3T a(T(x)/res,T(y)/res,T(z)/res);
        Vec3T pos=forward(a,NULL,NULL);
        Mat3T R=rotation(a,NULL);
        //vss
        vss.push_back(pos.template cast<double>());
        vss.push_back((pos+R.col(0)*len).template cast<double>());
        vss.push_back(pos.template cast<double>());
        vss.push_back((pos+R.col(1)*len).template cast<double>());
        vss.push_back(pos.template cast<double>());
        vss.push_back((pos+R.col(2)*len).template cast<double>());
        //css
        css.push_back(Eigen::Matrix<double,4,1>(1,0,0,1));
        css.push_back(Eigen::Matrix<double,4,1>(1,0,0,1));
        css.push_back(Eigen::Matrix<double,4,1>(0,1,0,1));
        css.push_back(Eigen::Matrix<double,4,1>(0,1,0,1));
        css.push_back(Eigen::Matrix<double,4,1>(0,0,1,1));
        css.push_back(Eigen::Matrix<double,4,1>(0,0,1,1));
      }
  //output
  VTKWriter<double> os("env",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCustomPointColorData("color",css.begin(),css.end());
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>(0,2,0),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>((int)vss.size()/2,2,0),
                 VTKWriter<double>::LINE);
}
template <typename T>
void DeformedEnvironment<T>::optimizeLocalGlobal(int heightRes,int volumeRes,T wACAP,T tol,bool writeHistory) {
  initializeNodes();
  std::vector<Eigen::Matrix<double,3,1>> vss;
  //we optimize two kinds of energies: height energy and ACAP energy
  //  \sum_i|pos_i(n)-pos_i^*|^2/2+\sum_j|J_j(n)-sR_j|^2/2
  //the optimizer is local-global type:
  //  local: fixing n and update pos_i^*, sR_j
  //  global: fixing pos_i^*, sR_j and update n
  STrips trips;
  wACAP=sqrt(wACAP);
  std::vector<T> RHS;
  OMP_PARALLEL_FOR_
  for(int x=0; x<=_nrCell[0]*heightRes; x++)
    for(int y=0; y<=_nrCell[1]*heightRes; y++) {
      int row=(x*(_nrCell[1]*heightRes+1)+y)*3;
      jacobian(row,Vec3T(T(x)/heightRes,T(y)/heightRes,0),trips);
    }
  OMP_PARALLEL_FOR_
  for(int x=0; x<=_nrCell[0]*volumeRes; x++)
    for(int y=0; y<=_nrCell[1]*volumeRes; y++)
      for(int z=0; z<=_nrCell[2]*volumeRes; z++) {
        int row=(_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1)*3;
        row+=(x*(_nrCell[1]*volumeRes+1)*(_nrCell[2]*volumeRes+1)+y*(_nrCell[2]*volumeRes+1)+z)*9;
        jacobianJ(row,Vec3T(T(x)/volumeRes,T(y)/volumeRes,T(z)/volumeRes),trips,wACAP);
      }
  //assemble and handle fixed nodes
  //invJTJ
  SMatT J,JTJ,JTJKKT;
  clock_t time=clock();
  J.resize((_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1)*3+
           (_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1)*(_nrCell[2]*heightRes+1)*9,_nodes.size());
  J.setFromTriplets(trips.begin(),trips.end());
  JTJ=J.transpose()*J;
  JTJKKT=buildKKT<T,0,int>(JTJ,_fixBasis,0);
  Eigen::SparseLU<SMatT> invJTJ(JTJKKT);
  ASSERT(invJTJ.info()==Eigen::Success)
  //main loop
  if(writeHistory)
    recreate(_historyPath);
  for(int iter=0;; iter++) {
    //local step
    RHS.resize(J.rows());
    vss.resize((_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1)*2);
    //height
    OMP_PARALLEL_FOR_
    for(int x=0; x<=_nrCell[0]*heightRes; x++)
      for(int y=0; y<=_nrCell[1]*heightRes; y++) {
        Vec3T pos=forward(Vec3T(T(x)/heightRes,T(y)/heightRes,0),NULL,NULL);
        int row=x*(_nrCell[1]*heightRes+1)+y;
        if(writeHistory)
          vss[row*2+0]=pos.template cast<double>();
        pos-=_env->phiGrad(pos,NULL)*_env->phi(pos,NULL);
        if(writeHistory)
          vss[row*2+1]=pos.template cast<double>();
        Eigen::Map<Vec3T>(RHS.data()+row*3)=pos;
      }
    if(writeHistory) {
      VTKWriter<double> os("env",_historyPath+"/iterP"+std::to_string(iter)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>(0,2,0),
                     VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,3,1>>(vss.size()/2,2,0),
                     VTKWriter<double>::LINE);
    }
    //ACAP
    OMP_PARALLEL_FOR_
    for(int x=0; x<=_nrCell[0]*volumeRes; x++)
      for(int y=0; y<=_nrCell[1]*volumeRes; y++)
        for(int z=0; z<=_nrCell[2]*volumeRes; z++) {
          Mat3T Jpos;
          forward(Vec3T(T(x)/volumeRes,T(y)/volumeRes,T(z)/volumeRes),&Jpos,NULL);
          int row=(_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1)*3;
          row+=(x*(_nrCell[1]*volumeRes+1)*(_nrCell[2]*volumeRes+1)+y*(_nrCell[2]*volumeRes+1)+z)*9;
          Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> SVD(Jpos.template cast<double>(),Eigen::ComputeFullU|Eigen::ComputeFullV);
          Mat3T R=(SVD.matrixU()*SVD.matrixV().transpose()).template cast<T>();
          R*=pow(Jpos.determinant(),T(1/3.0))*wACAP;
          Eigen::Map<Mat3T>(RHS.data()+row)=R;
        }
    //global step
    //we solve JTJ*(n+dn)=b and fixing some nodes
    //we have JTJ*dn=b-JTJ*n and dn=invJTJ*(b-JTJ*n)
    Vec b=J.transpose()*Eigen::Map<Vec>(RHS.data(),RHS.size())-JTJ*_nodes;
    _Dnodes=invJTJ.solve(concat<Vec>(b,Vec::Zero(_fixBasis.rows()))).segment(0,_nodes.size());
    T E=localInjectiveEnergy(heightRes,volumeRes,wACAP,0,NULL);
    std::cout << "iter=" << iter << " time=" << double(clock()-time)/CLOCKS_PER_SEC << " E=" << E <<
              " deltaNorm=" << _Dnodes.norm() << " deltaNormMax=" << _Dnodes.cwiseAbs().maxCoeff() << std::endl;
    if(writeHistory)
      writeVTK(_historyPath+"/iter"+std::to_string(iter)+".vtk");
    if(_Dnodes.cwiseAbs().maxCoeff()<tol)
      break;
    _nodes+=_Dnodes;
  }
}
template <typename T>
void DeformedEnvironment<T>::optimizeLocalInjective(int heightRes,int volumeRes,T wACAP,T mu,T tol,bool writeHistory) {
  initializeNodes();
  //main loop
  SMatT h;
  Vec HDnodes,nodesPrev;
  T alpha=1,alphaMin=1e-10,alphaMax=1e3,EPrev=0;
  if(writeHistory)
    recreate(_historyPath);
  clock_t time=clock();
  Eigen::SparseLU<SMatT> invH;
  Eigen::SimplicialLDLT<SMatT> invFBTFB(_fixBasis*_fixBasis.transpose());
  for(int iter=0; alpha>alphaMin; iter++) {
#define NEWTON
#ifdef NEWTON
    T E=localInjectiveEnergy(heightRes,volumeRes,wACAP,mu,&_Dnodes,&h);
    h=buildKKT<T,0,int>(h,_fixBasis,0);
    invH.compute(h);
    if(invH.info()!=Eigen::Success)
      break;
    HDnodes=invH.solve(concat<Vec>(_Dnodes,Vec::Zero(_fixBasis.rows()))).segment(0,_nodes.size());
    _Dnodes-=_fixBasis.transpose()*invFBTFB.solve(_fixBasis*_Dnodes);
#else
    T E=localInjectiveEnergy(heightRes,volumeRes,wACAP,mu,&_Dnodes,NULL);
    _Dnodes-=_fixBasis.transpose()*invFBTFB.solve(_fixBasis*_Dnodes);
    HDnodes=_Dnodes;
#endif
    //writeVTK
    if(writeHistory)
      writeVTK(_historyPath+"/iter"+std::to_string(iter)+".vtk");
    //terminate
    std::cout << "iter=" << iter << " time=" << double(clock()-time)/CLOCKS_PER_SEC <<
              " E=" << E << " alpha=" << alpha << " gNorm=" << _Dnodes.norm() <<
              " gNormMax=" << _Dnodes.cwiseAbs().maxCoeff() << std::endl;
    if(iter>0 && EPrev-E<tol)
      break;
    //line search
    EPrev=E;
    nodesPrev=_nodes;
    while(alpha>alphaMin) {
      _nodes=nodesPrev-HDnodes*alpha;
      T E2=localInjectiveEnergy(heightRes,volumeRes,wACAP,mu,NULL);
      if(isfinite(E2) && E2<E) {
        alpha=std::min<T>(alphaMax,alpha*1.1);
        break;
      } else {
        alpha*=0.5;
      }
    }
  }
}
template <typename T>
void DeformedEnvironment<T>::initializeNodes() {
  Eigen::Matrix<int,3,1> nrSample=_nrCell*3;
  STrips trips;
  std::vector<T> pos;
  BBoxExact bb=_env->getBB();
  Vec3T minC=bb.minCorner().template cast<T>();
  Vec3T maxC=bb.maxCorner().template cast<T>();
  T maxCX=minC[2]+(maxC[0]-minC[0])*_nrCell[2]/_nrCell[0];
  T maxCY=minC[2]+(maxC[1]-minC[1])*_nrCell[2]/_nrCell[1];
  maxC[2]=std::max(maxCX,maxCY);
  for(int x=0; x<=nrSample[0]; x++)
    for(int y=0; y<=nrSample[1]; y++)
      for(int z=0; z<=nrSample[2]; z++) {
        //alpha
        Vec3T a(x,y,z);
        a.array()*=_nrCell.array().template cast<T>();
        a.array()/=nrSample.array().template cast<T>();
        jacobian((int)pos.size(),a,trips);
        //pos
        Vec3T p;
        for(int d=0; d<3; d++)
          p[d]=interp1D<T,T>(minC[d],maxC[d],a[d]/_nrCell[d]);
        pos.push_back(p[0]);
        pos.push_back(p[1]);
        pos.push_back(p[2]);
      }
  //solve
  SMatT J;
  J.resize((int)pos.size(),_nodes.size());
  J.setFromTriplets(trips.begin(),trips.end());
  Vec RHS=J.transpose()*Eigen::Map<Vec>(pos.data(),pos.size());
  _nodes=Eigen::SimplicialLLT<SMatT>(J.transpose()*J).solve(RHS);
}
template <typename T>
void DeformedEnvironment<T>::debug(int N) {
  DEFINE_NUMERIC_DELTA_T(T)
  _nodes.setRandom();

  //local injective energy
  Vec g;
  _Dnodes.setRandom(_nodes.size());
  T E=localInjectiveEnergy(4,4,0.5,0,&g);
  _nodes+=_Dnodes*DELTA;
  T E2=localInjectiveEnergy(4,4,0.5,0,NULL);
  DEBUG_GRADIENT("localInjectiveEnergy",_Dnodes.dot(g),_Dnodes.dot(g)-(E2-E)/DELTA)

  for(int iter=0; iter<N; iter++) {
    Vec3T a=(Vec3T::Random()+Vec3T::Ones())*0.5;
    a.array()*=_nrCell.template cast<T>().array();
    Vec3T da=Vec3T::Random();
    //rotation
    Mat3T R,R2,w;
    R=rotation(a,&w);
    R2=rotation(a+da*DELTA,NULL);
    DEBUG_GRADIENT("rotation",(cross<T>(w*da)*R).norm(),(cross<T>(w*da)*R-(R2-R)/DELTA).norm())
  }

  for(int iter=0; iter<N; iter++) {
    Vec3T a=(Vec3T::Random()+Vec3T::Ones())*0.5;
    a.array()*=_nrCell.template cast<T>().array();
    Vec3T da=Vec3T::Random();
    //forward hessian vector
    Mat3T gradV,gradV2,hessV[3];
    forward(a,&gradV,hessV);
    forward(a+da*DELTA,&gradV2,NULL);
    DEBUG_GRADIENT("forward-hessV",(hessV[0]*da[0]+hessV[1]*da[1]+hessV[2]*da[2]).norm(),
                   (hessV[0]*da[0]+hessV[1]*da[1]+hessV[2]*da[2]-(gradV2-gradV)/DELTA).norm())
  }

  for(int iter=0; iter<N; iter++) {
    Vec3T a=(Vec3T::Random()+Vec3T::Ones())*0.5;
    a.array()*=_nrCell.template cast<T>().array();
    Vec3T da=Vec3T::Random();
    _nodes.setRandom();

    Mat3T hess;
    Vec3T grad,grad2;
    Eigen::Matrix<int,3,1> aFloor=aFloorSafe(a);
    int offset=_stride.dot(aFloor);
    for(int d=0; d<3; d++) {
      //forward
      NodeFetcher fetcher({_stride,offset,d,const_cast<Vec&>(_nodes)});
      T val=forward(a-aFloor.cast<T>(),fetcher,&grad,&hess);
      T val2=forward(a-aFloor.cast<T>()+da*DELTA,fetcher,&grad2,NULL);
      DEBUG_GRADIENT("forward-grad",da.dot(grad),da.dot(grad)-(val2-val)/DELTA)
      DEBUG_GRADIENT("forward-hess",(hess*da).norm(),(hess*da-(grad2-grad)/DELTA).norm())
      //backward
      _Dnodes.setZero();
      Vec DN=Vec::Random(_nodes.size());
      NodeFetcher Dfetcher({_stride,offset,d,const_cast<Vec&>(_Dnodes)});
      backward(a-aFloor.cast<T>(),Dfetcher,1,Vec3T::Zero());
      _nodes+=DN*DELTA;
      val2=forward(a-aFloor.cast<T>(),fetcher,NULL,NULL);
      DEBUG_GRADIENT("backward",_Dnodes.dot(DN),_Dnodes.dot(DN)-(val2-val)/DELTA)
    }
  }

  for(int iter=0; iter<N; iter++) {
    Vec3T a=(Vec3T::Random()+Vec3T::Ones())*0.5;
    a.array()*=_nrCell.template cast<T>().array();
    _nodes.setRandom();

    //jacobian
    Mat3T J;
    Vec3T pos=forward(a,&J,NULL);
    SMatT Jacobian=jacobian(a);
    SMatT JacobianJ=jacobianJ(a);
    DEBUG_GRADIENT("jacobian",(Jacobian*_nodes).norm(),(Jacobian*_nodes-pos).norm())
    DEBUG_GRADIENT("jacobianJ",(JacobianJ*_nodes).norm(),(JacobianJ*_nodes-Eigen::Map<Vec9T>(J.data())).norm())
  }
}
template <typename T>
typename DeformedEnvironment<T>::Vec3T DeformedEnvironment<T>::nrCell() const {
  return _nrCell.template cast<T>();
}
template <typename T>
std::shared_ptr<Environment<T>> DeformedEnvironment<T>::getOptimizedEnvironment(int res) const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  int stride=_nrCell[1]*res+1;
  for(int x=0; x<=_nrCell[0]*res; x++)
    for(int y=0; y<=_nrCell[1]*res; y++) {
      Vec3T a(T(x)/res,T(y)/res,0);
      vss.push_back(Vec3T(forward(a,NULL,NULL)).template cast<double>());
      if(x<_nrCell[0]*res && y<_nrCell[1]*res) {
        int off=x*stride+y;
        iss.push_back(Eigen::Matrix<int,3,1>(off,off+stride+1,off+1));
        iss.push_back(Eigen::Matrix<int,3,1>(off,off+stride,off+stride+1));
      }
    }
  MeshExact mesh(vss,iss);
  return std::shared_ptr<EnvironmentExact<T>>(new EnvironmentExact<T>(mesh));
}
template <typename T>
T DeformedEnvironment<T>::localInjectiveEnergy(int heightRes,int volumeRes,T wACAP,T mu,Vec* g,SMatT* h) {
  Mat3T H;
  STrips trips;
  H << 2./3.,-1./3.,-1./3.,
  -1./3., 2./3.,-1./3.,
  -1./3.,-1./3., 2./3.;
  Mat3T HTH=H.transpose()*H*wACAP;
  ParallelMatrix<Vec> EP(Vec::Zero(1));
  ParallelMatrix<Vec> GP(Vec::Zero(_nodes.size()));
  OMP_PARALLEL_FOR_
  for(int x=0; x<=_nrCell[0]*heightRes; x++)
    for(int y=0; y<=_nrCell[1]*heightRes; y++) {
      //std::cout << "heightEnergy(" << x << "," << y << ")" << std::endl;
      Vec3T a(T(x)/heightRes,T(y)/heightRes,0);
      Vec3T pos=forward(a,NULL,NULL),gPhi;
      T phi=_env->phi(pos,g?&gPhi:NULL);
      EP.getMatrixI()[0]+=phi*phi/2;
      SMatT JPhi;
      if(g || h)
        JPhi=jacobian(a).transpose()*toSparse<T,0,int,Vec3T>(gPhi);
      if(g)
        GP.getMatrixI()+=JPhi*phi;
      if(h) {
        int row=x*(_nrCell[1]*heightRes+1)+y;
        addBlock<T,0,int>(trips,row,0,JPhi.transpose());
      }
    }
  OMP_PARALLEL_FOR_
  for(int x=0; x<=_nrCell[0]*volumeRes; x++)
    for(int y=0; y<=_nrCell[1]*volumeRes; y++)
      for(int z=0; z<=_nrCell[2]*volumeRes; z++) {
        //std::cout << "ACAPEnergy(" << x << "," << y << "," << z << ")" << std::endl;
        Vec3T a(T(x)/volumeRes,T(y)/volumeRes,T(z)/volumeRes);
        Mat3T Jpos;
        SMatT JJac;
        forward(a,&Jpos,NULL);
        Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> SVD(Jpos.template cast<double>(),Eigen::ComputeFullU|Eigen::ComputeFullV);
        if(g || h)
          JJac=jacobianJ(a).transpose()*toSparse<T,0,int,MatT>(DSDJ(SVD)).transpose();
        if(g) {
          Vec3T DEDSV=HTH*SVD.singularValues().template cast<T>();
          if(mu>0)
            DEDSV-=mu*(1/SVD.singularValues().array()).matrix().template cast<T>();
          GP.getMatrixI()+=JJac*DEDSV;
        }
        if(h) {
          int row=(_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1);
          row+=(x*(_nrCell[1]*volumeRes+1)*(_nrCell[2]*volumeRes+1)+y*(_nrCell[2]*volumeRes+1)+z)*6;
          addBlock<T,0,int>(trips,row,0,toSparse<T,0,int,Mat3T>(H)*JJac.transpose()*sqrt(wACAP));
          if(mu>0) {
            Vec3T DDEDDSV=sqrt(mu)*(1/SVD.singularValues().array()).matrix().template cast<T>();
            addBlock<T,0,int>(trips,row+3,0,toSparse<T,0,int,Mat3T>(DDEDDSV.asDiagonal())*JJac.transpose());
          }
        }
        EP.getMatrixI()[0]+=SVD.singularValues().template cast<T>().dot(HTH*SVD.singularValues().template cast<T>())/2;
        if(mu>0) {
          EP.getMatrixI()[0]-=mu*log(SVD.singularValues()[0]);
          EP.getMatrixI()[0]-=mu*log(SVD.singularValues()[1]);
          EP.getMatrixI()[0]-=mu*log(SVD.singularValues()[2]);
        }
      }
  if(g)
    *g=GP.getMatrix();
  if(h) {
    SMatT J;
    int row=(_nrCell[0]*heightRes+1)*(_nrCell[1]*heightRes+1);
    row+=(_nrCell[0]*volumeRes+1)*(_nrCell[1]*volumeRes+1)*(_nrCell[2]*volumeRes+1)*6;
    J.resize(row,_nodes.size());
    J.setFromTriplets(trips.begin(),trips.end());
    *h=J.transpose()*J;
  }
  return EP.getMatrix()[0];
}
template <typename T>
typename DeformedEnvironment<T>::Mat3T DeformedEnvironment<T>::rotation(const Vec3T& alpha,Mat3T* w,Vec3T* posOut,Mat3T* JposOut) const {
  Mat3T Jpos,Hpos[3];
  Vec3T pos=forward(alpha,&Jpos,w?Hpos:NULL);
  if(posOut)
    *posOut=pos;
  if(JposOut)
    *JposOut=Jpos;
  Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> SVD(Jpos.template cast<double>(),Eigen::ComputeFullU|Eigen::ComputeFullV);
  Mat3T R=(SVD.matrixU()*SVD.matrixV().transpose()).template cast<T>();
  //T detR=R.determinant();
  if(R.determinant()<0)
    R*=-1;
  if(w) {
    Eigen::Matrix<double,3,1> invSigma=(1/(SVD.singularValues().sum()-SVD.singularValues().array())).matrix();
    for(int d=0; d<3; d++) {
      Eigen::Matrix<double,3,3> RTA=(R.transpose()*Hpos[d]).template cast<double>();
      w->col(d)=(SVD.matrixU()*(invSigma.asDiagonal()*(SVD.matrixV().transpose()*invCross<double>(RTA-RTA.transpose())))).template cast<T>();
      //if(detR<0)
      //  w->col(d)*=-1;
    }
  }
  return R;
}
template <typename T>
Eigen::Matrix<int,3,1> DeformedEnvironment<T>::aFloorSafe(const Vec3T& a) const {
  return a.template cast<int>().cwiseMin((_nrCell.array()-1).matrix()).cwiseMax(Eigen::Matrix<int,3,1>::Zero());
}
template <typename T>
void DeformedEnvironment<T>::jacobian(int row,const Vec3T& a,STrips& trips,int DIM) const {
  Vec ret=Vec::Zero(64*3);
  for(int d=0; d<3; d++)
    if(DIM==-1 || d==DIM) {
      Eigen::Matrix<int,3,1> aFloor=aFloorSafe(a);
      NodeFetcher fetcher({Eigen::Matrix<int,3,1>(1,4,16),0,d,ret});
      backward(a-aFloor.cast<T>(),fetcher,1,Vec3T::Zero());
      //assign
      int offset=_stride.dot(aFloor);
      for(int xyz=0; xyz<64; xyz++) {
        int xy=xyz%16;
        Eigen::Matrix<int,3,1> off(xy%4,xy/4,xyz/16);
        trips.push_back(STrip(row,(offset+_stride.dot(off))*3+d,ret[xyz*3+d]));
      }
      row++;
    }
}
template <typename T>
void DeformedEnvironment<T>::jacobianJ(int row,const Vec3T& a,STrips& trips,T wACAP) const {
  for(int c=0; c<3; c++) {
    Vec ret=Vec::Zero(64*3);
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,3,1> aFloor=aFloorSafe(a);
      NodeFetcher fetcher({Eigen::Matrix<int,3,1>(1,4,16),0,d,ret});
      backward(a-aFloor.cast<T>(),fetcher,0,Vec3T::Unit(c));
      //assign
      int offset=_stride.dot(aFloor);
      for(int xyz=0; xyz<64; xyz++) {
        int xy=xyz%16;
        Eigen::Matrix<int,3,1> off(xy%4,xy/4,xyz/16);
        trips.push_back(STrip(row+d+c*3,(offset+_stride.dot(off))*3+d,ret[xyz*3+d]*wACAP));
      }
    }
  }
}
template <typename T>
typename DeformedEnvironment<T>::SMatT DeformedEnvironment<T>::jacobian(const Vec3T& a) const {
  SMatT ret;
  STrips trips;
  jacobian(0,a,trips);
  ret.resize(3,_nodes.size());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T>
typename DeformedEnvironment<T>::SMatT DeformedEnvironment<T>::jacobianJ(const Vec3T& a,T wACAP) const {
  SMatT ret;
  STrips trips;
  jacobianJ(0,a,trips,wACAP);
  ret.resize(9,_nodes.size());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T>
typename DeformedEnvironment<T>::Vec3T DeformedEnvironment<T>::alpha(const Vec3T& pos,const T tol) const {
  //initial guess
  Mat3T h,J;
  T bestEnergy=0,energy;
  Vec3T ret,retCurr,g,d,posCurr;
  for(int x=0; x<_nrCell[0]; x++)
    for(int y=0; y<_nrCell[1]; y++)
      for(int z=0; z<_nrCell[2]; z++) {
        Vec3T a(x+0.5,y+0.5,z+0.5);
        energy=(forward(a,NULL,NULL)-pos).squaredNorm()/2;
        if(energy<bestEnergy || (x==0 && y==0 && z==0)) {
          bestEnergy=energy;
          ret=a;
        }
      }
  //solve
  T alpha=1,alphaMin=1e-10,alphaMax=1e3;
  while(true) {
    posCurr=forward(ret,&J,NULL);
    energy=(posCurr-pos).squaredNorm()/2;
    g=J.transpose()*(posCurr-pos);
    if(g.norm()<tol)
      break;
    h=J.transpose()*J;
    d=h.fullPivLu().solve(g);
    retCurr=ret;
    while(alpha>alphaMin) {
      ret=retCurr-d*alpha;
      if((forward(ret,NULL,NULL)-pos).squaredNorm()/2<energy) {
        alpha=std::min<T>(alphaMax,alpha*1.1);
        break;
      } else alpha*=0.5;
    }
    if(alpha<=alphaMin)
      break;
  }
  return ret;
}
template <typename T>
typename DeformedEnvironment<T>::Vec3T DeformedEnvironment<T>::forward(const Vec3T& a,Mat3T* J,Mat3T* H) const {
  Mat3T HRow;
  Vec3T ret,JRow;
  Eigen::Matrix<int,3,1> aFloor=aFloorSafe(a);
  int offset=_stride.dot(aFloor);
  for(int d=0; d<3; d++) {
    NodeFetcher fetcher({_stride,offset,d,const_cast<Vec&>(_nodes)});
    ret[d]=forward(a-aFloor.cast<T>(),fetcher,J?&JRow:NULL,H?&HRow:NULL);
    if(J)
      J->row(d)=JRow;
    if(H)
      for(int d2=0; d2<3; d2++)
        H[d2].row(d)=HRow.col(d2);
  }
  return ret;
}
template <typename T>
void DeformedEnvironment<T>::backward(const Vec3T& a,const Vec3T& dEde,const Mat3T& dEdJ) {
  Eigen::Matrix<int,3,1> aFloor=aFloorSafe(a);
  int offset=_stride.dot(aFloor);
  for(int d=0; d<3; d++) {
    NodeFetcher fetcher({_stride,offset,d,_Dnodes});
    backward(a,fetcher,dEde[d],dEdJ.row(d));
  }
}
template <typename T>
T DeformedEnvironment<T>::forward(const Vec3T& a,NodeFetcher n,Vec3T* grad,Mat3T* hess) {
  T x0=pow(a[2], 3);
  T x1=(1.0/6.0)*x0;
  T x2=1 - a[2];
  T x3=(1.0/6.0)*pow(x2, 3);
  T x4=pow(a[2], 2);
  T x5=(1.0/2.0)*x0;
  T x6=-x4 + x5 + 2.0/3.0;
  T x7=(1.0/2.0)*x4;
  T x8=(1.0/2.0)*a[2] - x5 + x7 + 1.0/6.0;
  T x9=n[15]*x3 + n[31]*x6 + n[47]*x8 + n[63]*x1;
  T x10=pow(a[1], 3);
  T x11=(1.0/6.0)*x10;
  T x12=n[19]*x6 + n[35]*x8 + n[3]*x3 + n[51]*x1;
  T x13=1 - a[1];
  T x14=(1.0/6.0)*pow(x13, 3);
  T x15=pow(a[1], 2);
  T x16=(1.0/2.0)*x10;
  T x17=-x15 + x16 + 2.0/3.0;
  T x18=n[23]*x6 + n[39]*x8 + n[55]*x1 + n[7]*x3;
  T x19=(1.0/2.0)*x15;
  T x20=(1.0/2.0)*a[1] - x16 + x19 + 1.0/6.0;
  T x21=n[11]*x3 + n[27]*x6 + n[43]*x8 + n[59]*x1;
  T x22=x11*x9 + x12*x14 + x17*x18 + x20*x21;
  T x23=pow(a[0], 3);
  T x24=(1.0/6.0)*x23;
  T x25=n[12]*x3 + n[28]*x6 + n[44]*x8 + n[60]*x1;
  T x26=n[0]*x3 + n[16]*x6 + n[32]*x8 + n[48]*x1;
  T x27=n[20]*x6 + n[36]*x8 + n[4]*x3 + n[52]*x1;
  T x28=n[24]*x6 + n[40]*x8 + n[56]*x1 + n[8]*x3;
  T x29=x11*x25 + x14*x26 + x17*x27 + x20*x28;
  T x30=1 - a[0];
  T x31=(1.0/6.0)*pow(x30, 3);
  T x32=pow(a[0], 2);
  T x33=(1.0/2.0)*x23;
  T x34=-x32 + x33 + 2.0/3.0;
  T x35=n[13]*x3 + n[29]*x6 + n[45]*x8 + n[61]*x1;
  T x36=n[17]*x6 + n[1]*x3 + n[33]*x8 + n[49]*x1;
  T x37=n[21]*x6 + n[37]*x8 + n[53]*x1 + n[5]*x3;
  T x38=n[25]*x6 + n[41]*x8 + n[57]*x1 + n[9]*x3;
  T x39=x11*x35 + x14*x36 + x17*x37 + x20*x38;
  T x40=(1.0/2.0)*x32;
  T x41=(1.0/2.0)*a[0] - x33 + x40 + 1.0/6.0;
  T x42=n[14]*x3 + n[30]*x6 + n[46]*x8 + n[62]*x1;
  T x43=n[18]*x6 + n[2]*x3 + n[34]*x8 + n[50]*x1;
  T x44=n[22]*x6 + n[38]*x8 + n[54]*x1 + n[6]*x3;
  T x45=n[10]*x3 + n[26]*x6 + n[42]*x8 + n[58]*x1;
  T x46=x11*x42 + x14*x43 + x17*x44 + x20*x45;
  T x47=(1.0/2.0)*pow(x30, 2);
  T x48=(3.0/2.0)*x32;
  T x49=a[0] - x48 + 1.0/2.0;
  T x50=-2*a[0] + x48;
  T x51=(1.0/2.0)*pow(x13, 2);
  T x52=(3.0/2.0)*x15;
  T x53=a[1] - x52 + 1.0/2.0;
  T x54=-2*a[1] + x52;
  T x55=-x12*x51 + x18*x54 + x19*x9 + x21*x53;
  T x56=x19*x25 - x26*x51 + x27*x54 + x28*x53;
  T x57=x19*x35 - x36*x51 + x37*x54 + x38*x53;
  T x58=x19*x42 - x43*x51 + x44*x54 + x45*x53;
  T x59=(1.0/2.0)*pow(x2, 2);
  T x60=(3.0/2.0)*x4;
  T x61=a[2] - x60 + 1.0/2.0;
  T x62=2*a[2];
  T x63=x60 - x62;
  T x64=-n[15]*x59 + n[31]*x63 + n[47]*x61 + n[63]*x7;
  T x65=n[19]*x63 + n[35]*x61 - n[3]*x59 + n[51]*x7;
  T x66=n[23]*x63 + n[39]*x61 + n[55]*x7 - n[7]*x59;
  T x67=-n[11]*x59 + n[27]*x63 + n[43]*x61 + n[59]*x7;
  T x68=x11*x64 + x14*x65 + x17*x66 + x20*x67;
  T x69=-n[12]*x59 + n[28]*x63 + n[44]*x61 + n[60]*x7;
  T x70=-n[0]*x59 + n[16]*x63 + n[32]*x61 + n[48]*x7;
  T x71=n[20]*x63 + n[36]*x61 - n[4]*x59 + n[52]*x7;
  T x72=n[24]*x63 + n[40]*x61 + n[56]*x7 - n[8]*x59;
  T x73=x11*x69 + x14*x70 + x17*x71 + x20*x72;
  T x74=-n[13]*x59 + n[29]*x63 + n[45]*x61 + n[61]*x7;
  T x75=n[17]*x63 - n[1]*x59 + n[33]*x61 + n[49]*x7;
  T x76=n[21]*x63 + n[37]*x61 + n[53]*x7 - n[5]*x59;
  T x77=n[25]*x63 + n[41]*x61 + n[57]*x7 - n[9]*x59;
  T x78=x11*x74 + x14*x75 + x17*x76 + x20*x77;
  T x79=-n[14]*x59 + n[30]*x63 + n[46]*x61 + n[62]*x7;
  T x80=n[18]*x63 - n[2]*x59 + n[34]*x61 + n[50]*x7;
  T x81=n[22]*x63 + n[38]*x61 + n[54]*x7 - n[6]*x59;
  T x82=-n[10]*x59 + n[26]*x63 + n[42]*x61 + n[58]*x7;
  T x83=x11*x79 + x14*x80 + x17*x81 + x20*x82;
  T x84=3*a[0];
  T x85=x40*x55 - x47*x56 + x49*x58 + x50*x57;
  T x86=x40*x68 - x47*x73 + x49*x83 + x50*x78;
  T x87=3*a[1];
  T x88=1 - x87;
  T x89=x87 - 2;
  T x90=x24*(x19*x64 - x51*x65 + x53*x67 + x54*x66) + x31*(x19*x69 - x51*x70 + x53*x72 + x54*x71) + x34*(x19*x74 - x51*x75 + x53*x77 + x54*x76) + x41*(x19*x79 - x51*x80 + x53*x82 + x54*x81);
  T x91=3*a[2];
  T x92=x91 - 2;
  T x93=1 - x91;
  T x94=(1.0/2.0)*x62 - 1;
  T E=x22*x24 + x29*x31 + x34*x39 + x41*x46;
  if(grad) {
    (*grad)[0]=x22*x40 - x29*x47 + x39*x50 + x46*x49;
    (*grad)[1]=x24*x55 + x31*x56 + x34*x57 + x41*x58;
    (*grad)[2]=x24*x68 + x31*x73 + x34*x78 + x41*x83;
  }
  if(hess) {
    (*hess)(0,0)=a[0]*x22 + x29*x30 + x39*(x84 - 2) + x46*(1 - x84);
    (*hess)(0,1)=x85;
    (*hess)(0,2)=x86;
    (*hess)(1,0)=x85;
    (*hess)(1,1)=x24*(a[1]*x9 + x12*x13 + x18*x89 + x21*x88) + x31*(a[1]*x25 + x13*x26 + x27*x89 + x28*x88) + x34*(a[1]*x35 + x13*x36 + x37*x89 + x38*x88) + x41*(a[1]*x42 + x13*x43 + x44*x89 + x45*x88);
    (*hess)(1,2)=x90;
    (*hess)(2,0)=x86;
    (*hess)(2,1)=x90;
    (*hess)(2,2)=x24*(x11*(a[2]*n[63] - n[15]*x94 + n[31]*x92 + n[47]*x93) + x14*(a[2]*n[51] + n[19]*x92 + n[35]*x93 - n[3]*x94) + x17*(a[2]*n[55] + n[23]*x92 + n[39]*x93 - n[7]*x94) + x20*(a[2]*n[59] - n[11]*x94 + n[27]*x92 + n[43]*x93)) + x31*(x11*(a[2]*n[60] - n[12]*x94 + n[28]*x92 + n[44]*x93) + x14*(a[2]*n[48] - n[0]*x94 + n[16]*x92 + n[32]*x93) + x17*(a[2]*n[52] + n[20]*x92 + n[36]*x93 - n[4]*x94) + x20*(a[2]*n[56] + n[24]*x92 + n[40]*x93 - n[8]*x94)) + x34*(x11*(a[2]*n[61] - n[13]*x94 + n[29]*x92 + n[45]*x93) + x14*(a[2]*n[49] + n[17]*x92 - n[1]*x94 + n[33]*x93) + x17*(a[2]*n[53] + n[21]*x92 + n[37]*x93 - n[5]*x94) + x20*(a[2]*n[57] + n[25]*x92 + n[41]*x93 - n[9]*x94)) + x41*(x11*(a[2]*n[62] - n[14]*x94 + n[30]*x92 + n[46]*x93) + x14*(a[2]*n[50] + n[18]*x92 - n[2]*x94 + n[34]*x93) + x17*(a[2]*n[54] + n[22]*x92 + n[38]*x93 - n[6]*x94) + x20*(a[2]*n[58] - n[10]*x94 + n[26]*x92 + n[42]*x93));
  }
  return E;
}
template <typename T>
void DeformedEnvironment<T>::backward(const Vec3T& a,NodeFetcher n,T dEde,const Vec3T& dEdJ) {
  T x0=1 - a[2];
  T x1=pow(x0, 3);
  T x2=(1.0/72.0)*x1;
  T x3=1 - a[0];
  T x4=pow(x3, 2);
  T x5=1 - a[1];
  T x6=pow(x5, 3);
  T x7=dEdJ[0]*x6;
  T x8=x4*x7;
  T x9=pow(x3, 3);
  T x10=dEdJ[1]*pow(x5, 2);
  T x11=x10*x2;
  T x12=x6*x9;
  T x13=dEdJ[2]*pow(x0, 2);
  T x14=(1.0/72.0)*x13;
  T x15=dEde*x1;
  T x16=(1.0/216.0)*x15;
  T x17=(1.0/36.0)*x6;
  T x18=pow(a[0], 2);
  T x19=(3.0/2.0)*x18;
  T x20=dEdJ[0]*(-2*a[0] + x19);
  T x21=x1*x20;
  T x22=pow(a[0], 3);
  T x23=(1.0/2.0)*x22;
  T x24=-x18 + x23 + 2.0/3.0;
  T x25=(1.0/12.0)*x24;
  T x26=x1*x10;
  T x27=x13*x6;
  T x28=x15*x24;
  T x29=dEdJ[0]*(a[0] - x19 + 1.0/2.0);
  T x30=x1*x17;
  T x31=(1.0/2.0)*x18;
  T x32=(1.0/2.0)*a[0] - x23 + x31 + 1.0/6.0;
  T x33=(1.0/12.0)*x32;
  T x34=dEde*x32;
  T x35=x18*x2;
  T x36=x22*x6;
  T x37=(1.0/36.0)*x9;
  T x38=pow(a[1], 2);
  T x39=(3.0/2.0)*x38;
  T x40=dEdJ[1]*(-2*a[1] + x39);
  T x41=x1*x40;
  T x42=dEdJ[0]*x4;
  T x43=pow(a[1], 3);
  T x44=(1.0/2.0)*x43;
  T x45=-x38 + x44 + 2.0/3.0;
  T x46=(1.0/12.0)*x45;
  T x47=x1*x46;
  T x48=x13*x46;
  T x49=x15*x45;
  T x50=(1.0/6.0)*x45;
  T x51=(1.0/6.0)*x41;
  T x52=x24*x45;
  T x53=(1.0/2.0)*x13;
  T x54=(1.0/6.0)*x15;
  T x55=x1*x50;
  T x56=x32*x53;
  T x57=(1.0/36.0)*x22;
  T x58=dEdJ[0]*x18;
  T x59=dEdJ[1]*(a[1] - x39 + 1.0/2.0);
  T x60=x1*x37;
  T x61=(1.0/2.0)*x38;
  T x62=(1.0/2.0)*a[1] - x44 + x61 + 1.0/6.0;
  T x63=(1.0/12.0)*x62;
  T x64=x1*x63;
  T x65=x13*x63;
  T x66=dEde*x62;
  T x67=(1.0/6.0)*x1;
  T x68=x59*x67;
  T x69=(1.0/6.0)*x62;
  T x70=x24*x62;
  T x71=x1*x29;
  T x72=x32*x66;
  T x73=x1*x57;
  T x74=dEdJ[1]*x38;
  T x75=x2*x74;
  T x76=x42*x43;
  T x77=x43*x9;
  T x78=(1.0/36.0)*x43;
  T x79=x1*x74;
  T x80=x13*x43;
  T x81=x34*x78;
  T x82=x22*x43;
  T x83=pow(a[2], 2);
  T x84=(3.0/2.0)*x83;
  T x85=dEdJ[2]*(-2*a[2] + x84);
  T x86=(1.0/36.0)*x12;
  T x87=pow(a[2], 3);
  T x88=(1.0/2.0)*x87;
  T x89=-x83 + x88 + 2.0/3.0;
  T x90=(1.0/12.0)*x89;
  T x91=x10*x90;
  T x92=dEde*x89;
  T x93=(1.0/6.0)*x6;
  T x94=x89*x93;
  T x95=x24*x85;
  T x96=x24*x89;
  T x97=(1.0/2.0)*x10;
  T x98=dEde*x93;
  T x99=x32*x85;
  T x100=x32*x89;
  T x101=x17*x22;
  T x102=x18*x7;
  T x103=(1.0/6.0)*x9;
  T x104=x40*x89;
  T x105=x50*x85;
  T x106=x45*x89;
  T x107=(1.0/2.0)*x42;
  T x108=dEde*x106;
  T x109=(1.0/6.0)*x22;
  T x110=dEdJ[0]*x31;
  T x111=x59*x89;
  T x112=x69*x85;
  T x113=x62*x89;
  T x114=dEde*x113;
  T x115=x37*x43;
  T x116=x74*x90;
  T x117=(1.0/6.0)*x43;
  T x118=x117*x89;
  T x119=dEdJ[1]*x61;
  T x120=dEde*x117;
  T x121=x43*x57;
  T x122=x43*x58;
  T x123=dEdJ[2]*(a[2] - x84 + 1.0/2.0);
  T x124=(1.0/2.0)*x83;
  T x125=(1.0/2.0)*a[2] + x124 - x88 + 1.0/6.0;
  T x126=(1.0/12.0)*x125;
  T x127=x10*x126;
  T x128=dEde*x125;
  T x129=x24*x93;
  T x130=x125*x93;
  T x131=x125*x97;
  T x132=x123*x32;
  T x133=x128*x32;
  T x134=x123*x50;
  T x135=x125*x40;
  T x136=x125*x45;
  T x137=dEde*x136;
  T x138=x125*x59;
  T x139=x123*x69;
  T x140=x125*x62;
  T x141=x125*x66;
  T x142=x126*x74;
  T x143=x117*x24;
  T x144=x117*x125;
  T x145=x119*x125;
  T x146=dEdJ[2]*x83;
  T x147=(1.0/72.0)*x146;
  T x148=(1.0/72.0)*x87;
  T x149=x10*x148;
  T x150=dEde*x87;
  T x151=(1.0/216.0)*x150;
  T x152=x17*x87;
  T x153=x146*x6;
  T x154=x10*x87;
  T x155=x150*x24;
  T x156=x40*x87;
  T x157=x146*x46;
  T x158=x46*x87;
  T x159=x150*x45;
  T x160=x50*x87;
  T x161=(1.0/6.0)*x156;
  T x162=dEdJ[2]*x124;
  T x163=(1.0/6.0)*x150;
  T x164=x162*x32;
  T x165=x37*x87;
  T x166=x146*x63;
  T x167=x63*x87;
  T x168=(1.0/6.0)*x87;
  T x169=x168*x59;
  T x170=x69*x87;
  T x171=x57*x87;
  T x172=x148*x74;
  T x173=x78*x87;
  T x174=x74*x87;
  T x175=x146*x43;
  n[0]+=-x11*x9 - x12*x14 + x12*x16 - x2*x8;
  n[1]+=x17*x21 + x17*x28 - x25*x26 - x25*x27;
  n[2]+=-x26*x33 - x27*x33 + x29*x30 + x30*x34;
  n[3]+=-x11*x22 - x14*x36 + x16*x36 + x35*x7;
  n[4]+=x37*x41 + x37*x49 - x42*x47 - x48*x9;
  n[5]+=x21*x50 + x24*x51 - x52*x53 + x52*x54;
  n[6]+=x29*x55 + x32*x51 + x34*x55 - x45*x56;
  n[7]+=-x22*x48 + x41*x57 + x47*x58 + x49*x57;
  n[8]+=-x42*x64 + x59*x60 + x60*x66 - x65*x9;
  n[9]+=x21*x69 + x24*x68 - x53*x70 + x54*x70;
  n[10]+=x32*x68 - x56*x62 + x67*x72 + x69*x71;
  n[11]+=-x22*x65 + x58*x64 + x59*x73 + x66*x73;
  n[12]+=-x14*x77 + x16*x77 - x2*x76 + x75*x9;
  n[13]+=x21*x78 + x25*x79 - x25*x80 + x28*x78;
  n[14]+=x1*x81 + x33*x79 - x33*x80 + x71*x78;
  n[15]+=dEdJ[0]*x35*x43 - x14*x82 + x16*x82 + x22*x75;
  n[16]+=-x8*x90 + x85*x86 + x86*x92 - x9*x91;
  n[17]+=x20*x94 + x93*x95 - x96*x97 + x96*x98;
  n[18]+=-x100*x97 + x100*x98 + x29*x94 + x93*x99;
  n[19]+=x101*x85 + x101*x92 + x102*x90 - x22*x91;
  n[20]+=x103*x104 + x103*x108 + x105*x9 - x106*x107;
  n[21]+=x106*x20 + x108*x24 + x40*x96 + x52*x85;
  n[22]+=x100*x40 + x106*x29 + x108*x32 + x45*x99;
  n[23]+=x104*x109 + x105*x22 + x106*x110 + x108*x109;
  n[24]+=x103*x111 + x103*x114 - x107*x113 + x112*x9;
  n[25]+=x113*x20 + x59*x96 + x66*x96 + x70*x85;
  n[26]+=x100*x59 + x100*x66 + x113*x29 + x62*x99;
  n[27]+=x109*x111 + x109*x114 + x110*x113 + x112*x22;
  n[28]+=x115*x85 + x115*x92 + x116*x9 - x76*x90;
  n[29]+=x117*x95 + x118*x20 + x119*x96 + x120*x96;
  n[30]+=x100*x119 + x100*x120 + x117*x99 + x118*x29;
  n[31]+=x116*x22 + x121*x85 + x121*x92 + x122*x90;
  n[32]+=x123*x86 - x126*x8 - x127*x9 + x128*x86;
  n[33]+=x123*x129 + x128*x129 + x130*x20 - x131*x24;
  n[34]+=x130*x29 - x131*x32 + x132*x93 + x133*x93;
  n[35]+=x101*x123 + x101*x128 + x102*x126 - x127*x22;
  n[36]+=x103*x135 + x103*x137 - x107*x136 + x134*x9;
  n[37]+=x123*x52 + x128*x52 + x135*x24 + x136*x20;
  n[38]+=x132*x45 + x135*x32 + x136*x29 + x136*x34;
  n[39]+=x109*x135 + x109*x137 + x110*x136 + x134*x22;
  n[40]+=x103*x138 + x103*x141 - x107*x140 + x139*x9;
  n[41]+=x123*x70 + x128*x70 + x138*x24 + x140*x20;
  n[42]+=x125*x72 + x132*x62 + x138*x32 + x140*x29;
  n[43]+=x109*x138 + x109*x141 + x110*x140 + x139*x22;
  n[44]+=x115*x123 + x115*x128 - x126*x76 + x142*x9;
  n[45]+=x123*x143 + x128*x143 + x144*x20 + x145*x24;
  n[46]+=x117*x132 + x117*x133 + x144*x29 + x145*x32;
  n[47]+=x121*x123 + x121*x128 + x122*x126 + x142*x22;
  n[48]+=x12*x147 + x12*x151 - x148*x8 - x149*x9;
  n[49]+=x152*x20 + x153*x25 - x154*x25 + x155*x17;
  n[50]+=x152*x29 + x152*x34 + x153*x33 - x154*x33;
  n[51]+=x102*x148 + x147*x36 - x149*x22 + x151*x36;
  n[52]+=x156*x37 + x157*x9 - x158*x42 + x159*x37;
  n[53]+=x160*x20 + x161*x24 + x162*x52 + x163*x52;
  n[54]+=x160*x29 + x160*x34 + x161*x32 + x164*x45;
  n[55]+=x156*x57 + x157*x22 + x158*x58 + x159*x57;
  n[56]+=x165*x59 + x165*x66 + x166*x9 - x167*x42;
  n[57]+=x162*x70 + x163*x70 + x169*x24 + x170*x20;
  n[58]+=x164*x62 + x168*x72 + x169*x32 + x170*x29;
  n[59]+=x166*x22 + x167*x58 + x171*x59 + x171*x66;
  n[60]+=x147*x77 - x148*x76 + x151*x77 + x172*x9;
  n[61]+=x155*x78 + x173*x20 + x174*x25 + x175*x25;
  n[62]+=x173*x29 + x174*x33 + x175*x33 + x81*x87;
  n[63]+=x122*x148 + x147*x82 + x151*x82 + x172*x22;
}
template <typename T>
typename DeformedEnvironment<T>::MatT DeformedEnvironment<T>::DSDJ(const Eigen::JacobiSVD<Eigen::Matrix<double,3,3>>& SVD) {
  MatT ret=MatT::Zero(3,9);
  for(int i=0; i<3; i++)
    for(int k=0; k<3; k++)
      for(int j=0; j<3; j++)
        ret(i,j+k*3)=SVD.matrixU()(j,i)*SVD.matrixV()(k,i);
  return ret;
}
template class DeformedEnvironment<FLOAT>;
#ifdef FORCE_ADD_DOUBLE_PRECISION
template struct DeformedEnvironment<double>;
#endif
}
