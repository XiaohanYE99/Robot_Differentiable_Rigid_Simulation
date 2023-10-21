#ifndef INTERP_H
#define INTERP_H

namespace PHYSICSMOTION {
//linear interpolation
template <typename T,typename TF>
T interp1D(const T& v0,const T& v1,
           const TF& px) {
  return v0*(TF)(1.0f-px)+v1*(TF)px;
}
template <typename T,typename TF>
T interp2D(const T& v0,const T& v1,
           const T& v2,const T& v3,
           const TF& px,const TF& py) {
  return interp1D(interp1D(v0,v1,px),interp1D(v2,v3,px),py);
}
template <typename T,typename TF>
T interp3D(const T& v0,const T& v1,const T& v2,const T& v3,
           const T& v4,const T& v5,const T& v6,const T& v7,
           const TF& px,const TF& py,const TF& pz) {
  return interp1D(interp2D(v0,v1,v2,v3,px,py),interp2D(v4,v5,v6,v7,px,py),pz);
}
//linear interpolation stencil
template <typename TF>
int stencil1D(TF* coefs,const TF& px) {
  coefs[0]=1.0f-px;
  coefs[1]=px;
  return 2;
}
template <typename TF>
int stencil2D(TF* coefs,const TF& px,const TF& py) {
  coefs[0]=(1.0f-px)*(1.0f-py);
  coefs[1]=px*(1.0f-py);
  coefs[2]=(1.0f-px)*py;
  coefs[3]=px*py;
  return 4;
}
template <typename TF>
int stencil3D(TF* coefs,const TF& px,const TF& py,const TF& pz) {
  coefs[0]=(1.0f-px)*(1.0f-py)*(1.0f-pz);
  coefs[1]=px*(1.0f-py)*(1.0f-pz);
  coefs[2]=(1.0f-px)*py*(1.0f-pz);
  coefs[3]=px*py*(1.0f-pz);
  coefs[4]=(1.0f-px)*(1.0f-py)*pz;
  coefs[5]=px*(1.0f-py)*pz;
  coefs[6]=(1.0f-px)*py*pz;
  coefs[7]=px*py*pz;
  return 8;
}
//linear interpolation grad
template <typename T,typename TF>
T interp1DGrad(const T& v0,const T& v1,const TF&) {
  return v1-v0;
}
template <typename T,typename TF>
std::tuple<T,T> interp2DGrad(const T& v0,const T& v1,const T& v2,const T& v3,const TF& px,const TF& py) {
  std::tuple<T,T> ret;
  std::get<0>(ret)=interp1D(interp1DGrad(v0,v1,px),
                            interp1DGrad(v2,v3,px),py);
  std::get<1>(ret)=interp1DGrad(interp1D(v0,v1,px),
                                interp1D(v2,v3,px),py);
  return ret;
}
template <typename T,typename TF>
std::tuple<T,T,T> interp3DGrad(const T& v0,const T& v1,const T& v2,const T& v3,
                               const T& v4,const T& v5,const T& v6,const T& v7,
                               const TF& px,const TF& py,const TF& pz) {
  std::tuple<T,T,T> ret;
  std::tuple<T,T> dxy0=interp2DGrad(v0,v1,v2,v3,px,py);
  std::tuple<T,T> dxy1=interp2DGrad(v4,v5,v6,v7,px,py);
  std::get<0>(ret)=interp1D(std::get<0>(dxy0),std::get<0>(dxy1),pz);
  std::get<1>(ret)=interp1D(std::get<1>(dxy0),std::get<1>(dxy1),pz);
  std::get<2>(ret)=interp1DGrad(interp2D(v0,v1,v2,v3,px,py),
                                interp2D(v4,v5,v6,v7,px,py),pz);
  return ret;
}
//linear interpolation grad stencil
template <typename TF>
int stencil1DGradStencil(TF* coefs,const TF&) {
  coefs[0]=-1;
  coefs[1]=1;
  return 2;
}
template <typename TF>
int stencil2DGradXStencil(TF* coefs,const TF&,const TF& py) {
  coefs[1]=-(1.0f-py);
  coefs[0]=(1.0f-py);
  coefs[2]=-py;
  coefs[3]=py;
  return 4;
}
template <typename TF>
int stencil2DGradYStencil(TF* coefs,const TF& px,const TF&) {
  coefs[0]=-(1.0f-px);
  coefs[1]=-px;
  coefs[2]=(1.0f-px);
  coefs[3]=px;
  return 4;
}
template <typename TF>
int stencil3DGradXStencil(TF* coefs,const TF&,const TF& py,const TF& pz) {
  coefs[0]=-(1.0f-py)*(1.0f-pz);
  coefs[1]=(1.0f-py)*(1.0f-pz);
  coefs[2]=-py*(1.0f-pz);
  coefs[3]=py*(1.0f-pz);
  coefs[4]=-(1.0f-py)*pz;
  coefs[5]=(1.0f-py)*pz;
  coefs[6]=-py*pz;
  coefs[7]=py*pz;
  return 8;
}
template <typename TF>
int stencil3DGradYStencil(TF* coefs,const TF& px,const TF&,const TF& pz) {
  coefs[0]=-(1.0f-px)*(1.0f-pz);
  coefs[1]=-px*(1.0f-pz);
  coefs[2]=(1.0f-px)*(1.0f-pz);
  coefs[3]=px*(1.0f-pz);
  coefs[4]=-(1.0f-px)*pz;
  coefs[5]=-px*pz;
  coefs[6]=(1.0f-px)*pz;
  coefs[7]=px*pz;
  return 8;
}
template <typename TF>
int stencil3DGradZStencil(TF* coefs,const TF& px,const TF& py,const TF&) {
  coefs[0]=-(1.0f-px)*(1.0f-py);
  coefs[1]=-px*(1.0f-py);
  coefs[2]=-(1.0f-px)*py;
  coefs[3]=-px*py;
  coefs[4]=(1.0f-px)*(1.0f-py);
  coefs[5]=px*(1.0f-py);
  coefs[6]=(1.0f-px)*py;
  coefs[7]=px*py;
  return 8;
}
//cubic interpolation
template <typename T>
T interp1DCubic(T t,std::function<T(int)> f) {
  T t2=t*t,t3=t2*t;
  T C0=-(t-2*t2+t3)/2;
  T C1=(2-5*t2+3*t3)/2;
  T C2=(t+4*t2-3*t3)/2;
  T C3=(t3-t2)/2;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
template <typename T>
T interp1DCubicDiff(T t,std::function<T(int)> f) {
  T t2=t*t;
  T C0=-(1-4*t+3*t2)/2;
  T C1=(9*t2-10*t)/2;
  T C2=(1+8*t-9*t2)/2;
  T C3=(3*t2-2*t)/2;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
template <typename T>
T interp1DCubicDDiff(T t,std::function<T(int)> f) {
  T C0=2-3*t;
  T C1=9*t-5;
  T C2=4-9*t;
  T C3=3*t-1;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
//bezier interpoltion
template <typename T>
T interp1DCubicBezier(T t,std::function<T(int)> f) {
  T C0=pow(1-t,3);
  T C1=3*pow(1-t,2)*t;
  T C2=3*(1-t)*pow(t,2);
  T C3=pow(t,3);
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
}

#endif
