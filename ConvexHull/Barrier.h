#ifndef BARRIER_H
#define BARRIER_H

#include <random>
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
struct Px {
  //P(d(x)-d0) in the paper
  template <typename T>
  static T power(T base, int pow) {
    T res=1;
    for(int i=0; i<pow; ++i) {
      res*=base;
    }
    return res;
  }
  template <typename T>
  T eval(T d,T* D,T* DD,T d0,T coef) const {
    T P;
    T x=d-d0;
    if(x<0)
      return std::numeric_limits<T>::infinity();
    else if(x<_x0) {
      P=coef*pow(_x0-x,3)/pow(x,4);
      if(D)
        *D=coef*(x-4*_x0)*pow(x-_x0,2)/pow(x,5);
//        *D=coef*(-4*pow(_x0-x,3)/pow(x,5)-3*pow(_x0-x,2)/pow(x,4));
      if(DD) {
        T firstPart=20;
        for(int i=0; i<3; ++i) {
          firstPart*=(_x0-x)/x;
        }
        T secondPart=24;
        for(int i=0; i<2; ++i) {
          secondPart*=(_x0-x)/x;
        }
        T thirdPart=6*(_x0-x)/x;
        T sum=firstPart+secondPart+thirdPart;
        for(int i=0; i<3; ++i) {
          sum/=x;
        }
        *DD=coef*sum;
//        *DD=coef*(20*pow(_x0-x,3)/pow(x,6)+24*pow(_x0-x,2)/pow(x,5)+6*(_x0-x)/pow(x,4));
      }
    } else {
      P=0;
      if(D)
        *D=0;
      if(DD)
        *DD=0;
    }
    return P;
  }
  template <typename T>
  void debug(T perturbRange) const {
    DEFINE_NUMERIC_DELTA_T(T)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.,_x0);
    T d=distribution(gen)+perturbRange;
    T D1=0,DD1=0;
    T D2=0,DD2=0;
    T coef=distribution(gen);
    T E1=eval<T>(d,&D1,&DD1,0,coef);
    T dx=distribution(gen);
    d+=dx*DELTA;
    T E2=eval<T>(d,&D2,&DD2,0,coef);
    DEBUG_GRADIENT("Px dE",D1*dx,abs(D1*dx-(E2-E1)/DELTA))
    DEBUG_GRADIENT("Px dG",DD1*dx,abs(DD1*dx-(D2-D1)/DELTA))
  }
  double _x0;
};
struct Logx {
  template <typename T>
  T eval(T d,T*D,T* DD,T d0,T coef) const {
    T P;
    T x=d-d0;
    if(x<0)
      return std::numeric_limits<T>::infinity();
    else {
      P=-coef*log(x);
      if(D)
        *D=-coef/x;
      if(DD)
        *DD=coef/(x*x);
    }
    return P;
  }
  template <typename T>
  void debug(T perturbRange) const {
    DEFINE_NUMERIC_DELTA_T(T)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.,1);
    T d=distribution(gen)+perturbRange;
    T D1=0,DD1=0;
    T D2=0,DD2=0;
    T E1=eval<T>(d,&D1,&DD1,0,1);
    T dx=distribution(gen);
    d+=dx*DELTA;
    T E2=eval<T>(d,&D2,&DD2,0,1);
    DEBUG_GRADIENT("dE",D1*dx,abs(D1*dx-(E2-E1)/DELTA))
    DEBUG_GRADIENT("dG",DD1*dx,abs(DD1*dx-(D2-D1)/DELTA))
  }
  double _x0;
};
struct CLogx {
  //P(d(x)-d0) in the paper
  template <typename T>
  T eval(T d,T* D,T* DD,T d0,T coef) const {
    T P;
    T x=d-d0;
    if(x<0)
      return std::numeric_limits<T>::infinity();
    else if(x<_x0) {
      T logVal=log(x/_x0);
      P=-logVal*(x-_x0)*(x-_x0)/x;
      if(D)
        *D=-2*logVal*(x-_x0)/x+
           (logVal-1)*(x-_x0)*(x-_x0)/x/x;
      if(DD)
        *DD=4*(logVal-1)*(x-_x0)/x/x-
            2*logVal*((x-_x0)*(x-_x0)/x/x/x+1/x)+
            3*(x-_x0)*(x-_x0)/x/x/x;
    } else {
      P=0;
      if(D)
        *D=0;
      if(DD)
        *DD=0;
    }
    return P;
  }
  template <typename T>
  void debug(T perturbRange) const {
    DEFINE_NUMERIC_DELTA_T(T)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.,1);
    T d=distribution(gen)+perturbRange;
    T D1=0,DD1=0;
    T D2=0,DD2=0;
    T E1=eval<T>(d,&D1,&DD1,0,1);
    T dx=distribution(gen);
    d+=dx*DELTA;
    T E2=eval<T>(d,&D2,&DD2,0,1);
    DEBUG_GRADIENT("dE",D1*dx,abs(D1*dx-(E2-E1)/DELTA))
    DEBUG_GRADIENT("dG",DD1*dx,abs(DD1*dx-(D2-D1)/DELTA))
  }
  double _x0;
};
struct InvQuadraticx {
  template <typename T>
  T eval(T d,T*D,T* DD,T d0,T coef) const {
    T P;
    T x=d-d0;
    if(x<0)
      return std::numeric_limits<T>::infinity();
    else {
      P=coef/(x*x);
      if(D)
        *D=-2*coef/(x*x*x);
      if(DD)
        *DD=6*coef/(x*x*x*x);
    }
    return P;
  }
  template <typename T>
  void debug(T perturbRange) const {
    DEFINE_NUMERIC_DELTA_T(T)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.,1);
    T d=distribution(gen)+perturbRange;
    T D1=0,DD1=0;
    T D2=0,DD2=0;
    T E1=eval<T>(d,&D1,&DD1,0,1);
    T dx=distribution(gen);
    d+=dx*DELTA;
    T E2=eval<T>(d,&D2,&DD2,0,1);
    DEBUG_GRADIENT("dE",D1*dx,abs(D1*dx-(E2-E1)/DELTA))
    DEBUG_GRADIENT("dG",DD1*dx,abs(DD1*dx-(D2-D1)/DELTA))
  }
  double _x0;
};
struct Cubicx {
  //(d(x)-d0)^3 for testing
  template <typename T>
  T eval(T d,T*D,T* DD,T d0,T coef) const {
    T P;
    T x=d-d0;
    P=x*x*x*coef;
    if(D)
      *D=x*x*3*coef;
    if(DD)
      *DD=x*6*coef;
    return P;
  }
  template <typename T>
  void debug(T perturbRange) const {
    DEFINE_NUMERIC_DELTA_T(T)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.,1);
    T d=distribution(gen)+perturbRange;
    T D1=0,DD1=0;
    T D2=0,DD2=0;
    T E1=eval<T>(d,&D1,&DD1,0,1);
    T dx=distribution(gen);
    d+=dx*DELTA;
    T E2=eval<T>(d,&D2,&DD2,0,1);
    DEBUG_GRADIENT("dE",D1*dx,abs(D1*dx-(E2-E1)/DELTA))
    DEBUG_GRADIENT("dG",DD1*dx,abs(DD1*dx-(D2-D1)/DELTA))
  }
  double _x0;
};
}
#endif
