//meta-template function to create linear system of equations for MLS solve
template <typename T>
struct Solve {
  typedef ParallelVector<Eigen::Triplet<typename ScalarOfT<T>::Type,int>> STrips;
  static void solve(int& nrEq,STrips& Lss,STrips& Rss,const T& L,const T& R) {}
};
template <typename T,char LABEL>
struct Solve<SOSPolynomial<T,LABEL>> {
  typedef SOSPolynomial<T,LABEL> Poly;
  typedef Eigen::Triplet<typename ScalarOfT<T>::Type,int> STrip;
  typedef ParallelVector<Eigen::Triplet<typename ScalarOfT<T>::Type,int>> STrips;
  static void addL(int nrEq,STrips& Lss,STrips& Rss,const SOSTerm<T,LABEL>& t) {
    if(t._id.empty())
      Rss.push_back(STrip(nrEq,0,-t._coef));
    else {
      ASSERT(t._id.size()==1 && t._order[0]==1)
      Lss.push_back(STrip(nrEq,t._id[0],t._coef));
    }
  }
  static void addR(int nrEq,STrips& Lss,STrips& Rss,const SOSTerm<T,LABEL>& t) {
    if(t._id.empty())
      Rss.push_back(STrip(nrEq,0,t._coef));
    else {
      ASSERT(t._id.size()==1 && t._order[0]==1)
      Lss.push_back(STrip(nrEq,t._id[0],-t._coef));
    }
  }
  static void solve(int& nrEq,STrips& Lss,STrips& Rss,const Poly& L,const Poly& R) {
    int i=0,j=0;
    while(i<(int)L._terms.size() && j<(int)R._terms.size())
      if(L._terms[i]<R._terms[j]) {
        addL(nrEq,Lss,Rss,L._terms[i]);
        i++;
      } else if(L._terms[i]>R._terms[j]) {
        addR(nrEq,Lss,Rss,R._terms[j]);
        j++;
      } else {
        addL(nrEq,Lss,Rss,L._terms[i]);
        addR(nrEq,Lss,Rss,R._terms[j]);
        i++;
        j++;
      }
    while(i<(int)L._terms.size()) {
      addL(nrEq,Lss,Rss,L._terms[i]);
      i++;
    }
    while(j<(int)R._terms.size()) {
      addR(nrEq,Lss,Rss,R._terms[j]);
      j++;
    }
    nrEq++;
  }
};
template <typename T,char LABEL2,char LABEL>
struct Solve<SOSPolynomial<SOSPolynomial<T,LABEL2>,LABEL>> {
  typedef SOSPolynomial<T,LABEL2> PolyT;
  typedef SOSPolynomial<PolyT,LABEL> Poly;
  typedef ParallelVector<Eigen::Triplet<typename ScalarOfT<T>::Type,int>> STrips;
  static void solve(int& nrEq,STrips& Lss,STrips& Rss,const Poly& L,const Poly& R) {
    int i=0,j=0;
    while(i<(int)L._terms.size() && j<(int)R._terms.size())
      if(L._terms[i]<R._terms[j]) {
        Solve<PolyT>::solve(nrEq,Lss,Rss,L._terms[i]._coef,ScalarOfT<PolyT>::convert(0));
        i++;
      } else if(L._terms[i]>R._terms[j]) {
        Solve<PolyT>::solve(nrEq,Lss,Rss,ScalarOfT<PolyT>::convert(0),R._terms[j]._coef);
        j++;
      } else {
        Solve<PolyT>::solve(nrEq,Lss,Rss,L._terms[i]._coef,R._terms[j]._coef);
        i++;
        j++;
      }
    while(i<(int)L._terms.size()) {
      Solve<PolyT>::solve(nrEq,Lss,Rss,L._terms[i]._coef,ScalarOfT<PolyT>::convert(0));
      i++;
    }
    while(j<(int)R._terms.size()) {
      Solve<PolyT>::solve(nrEq,Lss,Rss,ScalarOfT<PolyT>::convert(0),R._terms[j]._coef);
      j++;
    }
  }
};
