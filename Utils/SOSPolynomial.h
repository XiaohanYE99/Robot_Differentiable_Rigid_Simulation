#ifndef SOS_POLYNOMIAL_H
#define SOS_POLYNOMIAL_H

#include "Serializable.h"
#include "Pragma.h"
#include "ParallelVector.h"
#include <unordered_map>

namespace PHYSICSMOTION {
template <typename T,char LABEL>
class SOSTerm;
template <typename T,char LABEL>
class SOSPolynomial;
class SOSInfo {
 public:
  std::vector<int> _id,_order;
};
template <typename T2>
struct Zero<SOSPolynomial<T2,'a'>> {
  static T2 value() {
    return SOSPolynomial<T2,'a'>();
  }
  static T2 value(const T2&) {
    return SOSPolynomial<T2,'a'>();
  }
};
#include "internal/SOSPolynomialScalarOfT.h"
#include "internal/SOSPolynomialAffineTransXId.h"
#include "internal/SOSPolynomialContract.h"
#include "internal/SOSPolynomialRearrange.h"
#include "internal/SOSPolynomialCast.h"
#include "internal/SOSPolynomialNrVar.h"
//SOSTerm
template <typename T,char LABEL>
class SOSTerm : public SOSInfo, public SerializableBase {
 public:
  typedef Eigen::Matrix<T,-1,1> VEC;
  typedef Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,int>> STRIPS;
  typedef Eigen::Triplet<T,int> STRIP;
  SOSTerm();
  SOSTerm(T val);
  SOSTerm(T val,int id,int order);
  SOSTerm(const std::string& str,bool placeholder); //this last placeholder is to resolve some signature ambiguity bug under float128/mpfr
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual std::string type() const;
  //cmp
  bool operator<(const SOSTerm& other) const;
  bool operator>(const SOSTerm& other) const;
  bool operator==(const SOSTerm& other) const;
  bool operator!=(const SOSTerm& other) const;
  //op
  int nrVar() const;
  int order() const;
  bool hasId(int id) const;
  int order(int id) const;
  SOSTerm removeId(int id) const;
  SOSTerm operator*(T other) const;
  SOSTerm& operator*=(T other);
  SOSTerm operator*(const SOSTerm& other) const;
  SOSTerm& operator*=(const SOSTerm& other);
  SOSTerm operator+(const SOSTerm& other) const;
  SOSTerm& operator+=(const SOSTerm& other);
  SOSTerm operator-(const SOSTerm& other) const;
  SOSTerm& operator-=(const SOSTerm& other);
  SOSPolynomial<T,LABEL> integrate(int nrVar) const;
  void gradient(std::vector<SOSPolynomial<T,LABEL>>& grad) const;
  void gradient(int row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>>& grad) const;
  T integrate(const COLD& L,const COLD& U) const;
  template <typename HESS>
  T eval(const COLD& x,VEC* grad=NULL,HESS* hess=NULL) const;
  T evalGradDir(const COLD& x,const COLD& dir) const;
  template <typename JAC>
  T evalJacTpl(int row,const COLD& x,JAC* jac=NULL) const;
  T evalJac(int row,const COLD& x,STRIPS* jac=NULL) const;
  T evalJac(int row,const COLD& x,MAT* jac=NULL) const;
  //misc
  template <typename T2>
  typename Cast<SOSTerm,T2>::Type cast() const {
    typename Cast<SOSTerm,T2>::Type ret;
    ret._coef=Cast<T,T2>::cast(_coef);
    ret._id=_id;
    ret._order=_order;
    return ret;
  }
  //io
  std::string toString(const std::unordered_map<int,std::string>* varNames=NULL) const;
  std::string formattedString(const std::unordered_map<int,std::string>* varNames=NULL) const;
  void readFormattedString(std::string str);
  void add(int id,int order);
  void consistencyCheck() const;
  //data
  T _coef;
};
//SOSPolynomial
template <typename T,char LABEL>
class SOSPolynomial : public SerializableBase {
 public:
  typedef Eigen::Matrix<T,-1,1> VEC;
  typedef Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::SparseMatrix<T,0,int> SMAT;
  typedef Eigen::Matrix<SOSPolynomial<T,LABEL>,-1,1> VECP;
  typedef Eigen::Matrix<SOSPolynomial<T,LABEL>,-1,-1> MATP;
  typedef Eigen::SparseMatrix<SOSPolynomial<T,LABEL>,0,int> SMATP;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,int>> STRIPS;
  typedef Eigen::Triplet<T,int> STRIP;
  SOSPolynomial();
  SOSPolynomial(typename ScalarOfT<T>::Type other);
  SOSPolynomial(const SOSTerm<T,LABEL>& other);
  SOSPolynomial(const std::string& str,bool placeholder); //this last placeholder is to resolve some signature ambiguity bug under float128/mpfr
  template <typename T2,char LABEL2>
  SOSPolynomial(const SOSPolynomial<T2,LABEL2>& other) {
    Rearrange<T,LABEL,T2,LABEL2>::arrange(*this,other,std::unordered_map<char,SOSInfo>());
  }
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual std::string type() const;
  //cmp
  bool operator==(const SOSPolynomial& other) const;
  //op
  int nrVar() const;
  template <char LABEL2>
  int nrVar() const {
    return NrVar<T,LABEL,LABEL2>::nrVar(*this);
  }
  int order() const;
  template <char LABEL2>
  int order() const {
    return NrVar<T,LABEL,LABEL2>::order(*this);
  }
  int orderAll() const;
  template <char LABEL2>
  int orderAll() const {
    return NrVar<T,LABEL,LABEL2>::orderAll(*this);
  }
  SOSPolynomial& operator=(int other);
  SOSPolynomial& operator=(typename ScalarOfT<T>::Type other);
  SOSPolynomial& operator=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial& operator=(const std::string& str);
  SOSPolynomial operator*(const T& other) const;
  SOSPolynomial& operator*=(const T& other);
  SOSPolynomial operator*(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator*=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator*(const SOSPolynomial& other) const;
  SOSPolynomial& operator*=(const SOSPolynomial& other);
  SOSPolynomial operator+(const T& other) const;
  SOSPolynomial& operator+=(const T& other);
  SOSPolynomial operator+(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator+=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator+(const SOSPolynomial& other) const;
  SOSPolynomial& operator+=(const SOSPolynomial& other);
  SOSPolynomial operator-(const T& other) const;
  SOSPolynomial& operator-=(const T& other);
  SOSPolynomial operator-(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator-=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator-(const SOSPolynomial& other) const;
  SOSPolynomial& operator-=(const SOSPolynomial& other);
  SOSPolynomial integrate() const;
  std::vector<SOSPolynomial> gradient(int orderLimit=std::numeric_limits<int>::max()) const;
  std::vector<std::vector<SOSPolynomial>> hessian(int orderLimit=std::numeric_limits<int>::max()) const;
  VECP gradientV(int orderLimit=std::numeric_limits<int>::max()) const;
  MATP hessianM(int orderLimit=std::numeric_limits<int>::max()) const;
  void gradientSparse(int row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>>& grad,int orderLimit=std::numeric_limits<int>::max()) const;
  SMATP gradientSparse(int orderLimit=std::numeric_limits<int>::max()) const;
  SMATP hessianSparse(int orderLimit=std::numeric_limits<int>::max()) const;
  void sum(const std::vector<SOSPolynomial>& polys);
  VEC gradientCoef() const;
  MAT hessianCoef() const;
  MAT JTJCoef() const;
  T integrate(const COLD& L,const COLD& U) const;
  T eval(const COLD& x,VEC* grad=NULL,MAT* hess=NULL) const;
  T evalTrips(const COLD& x,VEC* grad=NULL,STRIPS* hess=NULL) const;
  T evalGradDir(const COLD& x,const COLD& dir) const;
  template <typename JAC>
  T evalJac(int row,const COLD& x,JAC* jac=NULL) const {
    T ret=ScalarOfT<T>::convert(0);
    for(int i=0; i<(int)_terms.size(); i++)
      ret+=_terms[i].evalJac(row,x,jac);
    return ret;
  }
  template <typename JAC>
  VEC evalJac(const COLD& x,JAC* jac=NULL) const {
    int nrT=(int)_terms.size();
    VEC vec=VEC::Constant(nrT,ScalarOfT<T>::convert(0));
    for(int row=0; row<nrT; row++)
      vec[row]=_terms[row].evalJac(row,x,jac);
    return vec;
  }
  SOSPolynomial retainOrder(int orderLimit) const;
  SOSPolynomial pruneOrder(int orderLimit) const;
  SOSPolynomial rename(const std::vector<int>& ids) const;
  SOSPolynomial varRemap(const std::vector<int>& ids) const;
  SOSPolynomial setAllCoef(T coef) const;
  SOSPolynomial removeZero(typename ScalarOfT<T>::Type eps) const;
  SOSPolynomial removeVariableId(int id) const;
  SOSPolynomial linearConstraint(int id,const SOSPolynomial& cons) const;
  SOSPolynomial linearTransform(const std::unordered_map<int,typename ScalarOfT<T>::Type>& cons) const;
  SOSPolynomial linearTransform(const std::unordered_map<int,SOSPolynomial>& cons,int reportInterval=0) const;
  static COLD solve(const std::vector<SOSPolynomial>& LHS,const std::vector<SOSPolynomial>& RHS);
  //misc
  static SMAT hessianGradDir(const SMATP& h,const COLD& x,const COLD& dir) {
    STRIPS trips;
    for(int k=0; k<h.outerSize(); ++k)
      for(typename SMATP::InnerIterator it(h,k); it; ++it)
        trips.push_back(STRIP((int)it.row(),(int)it.col(),it.value().evalGradDir(x,dir)));

    SMAT ret;
    ret.resize(h.rows(),h.cols());
    ret.setFromTriplets(trips.begin(),trips.end());
    return ret;
  }
  template <char LABEL2>
  typename Contract<SOSPolynomial,LABEL2>::Type eval(const COLD& x) const {
    return Contract<SOSPolynomial,LABEL2>::eval(*this,x);
  }
  template <char LABEL2>
  SOSPolynomial affineTransXId(int coef,int off) const {
    return AffineTransXId<SOSPolynomial,LABEL2>::transXId(*this,coef,off);
  }
  template <typename T2>
  typename Cast<SOSPolynomial,T2>::Type cast() const {
    typename Cast<SOSPolynomial,T2>::Type ret;
    ret._terms.resize(_terms.size());
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_terms.size(); i++)
      ret._terms[i]=_terms[i].template cast<T2>();
    return ret;
  }
  operator T() const;
  //string
  std::string toString(const std::unordered_map<int,std::string>* varNames=NULL) const;
  std::string formattedString(const std::unordered_map<int,std::string>* varNames=NULL) const;
  void readFormattedString(std::string str);
  void operator<<(std::string str);
  //debug
  void debugIntegrate() const;
  void debugGradient() const;
  void debugEval() const;
  void add(const SOSTerm<T,LABEL>& other);
  void consistencyCheck() const;
  void makeConsistent();
  //data
  std::vector<SOSTerm<T,LABEL>> _terms;
};
//common type
#define DECL_POLYTYPES(TT,POSTFIX)  \
typedef SOSTerm<TT,'x'> TermX##POSTFIX;  \
typedef SOSPolynomial<TT,'x'> PolyX##POSTFIX;  \
typedef SOSTerm<TT,'t'> TermT##POSTFIX;  \
typedef SOSPolynomial<TT,'t'> PolyT##POSTFIX;  \
typedef SOSTerm<TT,'a'> TermA##POSTFIX;  \
typedef SOSPolynomial<TT,'a'> PolyA##POSTFIX;  \
typedef SOSTerm<PolyX##POSTFIX,'a'> TermXA##POSTFIX;  \
typedef SOSPolynomial<PolyX##POSTFIX,'a'> PolyXA##POSTFIX;  \
typedef SOSTerm<PolyA##POSTFIX,'x'> TermAX##POSTFIX;  \
typedef SOSPolynomial<PolyA##POSTFIX,'x'> PolyAX##POSTFIX;  \
typedef SOSTerm<PolyT##POSTFIX,'x'> TermTX##POSTFIX;  \
typedef SOSPolynomial<PolyT##POSTFIX,'x'> PolyTX##POSTFIX;  \
typedef SOSTerm<PolyTX##POSTFIX,'a'> TermTXA##POSTFIX;  \
typedef SOSPolynomial<PolyTX##POSTFIX,'a'> PolyTXA##POSTFIX;
DECL_POLYTYPES(double,D)
#ifdef QUADMATH_SUPPORT
DECL_POLYTYPES(float128,Q)
#endif
DECL_POLYTYPES(mpfr_float,M)
#undef DECL_POLYTYPES
}

namespace Eigen {
namespace internal {
template <typename T>
struct scalar_product_traits<PHYSICSMOTION::SOSPolynomial<T,'a'>,T> {
  enum {
    Defined = 1
  };
  typedef PHYSICSMOTION::SOSPolynomial<T,'a'> ReturnType;
};
template <typename T>
struct scalar_product_traits<T,PHYSICSMOTION::SOSPolynomial<T,'a'>> {
  enum {
    Defined = 1
  };
  typedef PHYSICSMOTION::SOSPolynomial<T,'a'> ReturnType;
};
}
}

#endif
