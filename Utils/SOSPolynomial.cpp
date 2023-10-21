#include "SOSPolynomial.h"
#include "DebugGradient.h"
#include "SparseUtils.h"
#include "Utils.h"
#include "IO.h"
#include <stack>
#include <deque>
#include <string>
#include <functional>

namespace PHYSICSMOTION {
#define NO_CONSISTENCY_CHECK
#include "internal/SOSPolynomialLexer.h"
#include "internal/SOSPolynomialSolve.h"
#include "internal/SOSPolynomialConvert.h"
#include "internal/SOSPolynomialEvaluate.h"
#include "internal/SOSPolynomialRobustInversion.h"
#include "internal/SOSPolynomialIsZero.h"
//SOSTerm
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm():_coef(0) {}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(T val):_coef(val) {}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(T val,int id,int order):_coef(val) {
  if(order>0) {
    _id.assign(1,id);
    _order.assign(1,order);
  }
  consistencyCheck();
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(const std::string& str,bool) {
  readFormattedString(str);
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::read(std::istream& is) {
  readBinaryData(_id,is);
  readBinaryData(_order,is);
  readBinaryData(_coef,is);
  return is.good();
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::write(std::ostream& os) const {
  writeBinaryData(_id,os);
  writeBinaryData(_order,os);
  writeBinaryData(_coef,os);
  return os.good();
}
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::type() const {
  return typeid(SOSTerm).name();
}
//cmp
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator<(const SOSTerm& other) const {
  int i=0;
  while(i<(int)_id.size()&&i<(int)other._id.size()) {
    if(_id[i]<other._id[i])
      return true;
    else if(_id[i]>other._id[i])
      return false;
    else {
      if(_order[i]<other._order[i])
        return true;
      else if(_order[i]>other._order[i])
        return false;
    }
    i++;
  }
  if(_id.size()<other._id.size())
    return true;
  else if(_id.size()>other._id.size())
    return false;
  else return false;
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator>(const SOSTerm& other) const {
  return other<*this;
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator==(const SOSTerm& other) const {
  return !(*this<other) && !(other<*this);
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator!=(const SOSTerm& other) const {
  return !(*this==other);
}
//op
template <typename T,char LABEL>
int SOSTerm<T,LABEL>::nrVar() const {
  if(_id.empty())
    return 0;
  else return _id.back()+1;
}
template <typename T,char LABEL>
int SOSTerm<T,LABEL>::order() const {
  if(_order.empty())
    return 0;
  else {
    int order=0;
    for(int i=0; i<(int)_order.size(); i++)
      order+=_order[i];
    return order;
  }
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::hasId(int id) const {
  std::vector<int>::const_iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  return it!=_id.end() && *it==id;
}
template <typename T,char LABEL>
int SOSTerm<T,LABEL>::order(int id) const {
  std::vector<int>::const_iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  auto off=it-_id.begin();
  return _order[off];
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::removeId(int id) const {
  SOSTerm ret=*this;
  std::vector<int>::const_iterator it=std::lower_bound(ret._id.begin(),ret._id.end(),id);
  if(it!=ret._id.end() && *it==id) {
    auto off=it-ret._id.begin();
    ret._id.erase(ret._id.begin()+off);
    ret._order.erase(ret._order.begin()+off);
  }
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator*(T other) const {
  SOSTerm<T,LABEL> ret=*this;
  ret._coef*=other;
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator*=(T other) {
  return *this=*this*other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator*(const SOSTerm& other) const {
  SOSTerm<T,LABEL> ret;
  int i=0,j=0;
  while(i<(int)_id.size() && j<(int)other._id.size()) {
    if(_id[i]<other._id[j]) {
      ret.add(_id[i],_order[i]);
      i++;
    } else if(_id[i]>other._id[j]) {
      ret.add(other._id[j],other._order[j]);
      j++;
    } else {
      ret.add(_id[i],_order[i]+other._order[j]);
      i++;
      j++;
    }
  }
  while(i<(int)_id.size()) {
    ret.add(_id[i],_order[i]);
    i++;
  }
  while(j<(int)other._id.size()) {
    ret.add(other._id[j],other._order[j]);
    j++;
  }
  ret._coef=_coef*other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator*=(const SOSTerm& other) {
  return *this=*this*other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator+(const SOSTerm& other) const {
  ASSERT(*this==other)
  SOSTerm<T,LABEL> ret=*this;
  ret._coef+=other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator+=(const SOSTerm& other) {
  return *this=*this+other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator-(const SOSTerm& other) const {
  ASSERT(*this==other)
  SOSTerm<T,LABEL> ret=*this;
  ret._coef-=other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator-=(const SOSTerm& other) {
  return *this=*this-other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSTerm<T,LABEL>::integrate(int nrVar) const {
  int j=0;
  SOSPolynomial<T,LABEL> ret,p;
  for(int i=0; i<nrVar; i++) {
    int order=0;
    if(j<(int)_id.size()&&_id[j]==i)
      order=_order[j++];
    p =SOSTerm(ScalarOfT<T>::convert(1/double(order+1)),i+nrVar,order+1);
    p-=SOSTerm(ScalarOfT<T>::convert(1/double(order+1)),i,order+1);
    if(i==0)
      ret=p;
    else ret*=p;
  }
  ASSERT(j==(int)_id.size())
  return ret*SOSTerm(_coef);
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::gradient(std::vector<SOSPolynomial<T,LABEL>>& grad) const {
  for(int i=0; i<(int)_id.size(); i++) {
    SOSTerm<T,LABEL> t=*this;
    if(_id.empty())
      continue;
    else if(_order[i]>1) {
      t._order[i]--;
      t._coef*=ScalarOfT<T>::convert(_order[i]);
    } else {
      t._id.erase(t._id.begin()+i);
      t._order.erase(t._order.begin()+i);
    }
    grad[_id[i]]+=t;
  }
  for(int i=0; i<(int)grad.size(); i++)
    grad[i].consistencyCheck();
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::gradient(int row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>>& grad) const {
  for(int i=0; i<(int)_id.size(); i++) {
    SOSTerm<T,LABEL> t=*this;
    if(_id.empty())
      continue;
    else if(_order[i]>1) {
      t._order[i]--;
      t._coef*=ScalarOfT<T>::convert(_order[i]);
    } else {
      t._id.erase(t._id.begin()+i);
      t._order.erase(t._order.begin()+i);
    }
    grad.push_back(Eigen::Triplet<SOSPolynomial<T,LABEL>,int>(row,_id[i],t));
    //grad[_id[i]]+=t;
  }
  //for(int i=0; i<(int)grad.size(); i++)
  //  grad[i].consistencyCheck();
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::integrate(const COLD& L,const COLD& U) const {
  T ret=ScalarOfT<T>::convert(1);
  int j=0;
  for(int i=0; i<L.size(); i++) {
    int order=0;
    if(j<(int)_id.size()&&_id[j]==i)
      order=_order[j++];
    ret*=ScalarOfT<T>::convert((pow(U[i],order+1)-pow(L[i],order+1))/(order+1));
  }
  ASSERT(j==(int)_id.size())
  return _coef*ret;
}
template <typename T,char LABEL>
template <typename HESS>
T SOSTerm<T,LABEL>::eval(const COLD& x,VEC* grad,HESS* hess) const {
  T ret=ScalarOfT<T>::convert(1),fgradX,fhessX,tmp;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(int i=0; i<(int)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  if(grad)
    for(int i=(int)_id.size()-1; i>=0; i--) {
      fgradX=ScalarOfT<T>::convert(_order[i]*pow(x[_id[i]],_order[i]-1));
      ParallelEvaluate<T>::gradient(*grad,_id[i],_coef*fgradX*xx1[i]*xx2[i]);
      if(i>0)
        xx2[i-1]=xx2[i]*xx0[i];
    }
  //hessian
  if(hess)
    for(int i=0; i<(int)_id.size(); i++) {
      if(_order[i]>1) {
        fhessX=ScalarOfT<T>::convert(_order[i]*(_order[i]-1)*pow(x[_id[i]],_order[i]-2));
        ParallelEvaluate<T>::hessian(*hess,_id[i],_id[i],_coef*fhessX*xx1[i]*xx2[i]);
      }
      fgradX=ScalarOfT<T>::convert(_order[i]*pow(x[_id[i]],_order[i]-1));
      tmp=xx1[i]*fgradX;
      for(int j=i+1; j<(int)_id.size(); j++) {
        fgradX=ScalarOfT<T>::convert(_order[j]*pow(x[_id[j]],_order[j]-1));
        ParallelEvaluate<T>::hessianSym(*hess,_id[i],_id[j],_coef*tmp*fgradX*xx2[j]);
        tmp*=xx0[j];
      }
    }
  return _coef*ret;
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalGradDir(const COLD& x,const COLD& dir) const {
  T ret=ScalarOfT<T>::convert(0),fgradX;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(int i=0; i<(int)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  for(int i=(int)_id.size()-1; i>=0; i--) {
    fgradX=ScalarOfT<T>::convert(_order[i]*pow(x[_id[i]],_order[i]-1));
    ret+=ScalarOfT<T>::convert(dir[_id[i]])*_coef*fgradX*xx1[i]*xx2[i];
    if(i>0)
      xx2[i-1]=xx2[i]*xx0[i];
  }
  return ret;
}
template <typename T,char LABEL>
template <typename JAC>
T SOSTerm<T,LABEL>::evalJacTpl(int row,const COLD& x,JAC* jac) const {
  T ret=ScalarOfT<T>::convert(1),fgradX;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(int i=0; i<(int)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  if(jac)
    for(int i=(int)_id.size()-1; i>=0; i--) {
      fgradX=ScalarOfT<T>::convert(_order[i]*pow(x[_id[i]],_order[i]-1));
      ParallelEvaluate<T>::hessian(*jac,row,_id[i],_coef*fgradX*xx1[i]*xx2[i]);
      if(i>0)
        xx2[i-1]=xx2[i]*xx0[i];
    }
  return _coef*ret;
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalJac(int row,const COLD& x,STRIPS* jac) const {
  return evalJacTpl(row,x,jac);
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalJac(int row,const COLD& x,MAT* jac) const {
  return evalJacTpl(row,x,jac);
}
//io
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::toString(const std::unordered_map<int,std::string>* varNames) const {
  bool bracketNeeded=false;
  std::string ret=convert(_coef,bracketNeeded);
  if(bracketNeeded)
    ret="("+ret+")";
  for(int i=0; i<(int)_id.size(); i++) {
    if(varNames)
      ret+="*"+varNames->find(_id[i])->second;
    else ret+="*"+std::string(1,LABEL)+convert(_id[i],bracketNeeded);
    if(_order[i]>1)
      ret+="^"+convert(_order[i],bracketNeeded);
  }
  return ret;
}
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::formattedString(const std::unordered_map<int,std::string>* varNames) const {
  std::string ret="("+convertFormatted(_coef)+")";
  for(int i=0; i<(int)_id.size(); i++) {
    if(varNames)
      ret+="*"+varNames->find(_id[i])->second;
    else ret+="*"+std::string(1,LABEL)+convertFormatted(_id[i]);
    ret+="^"+convertFormatted(_order[i]);
  }
  return ret;
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::readFormattedString(std::string str) {
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  _id.clear();
  _order.clear();

  //coef part
  int j=(int)str.size();
  while(j>0 && str[j-1]!=')')
    j--;
  if(j>0)
    convertFormatted(str.substr(1,j-2),_coef);
  else {
    size_t pos=str.find_first_of('*');
    if(pos==std::string::npos)
      j=(int)str.size();
    else j=(int)pos;
    convertFormatted(str.substr(0,j),_coef);
  }
  if(j==(int)str.size())
    return;

  //other part
  std::istringstream iss(str.substr(j));
  char m,c,o;
  int id,order;
  while(!iss.eof()) {
    iss >> m >> c >> id >> o >> order;
    ASSERT(m=='*' && c==LABEL && o=='^')
    add(id,order);
  }
}
//helper
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::add(int id,int order) {
  std::vector<int>::iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  auto off=it-_id.begin();
  if(it==_id.end() || *it!=id) {
    _id.insert(off+_id.begin(),id);
    _order.insert(off+_order.begin(),order);
  } else {
    _order[off]+=order;
  }
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::consistencyCheck() const {
#ifndef NO_CONSISTENCY_CHECK
  for(int i=0; i<(int)_id.size()-1; i++) {
    ASSERT(_id[i]<_id[i+1])
  }
  for(int i=0; i<(int)_id.size(); i++) {
    ASSERT(_order[i]>0)
  }
#endif
}

//SOSPolynomial
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial() {}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(typename ScalarOfT<T>::Type other) {
  *this=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(const SOSTerm<T,LABEL>& other):_terms(1,other) {}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(const std::string& str,bool) {
  *this<<str;
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(std::istream& is) {
  readBinaryData(_terms,is);
  return is.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(std::ostream& os) const {
  writeBinaryData(_terms,os);
  return os.good();
}
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::type() const {
  return typeid(SOSPolynomial).name();
}
//cmp
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::operator==(const SOSPolynomial& other) const {
  return false;
}
//op
template <typename T,char LABEL>
int SOSPolynomial<T,LABEL>::nrVar() const {
  return nrVar<LABEL>();
}
template <typename T,char LABEL>
int SOSPolynomial<T,LABEL>::order() const {
  return order<LABEL>();
}
template <typename T,char LABEL>
int SOSPolynomial<T,LABEL>::orderAll() const {
  return orderAll<LABEL>();
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(int other) {
  *this=(typename ScalarOfT<T>::Type)other;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(typename ScalarOfT<T>::Type other) {
  *this=ScalarOfT<SOSPolynomial<T,LABEL>>::convert(other);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(const SOSTerm<T,LABEL>& other) {
  _terms.assign(1,other);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(const std::string& str) {
  *this<<str;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const T& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const T& other) {
  for(int i=0; i<(int)_terms.size(); i++)
    _terms[i]._coef*=other;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const SOSTerm<T,LABEL>& other) {
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_terms.size(); i++)
    _terms[i]*=other;
  std::sort(_terms.begin(),_terms.end());
  consistencyCheck();
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const SOSPolynomial& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const SOSPolynomial& other) {
  std::vector<SOSPolynomial> otherByTerms(other._terms.size());
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)other._terms.size(); i++)
    otherByTerms[i]=*this*other._terms[i];
  sum(otherByTerms);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator+(const T& other) const {
  SOSPolynomial ret=*this;
  return ret+=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator+=(const T& other) {
  return operator+=(SOSTerm<T,LABEL>(other));
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator+(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret+=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator+=(const SOSTerm<T,LABEL>& other) {
  typename std::vector<SOSTerm<T,LABEL>>::iterator it=std::lower_bound(_terms.begin(),_terms.end(),other);
  if(it==_terms.end() || *it!=other)
    _terms.insert(it,other);
  else it->_coef+=other._coef;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator+(const SOSPolynomial& other) const {
  //degenerate case
  if(_terms.empty())
    return other;
  else if(other._terms.empty())
    return *this;
  else if(other._terms.size()==1)
    return operator+(other._terms[0]);
  else if(_terms.size()==1)
    return other.operator+(_terms[0]);
  //merge sort
  int i=0,j=0;
  SOSPolynomial ret;
  while(i<(int)_terms.size() && j<(int)other._terms.size()) {
    if(_terms[i]<other._terms[j]) {
      ret.add(_terms[i]);
      i++;
    } else if(_terms[i]>other._terms[j]) {
      ret.add(other._terms[j]);
      j++;
    } else {
      ret.add(_terms[i]+other._terms[j]);
      i++;
      j++;
    }
  }
  while(i<(int)_terms.size()) {
    ret.add(_terms[i]);
    i++;
  }
  while(j<(int)other._terms.size()) {
    ret.add(other._terms[j]);
    j++;
  }
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator+=(const SOSPolynomial& other) {
  return *this=*this+other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator-(const T& other) const {
  SOSPolynomial ret=*this;
  return ret-=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator-=(const T& other) {
  return operator-=(SOSTerm<T,LABEL>(other));
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator-(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret-=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator-=(const SOSTerm<T,LABEL>& other) {
  SOSTerm<T,LABEL> otherNeg=other;
  otherNeg._coef*=ScalarOfT<T>::convert(-1);
  return operator+=(otherNeg);
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator-(const SOSPolynomial& other) const {
  SOSPolynomial ret=*this;
  return ret-=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator-=(const SOSPolynomial& other) {
  T coef=ScalarOfT<T>::convert(-1);
  SOSPolynomial<T,LABEL> ret=other*coef;
  return operator+=(ret);
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::integrate() const {
  if(_terms.empty())
    return SOSPolynomial();
  SOSPolynomial ret=_terms[0].integrate(nrVar());
  for(int i=1; i<(int)_terms.size(); i++)
    ret+=_terms[i].integrate(nrVar());
  return ret;
}
template <typename T,char LABEL>
std::vector<SOSPolynomial<T,LABEL>> SOSPolynomial<T,LABEL>::gradient(int orderLimit) const {
  std::vector<SOSPolynomial<T,LABEL>> ret(nrVar());
  for(int i=0; i<(int)_terms.size(); i++)
    if(_terms[i].order()<=orderLimit)
      _terms[i].gradient(ret);
  return ret;
}
template <typename T,char LABEL>
std::vector<std::vector<SOSPolynomial<T,LABEL>>> SOSPolynomial<T,LABEL>::hessian(int orderLimit) const {
  std::vector<std::vector<SOSPolynomial>> ret;
  std::vector<SOSPolynomial<T,LABEL>> g=gradient(orderLimit);
  for(int i=0; i<(int)g.size(); i++)
    ret.push_back(g[i].gradient());
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::VECP SOSPolynomial<T,LABEL>::gradientV(int orderLimit) const {
  VECP ret;
  ret.resize(nrVar());//=VECP::Zero(nrVar());
  std::vector<SOSPolynomial<T,LABEL>> retV=gradient(orderLimit);
  for(int i=0; i<ret.size(); i++)
    ret[i]=retV[i];
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MATP SOSPolynomial<T,LABEL>::hessianM(int orderLimit) const {
  MATP ret;
  ret.resize(nrVar(),nrVar());//=MATP::Zero(nrVar(),nrVar());
  std::vector<std::vector<SOSPolynomial<T,LABEL>>> hessV=hessian(orderLimit);
  for(int i=0; i<(int)hessV.size(); i++)
    for(int j=0; j<(int)hessV[i].size(); j++)
      ret(i,j)=hessV[i][j];
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::gradientSparse(int row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>>& grad,int orderLimit) const {
  for(int i=0; i<(int)_terms.size(); i++)
    if(_terms[i].order()<=orderLimit)
      _terms[i].gradient(row,grad);
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::SMATP SOSPolynomial<T,LABEL>::gradientSparse(int orderLimit) const {
  ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>> trips;
  gradientSparse(0,trips,orderLimit);

  SMATP ret;
  ret.resize(1,nrVar());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret.transpose();
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::SMATP SOSPolynomial<T,LABEL>::hessianSparse(int orderLimit) const {
  std::vector<SOSPolynomial<T,LABEL>> g=gradient(orderLimit);
  ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,int>> trips;
  for(int i=0; i<(int)g.size(); i++)
    g[i].gradientSparse(i,trips);

  SMATP ret;
  ret.resize(nrVar(),nrVar());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::sum(const std::vector<SOSPolynomial>& polys) {
  //merge
  _terms.clear();
  for(int i=0; i<(int)polys.size(); i++)
    _terms.insert(_terms.end(),polys[i]._terms.begin(),polys[i]._terms.end());
  if(_terms.empty())
    return;
  makeConsistent();
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::VEC SOSPolynomial<T,LABEL>::gradientCoef() const {
  VEC gc=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
  std::vector<SOSPolynomial<T,LABEL>> g=gradient();
  for(int i=0; i<nrVar(); i++)
    gc[i]=g.at(i).operator T();
  return gc;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MAT SOSPolynomial<T,LABEL>::hessianCoef() const {
  MAT hc=MAT::Constant(nrVar(),nrVar(),ScalarOfT<T>::convert(0));
  std::vector<std::vector<SOSPolynomial>> h=hessian();
  for(int i=0; i<nrVar(); i++)
    for(int j=0; j<(int)h.at(i).size(); j++)
      hc(i,j)=h.at(i).at(j).operator T();
  return hc;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MAT SOSPolynomial<T,LABEL>::JTJCoef() const {
  int nrT=(int)_terms.size();
  MAT JTJ=MAT::Constant(nrT,nrT,ScalarOfT<T>::convert(0));
  OMP_PARALLEL_FOR_
  for(int r=0; r<nrT; r++)
    for(int c=0; c<nrT; c++)
      JTJ(r,c)=_terms[r]._coef*_terms[c]._coef;
  return JTJ;
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::integrate(const COLD& L,const COLD& U) const {
  if(_terms.empty())
    return T();
  std::vector<T> termsI(_terms.size()),tmp;
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_terms.size(); i++)
    termsI[i]=_terms[i].integrate(L,U);
  while(termsI.size()>1) {
    tmp.resize((termsI.size()+1)/2);
    OMP_PARALLEL_FOR_
    for(int k=0; k<(int)tmp.size(); k++)
      if(k*2==(int)termsI.size()-1)
        tmp[k]=termsI[k*2];
      else tmp[k]=termsI[k*2]+termsI[k*2+1];
    termsI.swap(tmp);
  }
  return termsI[0];
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::eval(const COLD& x,VEC* grad,MAT* hess) const {
  return ParallelEvaluate<T>::eval(*this,x,grad,hess);
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::evalTrips(const COLD& x,VEC* grad,STRIPS* hess) const {
  return ParallelEvaluate<T>::eval(*this,x,grad,hess);
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::evalGradDir(const COLD& x,const COLD& dir) const {
  T ret=ScalarOfT<T>::convert(0);
  for(int i=0; i<(int)_terms.size(); i++)
    ret+=_terms[i].evalGradDir(x,dir);
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::retainOrder(int orderLimit) const {
  SOSPolynomial ret;
  for(int i=0; i<(int)_terms.size(); i++)
    if(_terms[i].order()<=orderLimit)
      ret._terms.push_back(_terms[i]);
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::pruneOrder(int orderLimit) const {
  SOSPolynomial ret;
  for(int i=0; i<(int)_terms.size(); i++)
    if(_terms[i].order()>orderLimit)
      ret._terms.push_back(_terms[i]);
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::rename(const std::vector<int>& ids) const {
  //collect ids
  std::unordered_map<int,int> idMap;
  for(int i=0; i<(int)ids.size(); i++)
    idMap[i]=ids[i];
  //remap
  SOSPolynomial ret=*this;
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)ret._terms.size(); i++) {
    const SOSTerm<T,LABEL>& t0=_terms[i];
    SOSTerm<T,LABEL>& t=ret._terms[i];
    t._id.clear();
    t._order.clear();
    for(int d=0; d<(int)t0._id.size(); d++)
      t.add(idMap.find(t0._id[d])->second,t0._order[d]);
  }
  std::sort(ret._terms.begin(),ret._terms.end());
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::varRemap(const std::vector<int>& ids) const {
  //remap
  ParallelVector<SOSTerm<T,LABEL>> terms;
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_terms.size(); i++) {
    const SOSTerm<T,LABEL>& term=_terms[i];
    SOSTerm<T,LABEL> termRet(term._coef);
    for(int d=0; d<(int)term._id.size(); d++) {
      int newId=ids.at(term._id[d]);
      ASSERT(newId>=0)
      termRet.add(newId,term._order[d]);
    }
    terms.push_back(termRet);
  }
  SOSPolynomial ret;
  ret._terms.insert(ret._terms.end(),terms.begin(),terms.end());
  ret.makeConsistent();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::setAllCoef(T coef) const {
  SOSPolynomial ret=*this;
  for(int i=0; i<(int)ret._terms.size(); i++)
    ret._terms[i]._coef=coef;
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::removeZero(typename ScalarOfT<T>::Type eps) const {
  return IsZero<SOSPolynomial<T,LABEL>>::removeZero(*this,eps);
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::removeVariableId(int id) const {
  SOSPolynomial ret=*this;
  for(int i=0; i<(int)ret._terms.size(); i++) {
    SOSTerm<T,LABEL>& t=ret._terms[i];
    for(int j=0; j<(int)t._id.size(); j++) {
      ASSERT(t._id[j]!=id)
      if(t._id[j]>id)
        t._id[j]--;
    }
  }
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearConstraint(int id,const SOSPolynomial& cons) const {
  SOSPolynomial ret;
  for(int i=0; i<(int)_terms.size(); i++) {
    SOSTerm<T,LABEL> t=_terms[i];
    if(t.hasId(id)) {
      int order=t.order(id);
      SOSTerm<T,LABEL> t2=t.removeId(id);
      SOSPolynomial p=cons;
      for(int i=1; i<order; i++)
        p*=cons;
      ret+=p*t2;
    } else {
      ret+=t;
    }
  }
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearTransform(const std::unordered_map<int,typename ScalarOfT<T>::Type>& cons) const {
  ParallelVector<SOSTerm<T,LABEL>> terms;
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_terms.size(); i++) {
    const SOSTerm<T,LABEL>& term=_terms[i];
    SOSTerm<T,LABEL> termRet(term._coef);
    for(int v=0; v<(int)term._id.size(); v++) {
      typename std::unordered_map<int,typename ScalarOfT<T>::Type>::const_iterator it=cons.find(term._id[v]);
      if(it!=cons.end()) {
        termRet._coef*=ScalarOfT<T>::convert(pow(it->second,term._order[v]));
      } else {
        termRet._id.push_back(term._id[v]);
        termRet._order.push_back(term._order[v]);
      }
    }
    terms.push_back(termRet);
  }
  SOSPolynomial ret;
  ret._terms.insert(ret._terms.end(),terms.begin(),terms.end());
  ret.makeConsistent();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearTransform(const std::unordered_map<int,SOSPolynomial>& cons,int reportInterval) const {
  SOSPolynomial ret;
  std::vector<std::vector<SOSPolynomial>> cache(nrVar(),std::vector<SOSPolynomial>(1,ScalarOfT<SOSPolynomial>::convert(1)));
  for(int i=0; i<(int)_terms.size(); i++) {
    if(reportInterval>0 && (i%reportInterval)==0) {
      std::cout << "linearTransform: " << i << "/" << (int)_terms.size() << std::endl;
    }
    const SOSTerm<T,LABEL>& t=_terms[i];
    SOSPolynomial sRet=SOSTerm<T,LABEL>(t._coef);
    for(int d=0; d<(int)t._id.size(); d++) {
      std::vector<SOSPolynomial>& cacheId=cache[t._id[d]];
      if((int)cacheId.size()<=t._order[d]) {
        SOSPolynomial tmpP=cacheId.back();
        const SOSPolynomial& p=cons.find(t._id[d])->second;
        for(int o=(int)cacheId.size()-1; o<(int)t._order[d]; o++) {
          tmpP*=p;
          cacheId.push_back(tmpP);
        }
      }
      sRet*=cacheId[t._order[d]];
    }
    ret+=sRet;
  }
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::COLD SOSPolynomial<T,LABEL>::solve(const std::vector<SOSPolynomial>& LHS,const std::vector<SOSPolynomial>& RHS) {
  int nrEQ=0;
  ParallelVector<Eigen::Triplet<typename ScalarOfT<T>::Type,int>> Lss,Rss;
  for(int i=0; i<(int)LHS.size(); i++)
    Solve<SOSPolynomial>::solve(nrEQ,Lss,Rss,LHS[i],RHS[i]);

  Eigen::SparseMatrix<typename ScalarOfT<T>::Type,0,int> LHSM,RHSM;
  LHSM.resize(nrEQ,nrEQ);
  LHSM.setFromTriplets(Lss.begin(),Lss.end());
  RHSM.resize(nrEQ,1);
  RHSM.setFromTriplets(Rss.begin(),Rss.end());

  COLD RHSD=RHSM.toDense();
  MATD LHSD=LHSM.toDense();
  //std::cout << LHSD << std::endl << RHSD << std::endl;
  //analyze nullspace rows
  int nrRow=0;
  for(int i=0; i<LHSD.rows(); i++)
    if(LHSD.row(i).cwiseAbs().maxCoeff()!=0.0) {
      LHSD.row(nrRow)=LHSD.row(i);
      RHSD[nrRow++]=RHSD[i];
    } else {
      ASSERT(RHSD[i]==0.0)
    }
  LHSD=LHSD.block(0,0,nrRow,nrRow).eval();
  RHSD=RHSD.segment(0,nrRow).eval();
  COLD ret=RobustInversion<T,LABEL>::eval(LHSD)*RHSD;
  //std::cout << (LHSD*ret-RHSD) << std::endl;
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::operator T() const {
  if(_terms.empty())
    return ScalarOfT<T>::convert(0);
  else {
    ASSERT_MSGV(_terms.size()==1 && _terms[0]._id.empty(),"Polynomial(%s) cannot be converted to a constant!",toString().c_str())
    return _terms[0]._coef;
  }
}
//string
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::toString(const std::unordered_map<int,std::string>* varNames) const {
  if(_terms.empty())
    return std::string();
  std::string ret=_terms[0].toString(varNames);
  for(int i=1; i<(int)_terms.size(); i++)
    ret+="+"+_terms[i].toString(varNames);
  return ret;
}
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::formattedString(const std::unordered_map<int,std::string>* varNames) const {
  if(_terms.empty())
    return std::string();
  std::string ret=_terms[0].formattedString(varNames);
  for(int i=1; i<(int)_terms.size(); i++)
    ret+="+"+_terms[i].formattedString(varNames);
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::readFormattedString(std::string str) {
  _terms.clear();
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  int depth=0;
  std::vector<int> cut;
  cut.push_back(-1);
  for(int i=0; i<(int)str.size(); i++) {
    if(str[i]=='(')
      depth++;
    else if(str[i]==')')
      depth--;
    else if(str[i]=='+'&&depth==0)
      cut.push_back(i);
  }
  cut.push_back((int)str.size());
  for(int i=1; i<(int)cut.size(); i++) {
    int p0=cut[i-1]+1,p1=cut[i];
    SOSTerm<T,LABEL> t;
    t.readFormattedString(str.substr(p0,p1-p0));
    _terms.push_back(t);
  }
  std::sort(_terms.begin(),_terms.end());
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::operator<<(std::string str) {
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  SOSPolynomialLexer<T,LABEL> parser;
  //parser.testLex(str);
  *this=parser.parse(str);
}
//debug
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugIntegrate() const {
  bool bracketNeeded=false;
  COLD L=COLD::Random(nrVar()),U=COLD::Random(nrVar()),LU=concat(L,U);
  T polyIEval=integrate().eval(LU),polyIEvalRef=polyIEval-integrate(L,U);
  std::cout << "Integrate: " << convert(polyIEval,bracketNeeded) << " err: " << convert(polyIEvalRef,bracketNeeded) << std::endl;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugGradient() const {
  bool bracketNeeded=false;
  COLD x=COLD::Random(nrVar());
  VEC g=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
  eval(x,&g);
  std::vector<SOSPolynomial> grad=gradient();
  for(int i=0; i<(int)grad.size(); i++) {
    T gradV=grad[i].eval(x),gradVRef=gradV-g[i];
    std::cout << "GradientA: " << convert(gradV,bracketNeeded) << " err: " << convert(gradVRef,bracketNeeded) << std::endl;
  }
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugEval() const {
  bool bracketNeeded=false;
  DEFINE_NUMERIC_DELTA_T(typename ScalarOfT<T>::Type)
  COLD x=COLD::Random(nrVar());
  COLD dx=COLD::Random(nrVar());
  VEC g=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0)),g2=g;
  MAT h=MAT::Constant(nrVar(),nrVar(),ScalarOfT<T>::convert(0)),h2=h;
  T f=eval(x,&g,&h);
  T f2=eval(x+dx*(typename ScalarOfT<T>::Type)(DELTA),&g2,&h2);
  SMATP hess=hessianSparse();
  {
    T gdxV=ScalarOfT<T>::convert(0),gdxVRef=ScalarOfT<T>::convert(0);
    for(int i=0; i<g.size(); i++) {
      gdxV+=g[i]*ScalarOfT<T>::convert(dx[i]);
      gdxVRef+=g[i]*ScalarOfT<T>::convert(dx[i]);
    }
    gdxVRef-=(f2-f)*ScalarOfT<T>::convert(1/DELTA);
    std::cout << "Gradient: " << convert(gdxV,bracketNeeded) << " err: " << convert(gdxVRef,bracketNeeded) << std::endl;
  }
  {
    VEC hdxV=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
    VEC hdxVRef=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
    for(int i=0; i<g.size(); i++)
      for(int j=0; j<g.size(); j++) {
        hdxV[i]+=h(i,j)*ScalarOfT<T>::convert(dx[j]);
        hdxVRef[i]+=h(i,j)*ScalarOfT<T>::convert(dx[j]);
      }
    hdxVRef-=(g2-g)*ScalarOfT<T>::convert(1/DELTA);
    T hdxVNorm=ScalarOfT<T>::convert(0),hdxVRefNorm=ScalarOfT<T>::convert(0);
    for(int i=0; i<hdxV.size(); i++) {
      hdxVNorm+=hdxV[i]*hdxV[i];
      hdxVRefNorm+=hdxVRef[i]*hdxVRef[i];
    }
    std::cout << "Hessian: " << convert(hdxVNorm,bracketNeeded) << " err: " << convert(hdxVRefNorm,bracketNeeded) << std::endl;
  }
  {
    MAT hdx=hessianGradDir(hess,x,dx).toDense();
    MAT hdxVRef=hdx-(h2-h)*ScalarOfT<T>::convert(1/DELTA);
    T hdxVNorm=ScalarOfT<T>::convert(0),hdxVRefNorm=ScalarOfT<T>::convert(0);
    for(int r=0; r<hdx.rows(); ++r)
      for(int c=0; c<hdx.cols(); ++c)
        hdxVNorm+=hdx(r,c)*hdx(r,c);
    for(int r=0; r<hdxVRef.rows(); ++r)
      for(int c=0; c<hdxVRef.cols(); ++c)
        hdxVRefNorm+=hdxVRef(r,c)*hdxVRef(r,c);
    std::cout << "HessianGradDir: " << convert(hdxVNorm,bracketNeeded) << " err: " << convert(hdxVRefNorm,bracketNeeded) << std::endl;
  }
}
//helper
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::add(const SOSTerm<T,LABEL>& other) {
  typename std::vector<SOSTerm<T,LABEL>>::iterator it=std::lower_bound(_terms.begin(),_terms.end(),other);
  if(it==_terms.end()||*it!=other)
    _terms.insert(it,other);
  else *it+=other;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::consistencyCheck() const {
#ifndef NO_CONSISTENCY_CHECK
  for(int i=0; i<(int)_terms.size()-1; i++) {
    ASSERT(_terms[i]<_terms[i+1])
  }
  for(int i=0; i<(int)_terms.size(); i++)
    _terms[i].consistencyCheck();
#endif
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::makeConsistent() {
  //sort
  std::sort(_terms.begin(),_terms.end());
  //compact
  int j=0;
  for(int i=1; i<(int)_terms.size(); i++)
    if(_terms[i]==_terms[j])
      _terms[j]._coef+=_terms[i]._coef;
    else _terms[++j]=_terms[i];
  _terms.resize(j+1);
  consistencyCheck();
}
template class SOSTerm<FLOAT,'x'>;
template class SOSPolynomial<FLOAT,'x'>;
#ifdef FORCE_ADD_DOUBLE_PRECISION
template class SOSTerm<double,'x'>;
template class SOSPolynomial<double,'x'>;
#endif
}
