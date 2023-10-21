//meta-template function to count max number of variables
template <typename T,char LABEL,char LABEL2>
struct NrVar {
  static int nrVar(const SOSPolynomial<T,LABEL>& s) {
    if(LABEL==LABEL2) {
      int nrVar=0;
      for(int i=0; i<(int)s._terms.size(); i++)
        nrVar=std::max<int>(nrVar,s._terms[i].nrVar());
      return nrVar;
    } else return 0;
  }
  static int order(const SOSPolynomial<T,LABEL>& s) {
    if(LABEL==LABEL2)
      return orderAll(s);
    else return 0;
  }
  static int orderAll(const SOSPolynomial<T,LABEL>& s) {
    int order=0;
    for(int i=0; i<(int)s._terms.size(); i++)
      order=std::max<int>(order,s._terms[i].order());
    return order;
  }
};
template <typename T,char LABEL3,char LABEL,char LABEL2>
struct NrVar<SOSPolynomial<T,LABEL3>,LABEL,LABEL2> {
  static int nrVar(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    if(LABEL==LABEL2) {
      int nrVar=0;
      for(int i=0; i<(int)s._terms.size(); i++)
        nrVar=std::max<int>(nrVar,s._terms[i].nrVar());
      return nrVar;
    } else {
      int nrVar=0;
      for(int i=0; i<(int)s._terms.size(); i++)
        nrVar=std::max<int>(nrVar,NrVar<T,LABEL3,LABEL2>::nrVar(s._terms[i]._coef));
      return nrVar;
    }
  }
  static int order(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    if(LABEL==LABEL2) {
      int order=0;
      for(int i=0; i<(int)s._terms.size(); i++)
        order=std::max<int>(order,s._terms[i].order());
      return order;
    } else {
      int order=0;
      for(int i=0; i<(int)s._terms.size(); i++)
        order=std::max<int>(order,NrVar<T,LABEL3,LABEL2>::order(s._terms[i]._coef));
      return order;
    }
  }
  static int orderAll(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    int order=0;
    for(int i=0; i<(int)s._terms.size(); i++)
      order=std::max<int>(order,s._terms[i].order()+NrVar<T,LABEL3,LABEL2>::orderAll(s._terms[i]._coef));
    return order;
  }
};
