#ifndef IO_H
#define IO_H

#include "Serializable.h"
#include <Eigen/Sparse>
#include <fstream>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace PHYSICSMOTION {
//io for basic types
extern std::ostream& writeBinaryData(const bool& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(bool& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const char& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(char& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const int& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(int& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const unsigned int& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(unsigned int& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const float& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(float& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const double& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(double& val,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const std::string& str,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(std::string& str,std::istream& is,IOData* dat=NULL);
extern std::ostream& writeBinaryData(const rational& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(rational& val,std::istream& is,IOData* dat=NULL);
#ifdef QUADMATH_SUPPORT
extern std::ostream& writeBinaryData(const float128& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(float128& val,std::istream& is,IOData* dat=NULL);
#endif
extern std::ostream& writeBinaryData(const mpfr_float& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(mpfr_float& val,std::istream& is,IOData* dat=NULL);
//io for matrix
template <typename type,int size1,int size2>
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat=NULL) {
  int d0=(int)v.rows();
  int d1=(int)v.cols();
  os.write((char*)&d0,sizeof(int));
  os.write((char*)&d1,sizeof(int));
  for(int r=0; r<d0; r++)
    for(int c=0; c<d1; c++)
      writeBinaryData(v(r,c),os);
  return os;
}
template <typename type,int size1,int size2>
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat=NULL) {
  int d0;
  int d1;
  is.read((char*)&d0,sizeof(int));
  is.read((char*)&d1,sizeof(int));
  v.resize(d0,d1);
  for(int r=0; r<d0; r++)
    for(int c=0; c<d1; c++)
      readBinaryData(v(r,c),is);
  return is;
}
//io for quaternion
template <typename type>
std::ostream& writeBinaryData(const Eigen::Quaternion<type>& v,std::ostream& os,IOData* dat=NULL) {
  writeBinaryData(v.w(),os);
  writeBinaryData(v.x(),os);
  writeBinaryData(v.y(),os);
  writeBinaryData(v.z(),os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Quaternion<type>& v,std::istream& is,IOData* dat=NULL) {
  readBinaryData(v.w(),is);
  readBinaryData(v.x(),is);
  readBinaryData(v.y(),is);
  readBinaryData(v.z(),is);
  return is;
}
template <typename type>
std::ostream& writeBinaryData(const Eigen::Translation<type,3>& v,std::ostream& os,IOData* dat=NULL) {
  writeBinaryData(v.x(),os);
  writeBinaryData(v.y(),os);
  writeBinaryData(v.z(),os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Translation<type,3>& v,std::istream& is,IOData* dat=NULL) {
  readBinaryData(v.x(),is);
  readBinaryData(v.y(),is);
  readBinaryData(v.z(),is);
  return is;
}
template <typename type>
std::ostream& writeBinaryData(const Eigen::Transform<type,3,Eigen::Affine>& v,std::ostream& os,IOData* dat=NULL) {
  Eigen::Matrix<type,3,3> r=v.translation();
  Eigen::Matrix<type,3,1> t=v.linear();
  writeBinaryData(r,os);
  writeBinaryData(t,os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Transform<type,3,Eigen::Affine>& v,std::istream& is,IOData* dat=NULL) {
  Eigen::Matrix<type,3,3> r;
  Eigen::Matrix<type,3,1> t;
  readBinaryData(r,is);
  readBinaryData(t,is);
  v.translation()=r;
  v.linear()=t;
  return is;
}

//io for advanced types
//io for pair
template <typename T1,typename T2>
inline std::ostream& writeBinaryData(const std::pair<T1,T2>& m,std::ostream& os,IOData* dat=NULL) {
  writeBinaryData(m.first,os,dat);
  writeBinaryData(m.second,os,dat);
  return os;
}
template <typename T1,typename T2>
inline std::istream& readBinaryData(std::pair<T1,T2>& m,std::istream& is,IOData* dat=NULL) {
  readBinaryData(m.first,is,dat);
  readBinaryData(m.second,is,dat);
  return is;
}
//io for tuple
template <std::size_t ID,typename... T>
struct IOBinaryDataTupleI {
  static inline std::ostream& write(const std::tuple<T...>& m,std::ostream& os,IOData* dat=NULL) {
    writeBinaryData(std::get<ID>(m),os,dat);
    return IOBinaryDataTupleI<ID-1,T...>::write(m,os,dat);
  }
  static inline std::istream& read(std::tuple<T...>& m,std::istream& is,IOData* dat=NULL) {
    readBinaryData(std::get<ID>(m),is,dat);
    return IOBinaryDataTupleI<ID-1,T...>::read(m,is,dat);
  }
};
template <typename... T>
struct IOBinaryDataTupleI<0,T...> {
  static inline std::ostream& write(const std::tuple<T...>& m,std::ostream& os,IOData* dat=NULL) {
    return writeBinaryData(std::get<0>(m),os,dat);
  }
  static inline std::istream& read(std::tuple<T...>& m,std::istream& is,IOData* dat=NULL) {
    return readBinaryData(std::get<0>(m),is,dat);
  }
};
template <typename... T>
inline std::ostream& writeBinaryData(const std::tuple<T...>& m,std::ostream& os,IOData* dat=NULL) {
  return IOBinaryDataTupleI<std::tuple_size<std::tuple<T...>>::value-1,T...>::write(m,os,dat);
}
template <typename... T>
inline std::istream& readBinaryData(std::tuple<T...>& m,std::istream& is,IOData* dat=NULL) {
  return IOBinaryDataTupleI<std::tuple_size<std::tuple<T...>>::value-1,T...>::read(m,is,dat);
}
//io for vector
template <typename T,typename ALLOC>inline std::istream& readBinaryData(std::vector<T,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  for(int i=0; i<size; i++) {
    v.push_back(T());
    readBinaryData(v.back(),is,dat);
  }
  return is;
}
template <typename T,typename ALLOC>inline std::ostream& writeBinaryData(const std::vector<T,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::vector<T,ALLOC>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++)
    writeBinaryData(*beg,os,dat);
  return os;
}
//io for deque
template <typename T,typename ALLOC>inline std::istream& readBinaryData(std::deque<T,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  for(int i=0; i<size; i++) {
    v.push_back(T());
    readBinaryData(v.back(),is,dat);
  }
  return is;
}
template <typename T,typename ALLOC>inline std::ostream& writeBinaryData(const std::deque<T,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::deque<T>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++)
    writeBinaryData(*beg,os,dat);
  return os;
}
//io for list
template <typename T,typename ALLOC>inline std::istream& readBinaryData(std::list<T,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  for(int i=0; i<size; i++) {
    v.push_back(T());
    readBinaryData(v.back(),is,dat);
  }
  return is;
}
template <typename T,typename ALLOC>inline std::ostream& writeBinaryData(const std::list<T,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::list<T>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++)
    writeBinaryData(*beg,os,dat);
  return os;
}
//io for set
template <typename T,typename CMP,typename ALLOC>inline std::istream& readBinaryData(std::set<T,CMP,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  T tmp;
  for(int i=0; i<size; i++) {
    readBinaryData(tmp,is,dat);
    v.insert(tmp);
  }
  return is;
}
template <typename T,typename CMP,typename ALLOC>inline std::ostream& writeBinaryData(const std::set<T,CMP,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::set<T,CMP,ALLOC>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++)
    writeBinaryData(*beg,os,dat);
  return os;
}
//io for map
template <typename K,typename T,typename CMP,typename ALLOC>inline std::istream& readBinaryData(std::map<K,T,CMP,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  K tmpK;
  T tmpV;
  for(int i=0; i<size; i++) {
    readBinaryData(tmpK,is,dat);
    readBinaryData(tmpV,is,dat);
    v[tmpK]=tmpV;
  }
  return is;
}
template <typename K,typename T,typename CMP,typename ALLOC>inline std::ostream& writeBinaryData(const std::map<K,T,CMP,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::map<K,T,CMP,ALLOC>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++) {
    writeBinaryData(beg->first,os,dat);
    writeBinaryData(beg->second,os,dat);
  }
  return os;
}
//io for hash set
template <typename T,typename H,typename P,typename ALLOC>inline std::istream& readBinaryData(std::unordered_set<T,H,P,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  T tmp;
  for(int i=0; i<size; i++) {
    readBinaryData(tmp,is,dat);
    v.insert(tmp);
  }
  return is;
}
template <typename T,typename H,typename P,typename ALLOC>inline std::ostream& writeBinaryData(const std::unordered_set<T,H,P,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::unordered_set<T,H,P,ALLOC>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++)
    writeBinaryData(*beg,os,dat);
  return os;
}
//io for hash map
template <typename K,typename T,typename H,typename P,typename ALLOC>inline std::istream& readBinaryData(std::unordered_map<K,T,H,P,ALLOC>& v,std::istream& is,IOData* dat=NULL,bool extend=false) {
  int size;
  readBinaryData(size,is,dat);
  if(!extend)
    v.clear();
  K tmpK;
  T tmpV;
  for(int i=0; i<size; i++) {
    readBinaryData(tmpK,is,dat);
    readBinaryData(tmpV,is,dat);
    v[tmpK]=tmpV;
  }
  return is;
}
template <typename K,typename T,typename H,typename P,typename ALLOC>inline std::ostream& writeBinaryData(const std::unordered_map<K,T,H,P,ALLOC>& v,std::ostream& os,IOData* dat=NULL) {
  int size=(int)v.size();
  writeBinaryData(size,os,dat);
  for(typename std::unordered_map<K,T,H,P,ALLOC>::const_iterator beg=v.begin(),end=v.end(); beg!=end; beg++) {
    writeBinaryData(beg->first,os,dat);
    writeBinaryData(beg->second,os,dat);
  }
  return os;
}
//sparseIO
template <typename T,int O,typename I>
bool readBinaryData(Eigen::SparseMatrix<T,O,I>& m,std::istream& is,IOData* dat=NULL) {
  int rows,cols;
  std::vector<I> r;
  std::vector<I> c;
  std::vector<T> v;
  readBinaryData(rows,is);
  readBinaryData(cols,is);
  readBinaryData(r,is);
  readBinaryData(c,is);
  readBinaryData(v,is);
  m.resize(rows,cols);
  m.reserve(v.size());
  for(int ci=0; ci<m.cols(); ++ci) {
    m.startVec(ci);
    for(int off=c[ci]; off<c[ci+1]; off++)
      m.insertBack(r[off],ci)=v[off];
  }
  m.finalize();
  return is.good();
}
template <typename T,int O,typename I>
bool writeBinaryData(const Eigen::SparseMatrix<T,O,I>& m,std::ostream& os,IOData* dat=NULL) {
  std::vector<I> r(m.nonZeros());
  std::vector<I> c(m.cols()+1);
  std::vector<T> v(m.nonZeros());
  for(int k=0,offr=0; k<m.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(m,k); it; ++it) {
      v[offr]=it.value();
      r[offr++]=it.row();
      c[k+1]++;
    }
  for(int k=0; k<m.outerSize(); ++k)
    c[k+1]+=c[k];
  writeBinaryData((int)m.rows(),os);
  writeBinaryData((int)m.cols(),os);
  writeBinaryData(r,os);
  writeBinaryData(c,os);
  writeBinaryData(v,os);
  return os.good();
}
}

#endif
