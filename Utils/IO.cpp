#include "IO.h"
#include "Pragma.h"

namespace PHYSICSMOTION {
//io for basic types
std::ostream& writeBinaryData(const bool& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(bool& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::ostream& writeBinaryData(const char& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(char& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::ostream& writeBinaryData(const int& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(int& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::ostream& writeBinaryData(const unsigned int& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(unsigned int& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::ostream& writeBinaryData(const float& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(float& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::ostream& writeBinaryData(const double& val,std::ostream& os,IOData* dat) {
  os.write((char*)&val,sizeof(val));
  return os;
}
std::istream& readBinaryData(double& val,std::istream& is,IOData* dat) {
  is.read((char*)&val,sizeof(val));
  return is;
}
std::istream& readBinaryData(std::string& str,std::istream& is,IOData* dat) {
  int len;
  readBinaryData(len,is,dat);
  str.assign(len,' ');
  is.read(&(str[0]),len);
  return is;
}
std::ostream& writeBinaryData(const std::string& str,std::ostream& os,IOData* dat) {
  const int len=(int)str.length();
  writeBinaryData(len,os,dat);
  return os.write(str.c_str(),len);
}
std::ostream& writeBinaryData(const rational& val,std::ostream& os,IOData* dat) {
  //std::cout << val.str() << std::endl;
  std::string s=val.str();
  writeBinaryData(s,os);
  return os;
}
std::istream& readBinaryData(rational& val,std::istream& is,IOData* dat) {
  std::string str;
  readBinaryData(str,is);
  val=rational(str);
  return is;
}
#ifdef QUADMATH_SUPPORT
std::ostream& writeBinaryData(const float128& val,std::ostream& os,IOData* dat) {
  //std::cout << val.str() << std::endl;
  writeBinaryData(val.str(),os);
  return os;
}
std::istream& readBinaryData(float128& val,std::istream& is,IOData* dat) {
  std::string str;
  readBinaryData(str,is);
  val=float128(str);
  return is;
}
#endif
std::ostream& writeBinaryData(const mpfr_float& val,std::ostream& os,IOData* dat) {
  //std::cout << val.str() << std::endl;
  writeBinaryData(val.str(),os);
  return os;
}
std::istream& readBinaryData(mpfr_float& val,std::istream& is,IOData* dat) {
  std::string str;
  readBinaryData(str,is);
  val=mpfr_float(str);
  return is;
}
}
