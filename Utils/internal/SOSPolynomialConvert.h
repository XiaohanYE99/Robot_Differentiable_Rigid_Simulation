//fromString
template <typename T>
void convertFormatted(const std::string& str,T& val) {
  val=std::stod(str);
}
#ifdef QUADMATH_SUPPORT
template <>
void convertFormatted(const std::string& str,float128& val) {
  val=float128(str);
}
#endif
template <>
void convertFormatted(const std::string& str,mpfr_float& val) {
  val=mpfr_float(str);
}
template <typename T,char LABEL>
void convertFormatted(const std::string& str,SOSTerm<T,LABEL>& val) {
  val.readFormattedString(str);
}
template <typename T,char LABEL>
void convertFormatted(const std::string& str,SOSPolynomial<T,LABEL>& val) {
  val.readFormattedString(str);
}
//fromString
template <typename T>
void convert(const std::string& str,T& val) {
  val=std::stod(str);
}
#ifdef QUADMATH_SUPPORT
template <>
void convert(const std::string& str,float128& val) {
  val=float128(str);
}
#endif
template <>
void convert(const std::string& str,mpfr_float& val) {
  val=mpfr_float(str);
}
template <typename T,char LABEL>
void convert(const std::string& str,SOSTerm<T,LABEL>& val) {
  val << str;
}
template <typename T,char LABEL>
void convert(const std::string& str,SOSPolynomial<T,LABEL>& val) {
  val << str;
}

//toString
template <typename T>
std::string convertFormatted(T val) {
  return std::to_string(val);
}
#ifdef QUADMATH_SUPPORT
template <>
std::string convertFormatted(float128 val) {
  return val.str();
}
#endif
template <>
std::string convertFormatted(mpfr_float val) {
  return val.str();
}
template <typename T,char LABEL>
std::string convertFormatted(const SOSTerm<T,LABEL>& val) {
  return val.formattedString();
}
template <typename T,char LABEL>
std::string convertFormatted(const SOSPolynomial<T,LABEL>& val) {
  return val.formattedString();
}
//toString
template <typename T>
std::string convert(T val,bool& bracketNeeded) {
  bracketNeeded=false;
  return std::to_string(val);
}
#ifdef QUADMATH_SUPPORT
template <>
std::string convert(float128 val,bool& bracketNeeded) {
  bracketNeeded=false;
  return val.str();
}
#endif
template <>
std::string convert(mpfr_float val,bool& bracketNeeded) {
  bracketNeeded=false;
  return val.str();
}
template <typename T,char LABEL>
std::string convert(const SOSTerm<T,LABEL>& val,bool& bracketNeeded) {
  bracketNeeded=false;
  return val.toString();
}
template <typename T,char LABEL>
std::string convert(const SOSPolynomial<T,LABEL>& val,bool& bracketNeeded) {
  bracketNeeded=(int)val._terms.size()>1;
  return val.toString();
}
