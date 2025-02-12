#ifndef EPSILON_H
#define EPSILON_H
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include <float.h>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
#ifdef QUADMATH_SUPPORT
#include <boost/multiprecision/float128.hpp>
#endif

namespace PHYSICSMOTION {
typedef boost::multiprecision::mpq_rational rational;
#ifdef QUADMATH_SUPPORT
typedef boost::multiprecision::float128 float128;
#endif
typedef boost::multiprecision::static_mpfr_float_100 mpfr_float;
template <typename T>
struct Epsilon {
  static T defaultEps() {
    return _defaultEps;
  }
  static T rotationEps() {
    return _rotationEps;
  }
  static T finiteDifferenceEps() {
    return _finiteDifferenceEps;
  }
  static void setDefaultEps(T defaultEps) {
    _defaultEps=defaultEps;
  }
  static void setRotationEps(T rotationEps) {
    _rotationEps=rotationEps;
  }
  static void setFiniteDifferenceEps(T finiteDifferenceEps) {
    _finiteDifferenceEps=finiteDifferenceEps;
  }
 private:
  static T _defaultEps;
  static T _rotationEps;
  static T _finiteDifferenceEps;
};
template <>
struct Epsilon<mpfr_float> {
  static mpfr_float defaultEps();
  static mpfr_float rotationEps();
  static mpfr_float finiteDifferenceEps();
  static void setDefaultEps(mpfr_float defaultEps);
  static void setRotationEps(mpfr_float rotationEps);
  static void setFiniteDifferenceEps(mpfr_float finiteDifferenceEps);
 private:
  static void initialize();
  static mpfr_float _defaultEps;
  static mpfr_float _rotationEps;
  static mpfr_float _finiteDifferenceEps;
  static bool _initialized;
};
}

#endif
