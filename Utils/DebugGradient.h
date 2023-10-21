#ifndef DEBUG_GRADIENT_H
#define DEBUG_GRADIENT_H

#include "Epsilon.h"
#include <iostream>

//numeric delta
#define DEFINE_NUMERIC_DELTA_T(T) T DELTA=Epsilon<T>::finiteDifferenceEps();

//gradient debug
#define DEBUG_GRADIENT(NAME,A,B) \
if(fabs(B) > sqrt(DELTA)) { \
  std::cout << "\033[1;31m" << NAME << ": " << A << " Err: " << B << "\033[0m" << std::endl; \
} else { \
  std::cout << NAME << ": " << A << " Err: " << B << std::endl; \
}

#define DEBUG_GRADIENT_REL(NAME,A,B) \
if(fabs(B) > sqrt(DELTA)*fabs(A)) { \
  std::cout << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m" << std::endl; \
} else {  \
  std::cout << NAME << ": " << A << " Err: " << B << std::endl; \
}

#define DEBUG_GRADIENT_ASSERT(NAME,A,B) \
if(fabs(B) > sqrt(DELTA)) { \
  std::cout << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m" << std::endl; \
  ASSERT(false) \
} else { \
  std::cout << NAME << ": " << A << " Err: " << B << std::endl; \
}

#define DEBUG_GRADIENT_REL_ASSERT(NAME,A,B) \
if(fabs(B) > sqrt(DELTA)*fabs(A)) { \
  std::cout << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m" << std::endl; \
  ASSERT(false) \
} else { \
  std::cout << NAME << ": " << A << " Err: " << B << std::endl; \
}

#endif
