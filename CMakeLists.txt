CMAKE_MINIMUM_REQUIRED(VERSION 3.21)
PROJECT(PhysicsMotion)

#COMPILE SETTINGS
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Compiler optimization level")
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
  MESSAGE(STATUS "Build type: Debug")
  IF(MSVC)
    SET(CMAKE_CXX_FLAGS "/openmp /bigobj /Wall /std:c++17 ${CMAKE_CXX_FLAGS}")
  ELSE(MSVC)
    ADD_COMPILE_DEFINITIONS(_GLIBCXX_DEBUG=1)
    ADD_COMPILE_DEFINITIONS(_GLIBCXX_USE_CXX11_ABI=1)
    SET(CMAKE_CXX_FLAGS "-fopenmp -O0 -Wall -fPIC -std=c++17 -Wno-unused-function -Wno-unused-local-typedefs ${CMAKE_CXX_FLAGS}")
  ENDIF(MSVC)
ELSE()
  MESSAGE(STATUS "Build type: Release")
  IF(MSVC)
    SET(CMAKE_CXX_FLAGS "/openmp /bigobj /Wall /std:c++17 ${CMAKE_CXX_FLAGS}")
  ELSE(MSVC)
    ADD_COMPILE_DEFINITIONS(_GLIBCXX_USE_CXX11_ABI=1)
    SET(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -fPIC -std=c++17 -Wno-unused-function -Wno-unused-local-typedefs ${CMAKE_CXX_FLAGS}")
  ENDIF(MSVC)
ENDIF()
SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})

#C++11 ABI
OPTION(USE_CXX11_ABI "Whether to use C++11's ABI" ON)
IF(NOT ${USE_CXX11_ABI})
  MESSAGE(STATUS "Disabling C++11's ABI")
  SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 ${CMAKE_CXX_FLAGS}")
ELSE()
  MESSAGE(STATUS "Enabling C++11's ABI")
ENDIF()

#LIBRARY SETTINGS
SET(FLOAT_PRECISION double CACHE STRING "Multi-precision setting")
IF(QUADMATH_FOUND)
  SET_PROPERTY(CACHE FLOAT_PRECISION PROPERTY STRINGS double float128 mpfr_float)
ELSE(QUADMATH_FOUND)
  SET_PROPERTY(CACHE FLOAT_PRECISION PROPERTY STRINGS double mpfr_float)
ENDIF(QUADMATH_FOUND)
ADD_DEFINITIONS(-DFLOAT=${FLOAT_PRECISION})
INCLUDE("${PROJECT_SOURCE_DIR}/cmake/LoadLibs.cmake")
IF(${FLOAT_PRECISION} STREQUAL "double")
ELSE()
  ADD_DEFINITIONS(-DFORCE_ADD_DOUBLE_PRECISION)
ENDIF()

#THIRD PARTY HEADERS
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
MACRO(SEARCH_SOURCE NAME DIRPRJ)
  FILE(GLOB_RECURSE header${NAME} ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.h ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.hpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.hh)
  FILE(GLOB_RECURSE source${NAME} ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cpp ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.cc ${PROJECT_SOURCE_DIR}/${DIRPRJ}/*.c)
ENDMACRO(SEARCH_SOURCE)
SEARCH_SOURCE(Articulated Articulated)
SEARCH_SOURCE(Environment Environment)
SEARCH_SOURCE(ConvexHull ConvexHull)
SEARCH_SOURCE(Simulator Simulator)
SEARCH_SOURCE(Utils Utils)
#BUILD LIBRARY
ADD_LIBRARY(PhysicsMotion STATIC
  ${sourceArticulated} ${headerArticulated}
  ${sourceEnvironment} ${headerEnvironment}
  ${sourceConvexHull} ${headerConvexHull}
  ${sourceSimulator} ${headerSimulator}
  ${sourceUtils} ${headerUtils})
TARGET_LINK_LIBRARIES(PhysicsMotion ${ALL_LIBRARIES} ${TinyVisualizer_LIBRARIES})
SET_TARGET_PROPERTIES(PhysicsMotion PROPERTIES PUBLIC_HEADER "${headerArticulated}" "${headerEnvironment}" "${headerConvexHull}" "${headerSimulator}" "${headerUtils}")

#INSTALL LIBRARY
INSTALL(TARGETS PhysicsMotion
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${CMAKE_BUILD_TYPE}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${CMAKE_BUILD_TYPE}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PhysicsMotion)
  
#Remove visualizer
MACRO(REMOVE_VISUALIZER SRC_FILES)
  STRING(CONCAT TAR_FILES ${SRC_FILES} "NoVis")
  FOREACH(ITEM ${${SRC_FILES}})
    STRING(FIND "${ITEM}" "Visualizer" RES)
    IF(${RES} LESS 0)
      LIST(APPEND ${TAR_FILES} "${ITEM}")
    ENDIF()
  ENDFOREACH()
ENDMACRO()
 
#Pybind11-Python
IF(${FLOAT_PRECISION} STREQUAL "double")
  IF(NOT PYTHON_BINDING STREQUAL "None")
    #We do not export visualizer functions
    REMOVE_VISUALIZER(sourceArticulated)
    REMOVE_VISUALIZER(headerArticulated)
    REMOVE_VISUALIZER(sourceEnvironment)
    REMOVE_VISUALIZER(headerEnvironment)
    REMOVE_VISUALIZER(sourceConvexHull)
    REMOVE_VISUALIZER(headerConvexHull)
    REMOVE_VISUALIZER(sourceSimulator)
    REMOVE_VISUALIZER(headerSimulator)
    REMOVE_VISUALIZER(sourceUtils)
    REMOVE_VISUALIZER(headerUtils)
    #Create pybind11 library
    ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/Python/pybind11)
    PYBIND11_ADD_MODULE(pyPBAD ${PROJECT_SOURCE_DIR}/Python/pyPBAD.cpp
      ${sourceArticulatedNoVis} ${headerArticulatedNoVis}
      ${sourceEnvironmentNoVis} ${headerEnvironmentNoVis}
      ${sourceConvexHullNoVis} ${headerConvexHullNoVis}
      ${sourceSimulatorNoVis} ${headerSimulatorNoVis}
      ${sourceUtilsNoVis} ${headerUtilsNoVis})
    TARGET_LINK_LIBRARIES(pyPBAD PRIVATE ${ALL_LIBRARIES})
    #INSTALL python
    INSTALL(TARGETS pyPBAD 
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  ENDIF()
ENDIF()

#EXE
MACRO(ADD_EXE NAME)
  ADD_EXECUTABLE(${NAME} Main/${FOLDER}/${NAME}.cpp)
  IF(MSVC)
    TARGET_LINK_LIBRARIES(${NAME} PhysicsMotion ${ALL_STATIC_LIBRARIES})
  ELSE(MSVC)
    TARGET_LINK_LIBRARIES(${NAME} PhysicsMotion gomp ${ALL_STATIC_LIBRARIES})
  ENDIF(MSVC)
ENDMACRO(ADD_EXE)

#DEBUG
SET(FOLDER DebugArticulated)
SET(BUILD_DEBUG_ARTICULATED OFF)
IF(BUILD_DEBUG_ARTICULATED)
  ADD_EXE(mainArticulatedLoadBasic)
  ADD_EXE(mainArticulatedLoadURDF)
  ADD_EXE(mainPBCentroidDynamics)
  ADD_EXE(mainArticulatedUtils)
  ADD_EXE(mainContactManifold)
  ADD_EXE(mainRotationUtils)
  ADD_EXE(mainSparseUtils)
  ADD_EXE(mainEnvironment)
  ADD_EXE(mainPolynomial)
  ADD_EXE(mainFacetEdge)
  ADD_EXE(mainConvex2D)
  ADD_EXE(mainUtils)
  ADD_EXE(mainGJK)
  ADD_EXE(mainSAT)
  ADD_EXE(mainIO)
ENDIF()

SET(FOLDER DebugSimulator)
SET(BUILD_DEBUG_SIMULATOR OFF)
IF(BUILD_DEBUG_SIMULATOR)
  ADD_EXE(mainDebugPBDMatrixSolve)
  ADD_EXE(mainDebugJointLimit)
  ADD_EXE(mainDebugSoftJoint)
  #ConvHullPBD
  ADD_EXE(mainDebugCC)
  ADD_EXE(mainDebugConvHullPBD)
  ADD_EXE(mainConvHullPBDChain)
  ADD_EXE(mainConvHullPBDTwoChain)
  ADD_EXE(mainConvexHullGrasper)
  ADD_EXE(mainConvexHullChainPD)
  ADD_EXE(mainConvexHullSpiderPD)
  #PBD
  ADD_EXE(mainPBDChainKinematic)
  ADD_EXE(mainPBDBox)
  ADD_EXE(mainPBDBoxEdge)
  ADD_EXE(mainPBDBoxStack)
  ADD_EXE(mainPBDBoxCrush)
  ADD_EXE(mainPBDLimit)
  #XPBD
  ADD_EXE(mainDebugXPBD)
  ADD_EXE(mainXPBDChain)
  ADD_EXE(mainXPBDChainPD)
  ADD_EXE(mainXPBDTwoChain)
  ADD_EXE(mainXPBDReadMJCF)
  ADD_EXE(mainXPBDBoxStack)
ENDIF()
ADD_EXE(mainDebugConvHullPBDDesign)
ADD_EXE(mainDebugConvHullPBD)
ADD_EXE(mainDebugCC)
ADD_EXE(mainDebugBaryDerivative)
ADD_EXE(mainConvHullPBDTwoChain)
ADD_EXE(mainDebugJacobianDeriv)
