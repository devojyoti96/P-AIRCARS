cmake_minimum_required(VERSION 2.8.8)

project(mwa-reduce)

# Add CMake subdirectory to module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CMAKE_CXX_FLAGS "-O3 -ggdb -Wvla -Wall -DNDEBUG -march=native -std=c++0x -pthread")

find_package(HDF5 COMPONENTS C CXX REQUIRED)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits scimath scimath_f)
set(CASA_LIBS ${CASACORE_LIBRARIES})

find_package(CFITSIO REQUIRED)

find_package(Threads REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
find_package(Boost COMPONENTS date_time filesystem system thread REQUIRED)
find_library(FITSIO_LIB cfitsio REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
find_path(GTS_INCLUDE_DIR NAMES gts.h)
find_path(GLIB_INCLUDE_DIR NAMES glib.h PATH_SUFFIXES glib-2.0)
find_package(LibXml2)

find_path(CASA_INCLUDE_DIR NAMES ms/MeasurementSets/MeasurementSet.h PATHS /usr/local/include/casacore PATH_SUFFIXES casacore)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CFITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})

set(BEAM_LIBS mwareducebeam ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  include_directories(${GSL_INCLUDE_DIR})
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

if(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)
  include_directories(${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR})
endif(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)

set(CMAKE_REQUIRED_FLAGS "-std=c++0x")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/uvector_11.h\"
int main(int argc, char* argv[])
{
return 0;
}" UVECTOR_COMPILES)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/lane_11.h\"
int main(int argc, char* argv[])
{
return 0;
}" LANE11_COMPILES)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

link_libraries(${CMAKE_THREAD_LIBS_INIT})

add_library(mwareducebeam fitsiochecker.cpp beam/joneslookupdipole.cpp beam/system.cpp beam/tilebeambase.cpp beam/tilebeam2013.cpp beam/tilebeam2014.cpp beam/tilebeam2016.cpp beam/beam2016implementation.cpp beam/tileimpedance.cpp beam/lnaimpedance.cpp alglib/ap.cpp alglib/alglibinternal.cpp alglib/alglibmisc.cpp alglib/integration.cpp alglib/interpolation.cpp alglib/linalg.cpp alglib/optimization.cpp alglib/solvers.cpp alglib/specialfunctions.cpp)
set(BEAMFILES "")

add_executable(applysolutions applysolutions.cpp)
add_executable(calibrate calibrate.cpp calibrator.cpp beamevaluator.cpp calibrationmethod.cpp fftwmanager.cpp model/model.cpp predicter.cpp mspredicter.cpp nlplfitter.cpp progressbar.cpp ${BEAMFILES})
add_executable(solutiontool solutiontool.cpp)

target_link_libraries(applysolutions ${CASA_LIBS})
target_link_libraries(calibrate ${BEAM_LIBS} ${CASA_LIBS} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${CFITSIO_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
add_executable(beam2016test EXCLUDE_FROM_ALL beam/beam2016implementation.cpp beam/beam2016test.cpp beam/system.cpp )
target_link_libraries(beam2016test ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${HDF5_LIBRARIES})

if(NOT GSL_CBLAS_LIB)
  message("GSL CBLAS lib was not found. GSL needs CBLAS: disabling GSL.")
endif(NOT GSL_CBLAS_LIB)

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  add_definitions(-DHAVE_GSL)
  
else(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  message("GNU scientific library (gsl) not found: some features will be disabled (this is not an error).")
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

if(UVECTOR_COMPILES)
  add_definitions(-DHAVE_UVECTOR)
else(UVECTOR_COMPILES)
  message("ao::uvector does not compile: will use slower std::vector (this is not an error).")
endif(UVECTOR_COMPILES)
if(LANE11_COMPILES)
  add_definitions(-DHAVE_LANE11)
else(LANE11_COMPILES)
  message("ao::lane does not compile with C++11 features: will only use C++03 (this is not an error).")
endif(LANE11_COMPILES)

