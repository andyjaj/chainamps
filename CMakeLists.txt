cmake_minimum_required(VERSION 2.6)
project(ChainAMPS)

option(TBB "Use tbb" ON)

#Requires c++11 compile flag
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_BUILD_TYPE Release)
add_compile_options(-Wno-return-type-c-linkage -m64 -march=native -mtune=native)

#Required dependencies are BLAS/LAPACK, cxsparse and ARPACK. All should be compiled to use 64 bit integers!
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message (STATUS "BLAS found:" ${BLAS_LIBRARIES})
message (STATUS "LAPACK found:" ${LAPACK_LIBRARIES})

get_filename_component(LIB_HINT ${BLAS_LIBRARIES} PATH)
#message (STATUS "LIBRARIES_DIR:" ${LIB_HINT})

set(INCLUDE_HINT "${LIB_HINT}/../include")
#message (STATUS "INCLUDE_DIR:" ${INCLUDE_HINT})

find_library(ARPACK_LIBRARIES NAMES arpack HINTS ${LIB_HINT} PATHS /usr/lib64/ /usr/local/ /opt/local/)
message (STATUS "ARPACK found:" ${ARPACK_LIBRARIES})
find_library(CXSPARSE_LIBRARIES cxsparse HINTS ${LIB_HINT} PATHS /usr/lib64/ /usr/local/ /opt/local/)
message (STATUS "CXSPARSE found:" ${CXSPARSE_LIBRARIES})
find_path(CXSPARSE_INCLUDE_DIRS cs.h HINTS ${INCLUDE_HINT} PATHS /usr/ /usr/local/ /opt/local/)

set(EXTERNAL_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CXSPARSE_LIBRARIES} ${ARPACK_LIBRARIES})

include_directories(${PROJECT_BINARY_DIR} ${CXSPARSE_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/vertex_generators ${PROJECT_SOURCE_DIR}/vertex_generator/ll)
include_directories(${PROJECT_SOURCE_DIR}/optionparser)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions(-DLAPACK_ILP64)
  #probably need gfortran too
  find_library(GFORTRAN_LIBRARIES gfortran HINTS ${LIB_HINT} PATHS /usr/lib64/ /usr/local/ /usr/)
  message (STATUS "GFORTRAN found:" ${GFORTRAN_LIBRARIES})	
  set (EXTERNAL_LIBS ${EXTERNAL_LIBS} ${GFORTRAN_LIBRARIES})
endif()

if(TBB)
  find_library(TBB_LIBRARIES tbb HINT ${LIB_HINT} PATHS /usr/lib64/ /usr/local/ /usr/ /opt/local/)
  message (STATUS "TBB found:" ${TBB_LIBRARIES})	
  find_path(TBB_INCLUDE_DIRS tbb.h HINT ${INCLUDE_HINT} PATHS /usr/local/include/tbb/ /usr/include/tbb/ /opt/local/include/tbb/)
  add_definitions( -DUSETBB)
  set (EXTERNAL_LIBS ${EXTERNAL_LIBS} ${TBB_LIBRARIES})
  include_directories(${TBB_INCLUDE_DIRS})
endif(TBB)

add_library(DMRG STATIC DMRG_routines.cpp measurement.cpp MPX.cpp sparse_interface.cpp dense_interface.cpp arpack_interface.cpp dense_matrix_functions.cpp states.cpp)
add_library(TIME STATIC TEBD_routines.cpp measurement.cpp MPX.cpp sparse_interface.cpp dense_interface.cpp arpack_interface.cpp dense_matrix_functions.cpp states.cpp)

add_executable(iDMRG_DRV.bin iDMRG_DRV.cpp)
target_link_libraries(iDMRG_DRV.bin DMRG ${EXTERNAL_LIBS})

add_executable(fDMRG_DRV.bin fDMRG_DRV.cpp)
target_link_libraries(fDMRG_DRV.bin DMRG ${EXTERNAL_LIBS})

add_executable(iTEBD_DRV.bin iTEBD_DRV.cpp)
target_link_libraries(iTEBD_DRV.bin TIME ${EXTERNAL_LIBS})

add_executable(TEBD_DRV.bin TEBD_DRV.cpp)
target_link_libraries(TEBD_DRV.bin TIME ${EXTERNAL_LIBS})
