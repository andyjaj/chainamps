cmake_minimum_required(VERSION 3.0)
project(ChainAMPS)

option(TBB "Use tbb" ON)
option(TIMING "Timings" OFF)

if(TIMING)
  add_definitions( -DTIMING)
endif(TIMING)

#Requires c++11 compile flag
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -ffp-contract=fast -fno-signed-zeros -freciprocal-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
add_compile_options(-Wno-return-type-c-linkage -m64 -march=native -mtune=native)

#Required dependencies are BLAS/LAPACK, cxsparse and ARPACK. All should be compiled to use 64 bit integers!
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message (STATUS "BLAS:" ${BLAS_LIBRARIES})
message (STATUS "LAPACK:" ${LAPACK_LIBRARIES})

get_filename_component(LIB_HINT ${BLAS_LIBRARIES} PATH)
message (STATUS "LIBRARIES_DIR HINT:" ${LIB_HINT})
set(INCLUDE_HINT "${LIB_HINT}/../include")
message (STATUS "INCLUDE_DIR HINT:" ${INCLUDE_HINT})

find_library(ARPACK_LIBRARIES NAMES arpack HINTS ${LIB_HINT} PATHS /usr/lib64/ /usr/lib/ /usr/local/ /opt/local/)
message (STATUS "ARPACK:" ${ARPACK_LIBRARIES})
find_library(CXSPARSE_LIBRARIES cxsparse HINTS ${LIB_HINT} PATHS /usr/lib64/ /usr/lib/ /usr/local/ /opt/local/)
message (STATUS "CXSPARSE:" ${CXSPARSE_LIBRARIES})
find_path(CXSPARSE_INCLUDE_DIRS cs.h HINTS ${INCLUDE_HINT} PATHS /usr/include/ /usr/include/suitesparse/ /usr/local/ /opt/local/)
message (STATUS "CXSPARSE_INCLUDE_DIRS:" ${CXSPARSE_INCLUDE_DIRS})

set(EXTERNAL_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CXSPARSE_LIBRARIES} ${ARPACK_LIBRARIES})

include_directories(${PROJECT_BINARY_DIR} ${CXSPARSE_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/vertex_generators ${PROJECT_SOURCE_DIR}/vertex_generator/ll)
include_directories(${PROJECT_SOURCE_DIR}/optionparser)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  #needs LAPACK_ILP64 defined for openblas ilp 64
  add_definitions(-DLAPACK_ILP64)

  find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${INCLUDE_HINT} PATHS /usr/include /usr/local /opt/OpenBLAS/include)
  include_directories(${LAPACKE_INCLUDE_DIRS})

  #probably needs gfortran too
  find_library(GFORTRAN_LIBRARIES gfortran HINTS ${LIB_HINT} PATHS /usr/lib64 /usr/lib/gcc/x86_64-linux-gnu/ /usr/lib /usr/local/ /usr/)
  message (STATUS "GFORTRAN:" ${GFORTRAN_LIBRARIES})	
  set (EXTERNAL_LIBS ${EXTERNAL_LIBS} ${GFORTRAN_LIBRARIES})
endif()

if(TBB)
  find_library(TBB_LIBRARIES tbb HINT ${LIB_HINT} PATHS /usr/lib/x86_64-linux-gnu /usr/lib64 /usr/lib /usr/local /usr /opt/local)
  message (STATUS "TBB found:" ${TBB_LIBRARIES})	
  find_path(TBB_INCLUDE_DIRS tbb.h HINT ${INCLUDE_HINT} PATHS /usr/local/include/tbb /usr/include/tbb /opt/local/include/tbb)
  add_definitions( -DUSETBB)
  set (EXTERNAL_LIBS ${EXTERNAL_LIBS} ${TBB_LIBRARIES})
  include_directories(${TBB_INCLUDE_DIRS})
endif(TBB)

add_library(MATRIX_LIB STATIC sparse_interface.cpp dense_interface.cpp arpack_interface.cpp dense_matrix_functions.cpp)
add_library(MPX_LIB STATIC FiniteMPS.cpp UnitCell.cpp MPS_matrix.cpp MPO_matrix.cpp MPX_matrix.cpp MPXIndex.cpp states.cpp)
add_library(DMRG STATIC DMRG_routines.cpp iDMRG_routines.cpp fDMRG_routines.cpp)
add_library(TIME STATIC TEBD_routines.cpp)
add_library(MEAS STATIC measurement.cpp)
add_library(TWOVE STATIC TwoVertexEvolution.cpp)

set (COMMON_LIBS MPX_LIB MATRIX_LIB ${EXTERNAL_LIBS})

add_executable(iDMRG_DRV.bin iDMRG_DRV.cpp)
target_link_libraries(iDMRG_DRV.bin DMRG MEAS ${COMMON_LIBS})

add_executable(fDMRG_DRV.bin fDMRG_DRV.cpp)
target_link_libraries(fDMRG_DRV.bin DMRG MEAS ${COMMON_LIBS})

add_executable(iTEBD_DRV.bin iTEBD_DRV.cpp)
target_link_libraries(iTEBD_DRV.bin TIME MEAS ${COMMON_LIBS})

add_executable(TEBD_DRV.bin TEBD_DRV.cpp)
target_link_libraries(TEBD_DRV.bin TIME MEAS ${COMMON_LIBS})

add_executable(2VE_DRV.bin 2VE_DRV.cpp)
target_link_libraries(2VE_DRV.bin TWOVE TIME MEAS ${COMMON_LIBS})

add_executable(UNITCELL_MEASURE.bin UNITCELL_MEASURE.cpp)
target_link_libraries(UNITCELL_MEASURE.bin MEAS ${COMMON_LIBS})

add_executable(FINITE_MEASURE.bin FINITE_MEASURE.cpp)
target_link_libraries(FINITE_MEASURE.bin TIME MEAS ${COMMON_LIBS})

add_executable(STORE_MODEL.bin STORE_MODEL.cpp)
target_link_libraries(STORE_MODEL.bin ${COMMON_LIBS})

add_executable(TEBD_DYN_MEASURE.bin TEBD_DYN_MEASURE.cpp)
target_link_libraries(TEBD_DYN_MEASURE.bin TIME MEAS ${COMMON_LIBS})
