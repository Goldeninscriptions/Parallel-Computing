cmake_minimum_required(VERSION 3.10)

include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/project_build_setup.cmake)
pc_include_machine_config()
pc_apply_machine_compilers()

project(matrix_free LANGUAGES CXX C)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

pc_set_default_build_flags()
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/basic_variable_log.cmake)
pc_setup_common_dependencies(
  USE_OPENMP
  USE_PETSC
  USE_SLEPC
  MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../conf"
)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(COMMON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common/include")
set(COMMON_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common/src")

include_directories(${INCLUDE_DIR})
include_directories(${COMMON_INCLUDE_DIR}/mesh ${COMMON_INCLUDE_DIR}/mf)

set(SOURCES
  ${COMMON_SOURCE_DIR}/mesh/BSplineBasis.cpp
  ${COMMON_SOURCE_DIR}/mesh/ControlPointGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/IENGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/IDGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/NURBSExtractionGenerator.cpp
  ${SOURCE_DIR}/Partition.cpp
  ${SOURCE_DIR}/FileManager.cpp
  ${SOURCE_DIR}/BernsteinBasis.cpp
  ${COMMON_SOURCE_DIR}/mesh/Element.cpp
  ${SOURCE_DIR}/RefElement.cpp
  ${COMMON_SOURCE_DIR}/mesh/QuadraturePoint.cpp
  ${COMMON_SOURCE_DIR}/mf/LocalAssembly.cpp
  ${SOURCE_DIR}/GlobalAssembly.cpp
  ${COMMON_SOURCE_DIR}/mesh/AbscissaeGenerator.cpp
  ${SOURCE_DIR}/ElementFEM.cpp
  ${SOURCE_DIR}/GlobalAssemblyMF.cpp
  ${SOURCE_DIR}/LocalAssemblyMF.cpp
  ${COMMON_SOURCE_DIR}/mf/LocalAssemblyMFSF.cpp
  ${SOURCE_DIR}/ElementMF.cpp
  ${COMMON_SOURCE_DIR}/mf/ElementMFSF.cpp
)

add_library(matrix_free_lib ${SOURCES})

add_executable(preprocess preprocess.cpp)
target_link_libraries(preprocess matrix_free_lib)

add_executable(preprocess_fem preprocess_fem.cpp)
target_link_libraries(preprocess_fem matrix_free_lib)

add_executable(iter_num_fem iteration_number_fempc.cpp)
target_link_libraries(iter_num_fem matrix_free_lib ${PETSC_LIB} ${SLEPC_LIB})

add_executable(iter_num_fem_mf iteration_number_fempc_mf.cpp)
target_link_libraries(iter_num_fem_mf matrix_free_lib ${PETSC_LIB})

add_executable(iter_num_fem_mfsf iteration_number_fempc_mfsf.cpp)
target_link_libraries(iter_num_fem_mfsf matrix_free_lib ${PETSC_LIB})

add_executable(mult_mf mult_mf.cpp)
target_link_libraries(mult_mf matrix_free_lib ${PETSC_LIB})

add_executable(mult_sf mult_sf.cpp)
target_link_libraries(mult_sf matrix_free_lib ${PETSC_LIB})

pc_enable_openmp_for_targets(matrix_free_lib preprocess preprocess_fem)
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/message_compiler_setup.cmake)
