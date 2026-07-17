cmake_minimum_required(VERSION 3.10)

include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/project_build_setup.cmake)
pc_include_machine_config()
pc_apply_machine_compilers()

project(eigenvalue_problem LANGUAGES CXX C)

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

include_directories(${INCLUDE_DIR})

set(SOURCES
  ${SOURCE_DIR}/BSplineBasis.cpp
  ${SOURCE_DIR}/ControlPointGenerator.cpp
  ${SOURCE_DIR}/IENGenerator.cpp
  ${SOURCE_DIR}/IDGenerator.cpp
  ${SOURCE_DIR}/NURBSExtractionGenerator.cpp
  ${SOURCE_DIR}/Partition.cpp
  ${SOURCE_DIR}/FileManager.cpp
  ${SOURCE_DIR}/BernsteinBasis.cpp
  ${SOURCE_DIR}/Element.cpp
  ${SOURCE_DIR}/RefElement.cpp
  ${SOURCE_DIR}/QuadraturePoint.cpp
  ${SOURCE_DIR}/LocalAssembly.cpp
  ${SOURCE_DIR}/GlobalAssembly.cpp
  ${SOURCE_DIR}/AbscissaeGenerator.cpp
  ${SOURCE_DIR}/ElementFEM.cpp
)

add_library(eigenvalue_problem_lib ${SOURCES})

add_executable(preprocess preprocess.cpp)
target_link_libraries(preprocess eigenvalue_problem_lib)

add_executable(preprocess_fem preprocess_fem.cpp)
target_link_libraries(preprocess_fem eigenvalue_problem_lib)

add_executable(solve driver.cpp)
target_link_libraries(solve eigenvalue_problem_lib ${PETSC_LIB} ${SLEPC_LIB})

add_executable(iter_num iteration_number.cpp)
target_link_libraries(iter_num eigenvalue_problem_lib ${PETSC_LIB})

add_executable(iter_num_fem iteration_number_fempc.cpp)
target_link_libraries(iter_num_fem eigenvalue_problem_lib ${PETSC_LIB})

add_executable(mult mult.cpp)
target_link_libraries(mult eigenvalue_problem_lib ${PETSC_LIB})

pc_enable_openmp_for_targets(eigenvalue_problem_lib preprocess preprocess_fem)
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/message_compiler_setup.cmake)
