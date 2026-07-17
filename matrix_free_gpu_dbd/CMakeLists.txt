cmake_minimum_required(VERSION 3.18)

include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/project_build_setup.cmake)
pc_include_machine_config()
pc_apply_machine_compilers()

project(matrix_free_gpu_dbd LANGUAGES CXX C CUDA)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

pc_set_default_build_flags()
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/basic_variable_log.cmake)
pc_setup_common_dependencies(
  USE_OPENMP
  USE_MPI
  USE_PETSC
  USE_CUDA
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
  ${SOURCE_DIR}/QuadraturePoint.cpp
  ${SOURCE_DIR}/LocalAssembly.cpp
  ${SOURCE_DIR}/GlobalAssemblyDR.cpp
  ${SOURCE_DIR}/AbscissaeGenerator.cpp
  ${SOURCE_DIR}/ElementFEM.cpp
  ${SOURCE_DIR}/Elem2COOGenerator.cpp
  ${SOURCE_DIR}/InvLM.cpp
)

set(CUDA_SOURCES
  ${SOURCE_DIR}/assembly.cu
  ${SOURCE_DIR}/mult.cu
  ${SOURCE_DIR}/GlobalAssembly.cu
  ${SOURCE_DIR}/GlobalAssemblyMF.cu
)

set(P_ORDER 3 CACHE STRING "Polynomial order p")
set(Q_ORDER 3 CACHE STRING "Polynomial order q")

add_library(matrix_free_gpu_dbd_lib ${SOURCES})
add_library(matrix_free_gpu_dbd_cuda_lib ${CUDA_SOURCES})
target_link_libraries(matrix_free_gpu_dbd_cuda_lib PUBLIC matrix_free_gpu_dbd_lib)
target_include_directories(matrix_free_gpu_dbd_cuda_lib PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_definitions(matrix_free_gpu_dbd_cuda_lib PRIVATE P_ORDER=${P_ORDER} Q_ORDER=${Q_ORDER})

add_executable(preprocess preprocess.cpp)
target_link_libraries(preprocess matrix_free_gpu_dbd_lib)

add_executable(preprocess_fem preprocess_fem.cpp)
target_link_libraries(preprocess_fem matrix_free_gpu_dbd_lib)

add_executable(dry_run dry_run.cpp)
target_link_libraries(dry_run matrix_free_gpu_dbd_lib ${PETSC_LIB})

add_executable(invlmtest InvLMtest.cpp)
target_link_libraries(invlmtest matrix_free_gpu_dbd_lib ${PETSC_LIB})

add_executable(multtest multtest.cu)
target_link_libraries(multtest matrix_free_gpu_dbd_lib matrix_free_gpu_dbd_cuda_lib ${PETSC_LIB})
target_include_directories(multtest PRIVATE ${PETSC_INC} ${CUDA_INC} ${MPI_C_INCLUDE_PATH})
target_compile_definitions(multtest PRIVATE P_ORDER=${P_ORDER} Q_ORDER=${Q_ORDER})

add_executable(matrixfreegpu matrixfreegpu.cu)
target_link_libraries(matrixfreegpu matrix_free_gpu_dbd_lib matrix_free_gpu_dbd_cuda_lib ${PETSC_LIB})
target_include_directories(matrixfreegpu PRIVATE ${PETSC_INC} ${CUDA_INC} ${MPI_C_INCLUDE_PATH})
target_compile_definitions(matrixfreegpu PRIVATE P_ORDER=${P_ORDER} Q_ORDER=${Q_ORDER})

pc_enable_openmp_for_targets(matrix_free_gpu_dbd_lib preprocess preprocess_fem)
pc_configure_cuda_targets(matrix_free_gpu_dbd_cuda_lib multtest matrixfreegpu)

include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/message_compiler_setup.cmake)
