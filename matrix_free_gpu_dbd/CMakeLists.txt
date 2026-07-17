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
set(COMMON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common/include")
set(COMMON_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../common/src")

include_directories(${INCLUDE_DIR})
include_directories(${COMMON_INCLUDE_DIR}/mesh ${COMMON_INCLUDE_DIR}/mf ${COMMON_INCLUDE_DIR}/gpu)

set(SOURCES
  ${COMMON_SOURCE_DIR}/mesh/BSplineBasis.cpp
  ${COMMON_SOURCE_DIR}/mesh/ControlPointGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/IENGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/IDGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/NURBSExtractionGenerator.cpp
  ${SOURCE_DIR}/Partition.cpp
  ${SOURCE_DIR}/FileManager.cpp
  ${COMMON_SOURCE_DIR}/mesh/BernsteinBasis.cpp
  ${COMMON_SOURCE_DIR}/mesh/QuadraturePoint.cpp
  ${COMMON_SOURCE_DIR}/mf/LocalAssembly.cpp
  ${COMMON_SOURCE_DIR}/mf/ElementMFSF.cpp
  ${COMMON_SOURCE_DIR}/gpu/GlobalAssemblyDR.cpp
  ${COMMON_SOURCE_DIR}/mesh/AbscissaeGenerator.cpp
  ${COMMON_SOURCE_DIR}/mesh/ElementFEM.cpp
  ${COMMON_SOURCE_DIR}/mesh/Elem2COOGenerator.cpp
  ${SOURCE_DIR}/InvLM.cpp
)

set(CUDA_SOURCES
  ${COMMON_SOURCE_DIR}/gpu/assembly.cu
  ${SOURCE_DIR}/mult.cu
  ${COMMON_SOURCE_DIR}/gpu/GlobalAssembly.cu
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
