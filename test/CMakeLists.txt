cmake_minimum_required(VERSION 3.18)

include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/project_build_setup.cmake)
pc_include_machine_config()
pc_apply_machine_compilers()

project(parallel_computing_test LANGUAGES CXX C CUDA)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

pc_set_default_build_flags()
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/basic_variable_log.cmake)
pc_setup_common_dependencies(
  USE_MPI
  USE_PETSC
  USE_CUDA
  MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../conf"
)

add_executable(test ghost_with_cuda.cu)
target_link_libraries(test ${PETSC_LIB} ${MPI_CXX_LIBRARIES})
target_include_directories(test PRIVATE ${PETSC_INC} ${CUDA_INC})

add_executable(test_func cuda_function_test.cu)
target_link_libraries(test_func ${PETSC_LIB} ${MPI_CXX_LIBRARIES})
target_include_directories(test_func PRIVATE ${PETSC_INC} ${CUDA_INC})

pc_configure_cuda_targets(test test_func)
include(${CMAKE_CURRENT_SOURCE_DIR}/../conf/message_compiler_setup.cmake)
