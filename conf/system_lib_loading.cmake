# Machine-local external library setup for Parallel-Computing.
# This file is intended to be edited per machine and loaded by
# Parallel-Computing/conf/project_build_setup.cmake.

# ========================================================
# Core library roots
# ========================================================
set(PC_PETSC_DIR /home/xuanming/lib/petsc-3.24.0-opt)
set(PC_PETSC_ARCH .)
set(PC_SLEPC_DIR /home/xuanming/lib/slepc-3.24.0)
set(PC_EIGEN3_INCLUDE_DIR /home/xuanming/lib/eigen-3.4.1)

# ========================================================
# Compiler setup
# ========================================================
set(PC_MPI_C_COMPILER /home/xuanming/lib/mpich-4.3.2/bin/mpicc)
set(PC_MPI_CXX_COMPILER /home/xuanming/lib/mpich-4.3.2/bin/mpicxx)
set(PC_MPI_INCLUDE_DIRS /home/xuanming/lib/mpich-4.3.2/include)
set(PC_MPI_LIBRARIES /home/xuanming/lib/mpich-4.3.2/lib/libmpi.so)

# ========================================================
# CUDA setup
# ========================================================
set(PC_CUDA_ROOT /usr/local/cuda-13)
set(PC_CUDA_COMPILER /usr/local/cuda-13/bin/nvcc)
set(PC_CUDA_ARCHITECTURES 75)

# ========================================================
# Optional build knobs
# ========================================================
set(PC_OPENMP_CXX_INCLUDE_DIR "/usr/lib/gcc/x86_64-linux-gnu/13/include")
set(PC_VERBOSE_MAKEFILE OFF)
