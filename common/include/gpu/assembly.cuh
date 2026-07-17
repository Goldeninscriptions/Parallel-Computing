#pragma once

#include <cuda_runtime.h>
#include <stdio.h>
#include <memory.cuh>

void AssembleStiffnessCUDA(const int nqp1, const int nqp2,
    const int nlocalelemx, const int nlocalelemy,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight, const int * d_IEN,
    const double * d_CP, const int * d_elem2coo,
    double * d_val);

void DirichletBCKCUDA(const int * d_dir2coo, const int dir_size, double * d_val);