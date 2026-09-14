#pragma once

#include <cuda_runtime.h>
#include <stdio.h>
#include <memory.cuh>

void AssembleLoadCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2, double * d_x_array);

void MatrixFreeMatMultCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2,
    const double * d_F_array_in, double * d_F_array_out);

void DirichletBCCUDA(const int * d_Dir,
    const int dir_size, double * d_x_array, const double value);