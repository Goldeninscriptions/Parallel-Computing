#include "assembly.cuh"

__device__ void compute_jacobian_basis_derivative(
    const int nLocBas, const double * dN_dxi, const double * dN_deta,
    const double * eCP, double &jacobian, double * dR_dx, double * dR_dy)
{
    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    jacobian = 0.0;

    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    for (int ii = 0; ii < nLocBas; ++ii)
    {
        dx_dxi += eCP[ii * 2] * dN_dxi[ii];
        dx_deta += eCP[ii * 2] * dN_deta[ii];
        dy_dxi += eCP[ii * 2 + 1] * dN_dxi[ii];
        dy_deta += eCP[ii * 2 + 1] * dN_deta[ii];
    }

    jacobian = dx_dxi*dy_deta - dx_deta*dy_dxi;

    dxi_dx = dy_deta/jacobian;
    dxi_dy = -dy_dxi/jacobian;
    deta_dx = -dx_deta/jacobian;
    deta_dy = dx_dxi/jacobian;

    for (int ii = 0; ii < nLocBas; ++ii)
    {
        dR_dx[ii] = dxi_dx * dN_dxi[ii] + dxi_dy * dN_deta[ii];
        dR_dy[ii] = deta_dx * dN_dxi[ii] + deta_dy * dN_deta[ii];
    }
}

__global__ void AssembleStiffnessKernel(const int nqp, const int nelem,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight, const int * d_IEN,
    const double * d_CP, const int * elem2coo, double * d_val)
{
    int elemIndex = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIndex >= nelem) return;

    double eCP[2 * 4];

    for (int j = 0; j < 4; ++j)
    {
        int ien = d_IEN[elemIndex * 4 + j];
        eCP[2 * j]     = d_CP[2 * ien];
        eCP[2 * j + 1] = d_CP[2 * ien + 1];
    }

    double dN_dxi_q[4];
    double dN_deta_q[4];

    for (int qp = 0; qp < nqp; ++qp)
    {
        for (int i = 0; i < 4; ++i)
        {
            dN_dxi_q[i] = d_dN_dxi[qp * 4 + i];
            dN_deta_q[i] = d_dN_deta[qp * 4 + i];
        }

        double jacobian;
        double dR_dx[4];
        double dR_dy[4];

        compute_jacobian_basis_derivative(4, dN_dxi_q, dN_deta_q, eCP, jacobian, dR_dx, dR_dy);

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                int elem2coo_index = elem2coo[elemIndex * 16 + i * 4 + j];
                if (elem2coo_index >= 0)
                {
                    double val = -(dR_dx[i] * dR_dx[j] + dR_dy[i] * dR_dy[j]) * jacobian * d_weight[qp];
                    atomicAdd(&d_val[elem2coo_index], val);
                }
            }
        }
    }
}

__global__ void DirichletBCKKernel(const int * dir2coo, const int dirsize, double * d_val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= dirsize) return;

    int coo_index = dir2coo[idx];
    d_val[coo_index] = 1.0;
}

void AssembleStiffnessCUDA(const int nqp1, const int nqp2,
    const int nlocalelemx, const int nlocalelemy,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight, const int * d_IEN,
    const double * d_CP, const int * d_elem2coo,
    double * d_val)
{
    int nelem = nlocalelemx * nlocalelemy;
    int blocksize = 256;
    int gridsize = (nelem + blocksize - 1) / blocksize;

    AssembleStiffnessKernel<<<gridsize, blocksize>>>(nqp1*nqp2,
        nlocalelemx*nlocalelemy,
        d_N, d_dN_dxi, d_dN_deta,
        d_weight, d_IEN, d_CP, d_elem2coo, d_val);
    
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error launching kernel: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void DirichletBCKCUDA(const int * d_dir2coo, const int dir_size, double * d_val)
{
    int blockSize = 256;
    int numBlocks = (dir_size + blockSize - 1) / blockSize;

    DirichletBCKKernel<<<numBlocks, blockSize>>>(d_dir2coo, dir_size, d_val);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error launching DirichletBC kernel: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}