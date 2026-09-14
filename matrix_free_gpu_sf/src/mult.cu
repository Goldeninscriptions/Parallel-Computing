#ifndef P_ORDER
#define P_ORDER 3
#endif
#ifndef Q_ORDER
#define Q_ORDER 3
#endif
#define NLOC_MOD2 (((P_ORDER + 1) * (Q_ORDER + 1)) % 2)
#define MYOFFSET (NLOC_MOD2 ? 4 : 0)

#include "mult.cuh"

__device__ void compute_jacobian_basis(
    const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP,
    double &jacobian,
    double *R)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    double N1[nx];
    double dN1[nx];
    for (int i = 0; i < nx; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[ny];
    double dN2[ny];
    for (int j = 0; j < ny; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < nx; ++jj)
    {
        for (int kk = 0; kk < nx; ++kk)
        {
            N1[jj] += s_nurbs_extraction1[jj * (nx) + kk] * d_B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * (nx) + kk] * d_dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < ny; ++jj)
    {
        for (int kk = 0; kk < ny; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * (ny) + kk] * d_B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * (ny) + kk] * d_dB2[kk];
        }
        dN2[jj] /= h2;
    }

    const int nLocBas = nx * ny;
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            N[j * (nx) + i] = N1[i] * N2[j];
            w += N[j * (nx) + i];
            dN_dxi[j * (nx) + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * (nx) + i];
            dN_deta[j * (nx) + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * (nx) + i];
        }
    }

    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            R[j * (nx) + i] = N[j * (nx) + i] / w;
            dR_dxi[j * (nx) + i] = (dN_dxi[j * (nx) + i] - dw_dxi * R[j * (nx) + i]) / w;
            dR_deta[j * (nx) + i] = (dN_deta[j * (nx) + i] - dw_deta * R[j * (nx) + i]) / w;
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dx_dxi += eCP[2 * (j * (nx) + i)] * dR_dxi[j * (nx) + i];
            dx_deta += eCP[2 * (j * (nx) + i)] * dR_deta[j * (nx) + i];
            dy_dxi += eCP[2 * (j * (nx) + i) + 1] * dR_dxi[j * (nx) + i];
            dy_deta += eCP[2 * (j * (nx) + i) + 1] * dR_deta[j * (nx) + i];
        }
    }
    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;
    jacobian *= h1*h2;
}

__device__ void compute_jacobian_derivative(
    const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP, 
    double &jacobian,
    double *dR_dx, 
    double *dR_dy)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    double N1[nx];
    double dN1[nx];
    for (int i = 0; i < nx; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[ny];
    double dN2[ny];
    for (int j = 0; j < ny; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < nx; ++jj)
    {
        for (int kk = 0; kk < nx; ++kk)
        {
            N1[jj] +=  s_nurbs_extraction1[jj * nx + kk] * d_B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * nx + kk] * d_dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < ny; ++jj)
    {
        for (int kk = 0; kk < ny; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * ny + kk] * d_B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * ny + kk] * d_dB2[kk];
        }
        dN2[jj] /= h2;
    }

    const int nLocBas = nx * ny;
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            N[j * nx + i] = N1[i] * N2[j];
            w += N[j * nx + i];
            dN_dxi[j * nx + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * nx + i];
            dN_deta[j * nx + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * nx + i];
        }
    }

    double R[nLocBas];
    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            R[j * nx + i] = N[j * nx + i] / w;
            dR_dxi[j * nx + i] = (dN_dxi[j * nx + i] - dw_dxi * R[j * nx + i]) / w;
            dR_deta[j * nx + i] = (dN_deta[j * nx + i] - dw_deta * R[j * nx + i]) / w;
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dx_dxi += eCP[2 * (j * nx + i)] * dR_dxi[j * nx + i];
            dx_deta += eCP[2 * (j * nx + i)] * dR_deta[j * nx + i];
            dy_dxi += eCP[2 * (j * nx + i) + 1] * dR_dxi[j * nx + i];
            dy_deta += eCP[2 * (j * nx + i) + 1] * dR_deta[j * nx + i];
        }
    }

    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;

    dxi_dx = dy_deta / jacobian;
    dxi_dy = -dx_deta / jacobian;
    deta_dx = -dy_dxi / jacobian;
    deta_dy = dx_dxi / jacobian;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dR_dx[j * nx + i] = dxi_dx * dR_dxi[j * nx + i] + deta_dx * dR_deta[j * nx + i];
            dR_dy[j * nx + i] = dxi_dy * dR_dxi[j * nx + i] + deta_dy * dR_deta[j * nx + i];
        }
    }

    jacobian *= h1 * h2;
}

__device__ double get_force(double x, double y)
{
    return x * (1.0 - x) * y * (1.0 - y);
}

__device__ void compute_extracted_basis_1d(
    const int n,
    const double *basis_ref,
    const double *basis_der_ref,
    const double *extraction,
    const double h,
    double *basis,
    double *basis_der)
{
    for (int jj = 0; jj < n; ++jj)
    {
        basis[jj] = 0.0;
        basis_der[jj] = 0.0;
        for (int kk = 0; kk < n; ++kk)
        {
            basis[jj] += extraction[jj * n + kk] * basis_ref[kk];
            basis_der[jj] += extraction[jj * n + kk] * basis_der_ref[kk];
        }
        basis_der[jj] /= h;
    }
}

__device__ void compute_tensor_load_scale(
    const double h1,
    const double h2,
    const double *basis1_ref,
    const double *basis2_ref,
    const double *basis_der1_ref,
    const double *basis_der2_ref,
    const double *extraction1,
    const double *extraction2,
    const double *eCP,
    double &weight_sum,
    double &jacobian,
    double &x,
    double &y,
    double *basis1,
    double *basis2)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    double basis_der1[P_ORDER + 1];
    double basis_der2[Q_ORDER + 1];
    compute_extracted_basis_1d(nx, basis1_ref, basis_der1_ref, extraction1, h1, basis1, basis_der1);
    compute_extracted_basis_1d(ny, basis2_ref, basis_der2_ref, extraction2, h2, basis2, basis_der2);

    double dN_dxi[(P_ORDER + 1) * (Q_ORDER + 1)];
    double dN_deta[(P_ORDER + 1) * (Q_ORDER + 1)];
    double N[(P_ORDER + 1) * (Q_ORDER + 1)];
    weight_sum = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;
    x = 0.0;
    y = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const int idx = j * nx + i;
            N[idx] = basis1[i] * basis2[j];
            dN_dxi[idx] = basis_der1[i] * basis2[j];
            dN_deta[idx] = basis1[i] * basis_der2[j];
            weight_sum += N[idx];
            dw_dxi += dN_dxi[idx];
            dw_deta += dN_deta[idx];
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const int idx = j * nx + i;
            const double R = N[idx] / weight_sum;
            const double dR_dxi = (dN_dxi[idx] - dw_dxi * R) / weight_sum;
            const double dR_deta = (dN_deta[idx] - dw_deta * R) / weight_sum;
            x += eCP[2 * idx] * R;
            y += eCP[2 * idx + 1] * R;
            dx_dxi += eCP[2 * idx] * dR_dxi;
            dx_deta += eCP[2 * idx] * dR_deta;
            dy_dxi += eCP[2 * idx + 1] * dR_dxi;
            dy_deta += eCP[2 * idx + 1] * dR_deta;
        }
    }

    jacobian = (dx_dxi * dy_deta - dx_deta * dy_dxi) * h1 * h2;
}

__device__ void compute_tensor_geometry_data(
    const double h1,
    const double h2,
    const double *basis1_ref,
    const double *basis2_ref,
    const double *basis_der1_ref,
    const double *basis_der2_ref,
    const double *extraction1,
    const double *extraction2,
    const double *eCP,
    double &weight_sum,
    double &dw_dxi,
    double &dw_deta,
    double &jacobian,
    double &dxi_dx,
    double &dxi_dy,
    double &deta_dx,
    double &deta_dy,
    double *basis1,
    double *basis2,
    double *basis_der1,
    double *basis_der2)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;

    compute_extracted_basis_1d(nx, basis1_ref, basis_der1_ref, extraction1, h1, basis1, basis_der1);
    compute_extracted_basis_1d(ny, basis2_ref, basis_der2_ref, extraction2, h2, basis2, basis_der2);

    double N[(P_ORDER + 1) * (Q_ORDER + 1)];
    double dN_dxi[(P_ORDER + 1) * (Q_ORDER + 1)];
    double dN_deta[(P_ORDER + 1) * (Q_ORDER + 1)];
    weight_sum = 0.0;
    dw_dxi = 0.0;
    dw_deta = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const int idx = j * nx + i;
            N[idx] = basis1[i] * basis2[j];
            dN_dxi[idx] = basis_der1[i] * basis2[j];
            dN_deta[idx] = basis1[i] * basis_der2[j];
            weight_sum += N[idx];
            dw_dxi += dN_dxi[idx];
            dw_deta += dN_deta[idx];
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const int idx = j * nx + i;
            const double R = N[idx] / weight_sum;
            const double dR_dxi = (dN_dxi[idx] - dw_dxi * R) / weight_sum;
            const double dR_deta = (dN_deta[idx] - dw_deta * R) / weight_sum;
            dx_dxi += eCP[2 * idx] * dR_dxi;
            dx_deta += eCP[2 * idx] * dR_deta;
            dy_dxi += eCP[2 * idx + 1] * dR_dxi;
            dy_deta += eCP[2 * idx + 1] * dR_deta;
        }
    }

    const double jacobian_param = dx_dxi * dy_deta - dx_deta * dy_dxi;
    dxi_dx = dy_deta / jacobian_param;
    dxi_dy = -dx_deta / jacobian_param;
    deta_dx = -dy_dxi / jacobian_param;
    deta_dy = dx_dxi / jacobian_param;
    jacobian = jacobian_param * h1 * h2;
}

__global__ void AssembleKernel(
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    double *d_x_array
    )
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    extern __shared__ char shared_data[];

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += nx * ny * sizeof(int) + MYOFFSET;
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * nx * ny * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);

    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;
    const int nLocBas = nx * ny;

    for (int j = 0; j < nLocBas; ++j)
    {
        s_eID[j] = d_ID[d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
    }

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < nx * nx; ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * nx * nx + i];
    for (int i = 0; i < ny * ny; ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * ny * ny + i];

    double h1 = d_elem_size1[blockIdx.x];
    double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    int qpx = threadIdx.x;
    int qpy = threadIdx.y;
    int qp = threadIdx.y * blockDim.x + threadIdx.x;

    double B1[nx];
    double dB1[nx];
    double B2[ny];
    double dB2[ny];

    if (qp < nx*ny)
    {
        for (int i = 0; i < nx; ++i)
        {
            B1[i] = d_B1[qpx * nx + i];
            dB1[i] = d_dB1[qpx * nx + i];
        }
        for (int j = 0; j < ny; ++j)
        {
            B2[j] = d_B2[qpy * ny + j];
            dB2[j] = d_dB2[qpy * ny + j];
        }

        double jacobian;
        double R[nLocBas];

        compute_jacobian_basis(h1, h2, B1, B2, dB1, dB2,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP, jacobian, R);
        
        double x = 0.0;
        double y = 0.0;

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            x += s_eCP[2 * ii] * R[ii];
            y += s_eCP[2 * ii + 1] * R[ii];
        }

        double force = get_force(x, y);

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            double val = R[ii] * force * jacobian * s_qw[qp];
            int coo_index = s_eID[ii];
            if (coo_index >= 0)
            {
                atomicAdd(&d_x_array[coo_index], val);
            }
        }
    }
}

__global__ void SFAssembleKernel(
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    double *d_x_array
    )
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    const int nLocBas = nx * ny;
    extern __shared__ char shared_data[];

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += nLocBas * sizeof(int) + MYOFFSET;
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * nLocBas * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_basis1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_basis2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_scale = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_tmp = (double*)(shared_data + offset);

    const int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;
    const int qpx = threadIdx.x;
    const int qpy = threadIdx.y;
    const int qp = qpy * nx + qpx;

    for (int j = 0; j < nLocBas; ++j)
    {
        const int ien = d_IEN[elemIndex * nLocBas + j];
        s_eID[j] = d_ID[ien];
        s_eCP[2 * j] = d_CP[2 * ien];
        s_eCP[2 * j + 1] = d_CP[2 * ien + 1];
    }

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < nx * nx; ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * nx * nx + i];
    for (int i = 0; i < ny * ny; ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * ny * ny + i];

    const double h1 = d_elem_size1[blockIdx.x];
    const double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    if (qp < nLocBas)
    {
        double basis1[P_ORDER + 1];
        double basis2[Q_ORDER + 1];
        double weight_sum;
        double jacobian;
        double x;
        double y;

        compute_tensor_load_scale(
            h1, h2,
            d_B1 + qpx * nx, d_B2 + qpy * ny,
            d_dB1 + qpx * nx, d_dB2 + qpy * ny,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP,
            weight_sum, jacobian, x, y, basis1, basis2);

        for (int i = 0; i < nx; ++i)
            s_basis1[qpx * nx + i] = basis1[i];
        for (int j = 0; j < ny; ++j)
            s_basis2[qpy * ny + j] = basis2[j];

        const double force = get_force(x, y);
        s_scale[qp] = force * jacobian * s_qw[qp] / weight_sum;
    }

    __syncthreads();

    if (qp < nLocBas)
    {
        const int basis_x = qpx;
        const int basis_y = qpy;

        if (basis_y == 0)
        {
            for (int qy = 0; qy < ny; ++qy)
            {
                double tmp = 0.0;
                for (int qx = 0; qx < nx; ++qx)
                {
                    tmp += s_basis1[qx * nx + basis_x] * s_scale[qy * nx + qx];
                }
                s_tmp[qy * nx + basis_x] = tmp;
            }
        }

        __syncthreads();

        double val = 0.0;
        for (int qy = 0; qy < ny; ++qy)
        {
            val += s_basis2[qy * ny + basis_y] * s_tmp[qy * nx + basis_x];
        }

        const int coo_index = s_eID[basis_y * nx + basis_x];
        if (coo_index >= 0)
        {
            atomicAdd(&d_x_array[coo_index], val);
        }
    }
}

__global__ void MatrixFreeMatMultKernel(
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    const double *d_F_array_in,
    double *d_F_array_out
    )
{
    extern __shared__ char shared_data[];

    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    const int nLocBas = nx * ny;

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += nLocBas * sizeof(int) + MYOFFSET;
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * nLocBas * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *Floc_in = (double*)(shared_data + offset);
    
    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;

    for (int j = 0; j < nLocBas; ++j)
    {
        s_eID[j] = d_ID[d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
    }

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < nx * nx; ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * nx * nx + i];
    for (int i = 0; i < ny * ny; ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * ny * ny + i];

    for (int i = 0; i < nLocBas; ++i)
    {
        int coo_index = d_IEN[elemIndex * nLocBas + i];
        Floc_in[i] = d_F_array_in[coo_index];
    }

    double h1 = d_elem_size1[blockIdx.x];
    double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    double Floc_out[nLocBas];
    for (int i = 0; i < nLocBas; ++i)
    {
        Floc_out[i] = 0.0;
    }

    int qpx = threadIdx.x;
    int qpy = threadIdx.y;
    int qp = threadIdx.y * blockDim.x + threadIdx.x;

    double B1[nx];
    double dB1[nx];
    double B2[ny];
    double dB2[ny];

    if (qp < nx * ny)
    {
        for (int i = 0; i < nx; ++i)
        {
            B1[i] = d_B1[qpx * nx + i];
            dB1[i] = d_dB1[qpx * nx + i];
        }
        for (int j = 0; j < ny; ++j)
        {
            B2[j] = d_B2[qpy * ny + j];
            dB2[j] = d_dB2[qpy * ny + j];
        }

        double jacobian;
        double dR_dx[nLocBas];
        double dR_dy[nLocBas];

        compute_jacobian_derivative(h1, h2, B1, B2, dB1, dB2,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP, jacobian, dR_dx, dR_dy);
        
        double temp_x = 0.0;
        double temp_y = 0.0;

        for (int jj = 0; jj < nLocBas; ++jj)
        {
            temp_x += dR_dx[jj] * Floc_in[jj];
            temp_y += dR_dy[jj] * Floc_in[jj];
        }

        temp_x *= -s_qw[qp]*jacobian;
        temp_y *= -s_qw[qp]*jacobian;

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            Floc_out[ii] += (dR_dx[ii] * temp_x + dR_dy[ii] * temp_y);
        }

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            int coo_index = s_eID[ii];
            if (coo_index >= 0)
            {
                atomicAdd(&d_F_array_out[coo_index], Floc_out[ii]);
            }
        }
    }
}

__global__ void SFMatrixFreeMatMultKernel(
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    const double *d_F_array_in,
    double *d_F_array_out
    )
{
    extern __shared__ char shared_data[];

    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    const int nLocBas = nx * ny;

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += nLocBas * sizeof(int) + MYOFFSET;
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * nLocBas * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *Floc_in = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_basis1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_basis2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_basis_der1 = (double*)(shared_data + offset);
    offset += nx * nx * sizeof(double);
    double *s_basis_der2 = (double*)(shared_data + offset);
    offset += ny * ny * sizeof(double);
    double *s_scale1 = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_scale2 = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_scale3 = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_tmp1 = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *s_tmp2 = (double*)(shared_data + offset);

    const int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;
    const int qpx = threadIdx.x;
    const int qpy = threadIdx.y;
    const int qp = qpy * nx + qpx;

    for (int j = 0; j < nLocBas; ++j)
    {
        const int ien = d_IEN[elemIndex * nLocBas + j];
        s_eID[j] = d_ID[ien];
        s_eCP[2 * j] = d_CP[2 * ien];
        s_eCP[2 * j + 1] = d_CP[2 * ien + 1];
        Floc_in[j] = d_F_array_in[ien];
    }

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < nx * nx; ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * nx * nx + i];
    for (int i = 0; i < ny * ny; ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * ny * ny + i];

    const double h1 = d_elem_size1[blockIdx.x];
    const double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    if (qp < nLocBas)
    {
        double basis1[P_ORDER + 1];
        double basis2[Q_ORDER + 1];
        double basis_der1[P_ORDER + 1];
        double basis_der2[Q_ORDER + 1];
        double weight_sum;
        double dw_dxi;
        double dw_deta;
        double jacobian;
        double dxi_dx;
        double dxi_dy;
        double deta_dx;
        double deta_dy;

        compute_tensor_geometry_data(
            h1, h2,
            d_B1 + qpx * nx, d_B2 + qpy * ny,
            d_dB1 + qpx * nx, d_dB2 + qpy * ny,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP,
            weight_sum, dw_dxi, dw_deta, jacobian,
            dxi_dx, dxi_dy, deta_dx, deta_dy,
            basis1, basis2, basis_der1, basis_der2);

        for (int i = 0; i < nx; ++i)
        {
            s_basis1[qpx * nx + i] = basis1[i];
            s_basis_der1[qpx * nx + i] = basis_der1[i];
        }
        for (int j = 0; j < ny; ++j)
        {
            s_basis2[qpy * ny + j] = basis2[j];
            s_basis_der2[qpy * ny + j] = basis_der2[j];
        }

        double A[Q_ORDER + 1];
        double Ad[Q_ORDER + 1];
        for (int j = 0; j < ny; ++j)
        {
            A[j] = 0.0;
            Ad[j] = 0.0;
            for (int i = 0; i < nx; ++i)
            {
                const int idx = j * nx + i;
                A[j] += basis1[i] * Floc_in[idx];
                Ad[j] += basis_der1[i] * Floc_in[idx];
            }
        }

        double uval = 0.0;
        double dudxi = 0.0;
        double dudeta = 0.0;
        for (int j = 0; j < ny; ++j)
        {
            uval += basis2[j] * A[j];
            dudxi += basis2[j] * Ad[j];
            dudeta += basis_der2[j] * A[j];
        }
        uval /= weight_sum;
        dudxi = (dudxi - dw_dxi * uval) / weight_sum;
        dudeta = (dudeta - dw_deta * uval) / weight_sum;

        const double dudx = dxi_dx * dudxi + deta_dx * dudeta;
        const double dudy = dxi_dy * dudxi + deta_dy * dudeta;

        const double coef_xi = dxi_dx * dudx + dxi_dy * dudy;
        const double coef_eta = deta_dx * dudx + deta_dy * dudy;

        const double scale = jacobian * s_qw[qp];
        s_scale1[qp] = coef_xi * scale / weight_sum;
        s_scale2[qp] = coef_eta * scale / weight_sum;
        s_scale3[qp] = -(dw_dxi * coef_xi + dw_deta * coef_eta) * scale / (weight_sum * weight_sum);
    }

    __syncthreads();

    if (qp < nLocBas)
    {
        const int basis_x = qpx;
        const int basis_y = qpy;

        if (basis_y == 0)
        {
            for (int qy = 0; qy < ny; ++qy)
            {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                for (int qx = 0; qx < nx; ++qx)
                {
                    tmp1 += s_basis_der1[qx * nx + basis_x] * s_scale1[qy * nx + qx]
                          + s_basis1[qx * nx + basis_x] * s_scale3[qy * nx + qx];
                    tmp2 += s_basis1[qx * nx + basis_x] * s_scale2[qy * nx + qx];
                }
                s_tmp1[qy * nx + basis_x] = tmp1;
                s_tmp2[qy * nx + basis_x] = tmp2;
            }
        }

        __syncthreads();

        double val = 0.0;
        for (int qy = 0; qy < ny; ++qy)
        {
            val += s_basis2[qy * ny + basis_y] * s_tmp1[qy * nx + basis_x]
                 + s_basis_der2[qy * ny + basis_y] * s_tmp2[qy * nx + basis_x];
        }

        const int coo_index = s_eID[basis_y * nx + basis_x];
        if (coo_index >= 0)
        {
            atomicAdd(&d_F_array_out[coo_index], -val);
        }
    }
}

__global__ void DirichletBCKernel(const int * d_Dir, const int dirsize, double * d_val, double value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dirsize)
    {
        int coo_index = d_Dir[idx];
        if (coo_index >= 0)
        {
            d_val[coo_index] = value;
        }
    }
}

void AssembleLoadCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2, double * d_F_array)
{
    if (p != P_ORDER || q != Q_ORDER)
    {
        printf("Error: p and q must match the defined P_ORDER and Q_ORDER.\n");
        exit(EXIT_FAILURE);
    }

    size_t shared_size = (p + 1) * (q + 1) * sizeof(int)
                + 2 * (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + MYOFFSET;

    size_t sf_shared_size = shared_size
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double);

    SFAssembleKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(p+1, q+1), sf_shared_size>>>(
        d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2, d_F_array);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in AssembleLoadCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void MatrixFreeMatMultCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2,
    const double * d_F_array_in, double * d_F_array_out)
{
    if (p != P_ORDER || q != Q_ORDER)
    {
        printf("Error: p and q must match the defined P_ORDER and Q_ORDER.\n");
        exit(EXIT_FAILURE);
    }

    size_t shared_size = (p + 1) * (q + 1) * sizeof(int)
                + 2 * (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + MYOFFSET;

    size_t sf_shared_size = shared_size
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double);

    SFMatrixFreeMatMultKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(p+1, q+1), sf_shared_size>>>(
        d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_F_array_in, d_F_array_out);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in MatrixFreeMatMultCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void DirichletBCCUDA(const int * d_Dir, const int dirsize, double * d_x_array, double value)
{
    int blocksize = 256;
    int nblocks = (dirsize + blocksize - 1) / blocksize;

    DirichletBCKernel<<<nblocks, blocksize>>>(d_Dir, dirsize, d_x_array, value);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in DirichletBCCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}
