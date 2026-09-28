#include "LocalAssemblyMFSF.hpp"

void LocalAssemblyMFSF::AssemLocalLoad(ElementMFSF * const &elem,
    const std::vector<double> &eCP)
{
    std::vector<double> B1, B2, dB1, dB2, W, J, dW_dxi, dW_deta;
    std::vector<double> dxi_dx, dxi_dy, deta_dx, deta_dy;
    elem->GenerateElement(quad1, quad2, eCP, B1, B2, dB1, dB2,
        W, J, dW_dxi, dW_deta, dxi_dx, dxi_dy, deta_dx, deta_dy);

    const std::vector<double> qw1 = quad1->GetWeight();
    const std::vector<double> qw2 = quad2->GetWeight();

    // Forward tensor contraction evaluates the physical coordinates at quadrature points.
    std::vector<double> x_stage(nx*ny, 0.0);
    std::vector<double> y_stage(nx*ny, 0.0);
    for (int qx = 0; qx < nx; ++qx)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const int basis = j*nx+i;
                x_stage[qx*ny+j] += B1[qx*nx+i] * eCP[2*basis];
                y_stage[qx*ny+j] += B1[qx*nx+i] * eCP[2*basis+1];
            }
        }
    }

    std::vector<double> source(nx*ny, 0.0);
    for (int qy = 0; qy < ny; ++qy)
    {
        for (int qx = 0; qx < nx; ++qx)
        {
            const int qp = qy*nx+qx;
            double x = 0.0;
            double y = 0.0;
            for (int j = 0; j < ny; ++j)
            {
                x += B2[qy*ny+j] * x_stage[qx*ny+j];
                y += B2[qy*ny+j] * y_stage[qx*ny+j];
            }
            x /= W[qp];
            y /= W[qp];
            source[qp] = Getf(x, y) * qw1[qx] * qw2[qy] * J[qp] / W[qp];
        }
    }

    // Transpose tensor contraction applies the test basis without a local matrix.
    std::vector<double> backward_x(ny*nx, 0.0);
    for (int qy = 0; qy < ny; ++qy)
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int qx = 0; qx < nx; ++qx)
                backward_x[qy*nx+i] += B1[qx*nx+i] * source[qy*nx+qx];
        }
    }

    ResetLoad();
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int qy = 0; qy < ny; ++qy)
                Floc[j*nx+i] += B2[qy*ny+j] * backward_x[qy*nx+i];
        }
    }
}

void LocalAssemblyMFSF::LocalMatMulMF(ElementMFSF * const &elem,
    const std::vector<double> &eCP)
{
    std::vector<double> B1, B2, dB1, dB2, W, J, dW_dxi, dW_deta;
    std::vector<double> dxi_dx, dxi_dy, deta_dx, deta_dy;
    elem->GenerateElement(quad1, quad2, eCP, B1, B2, dB1, dB2,
        W, J, dW_dxi, dW_deta, dxi_dx, dxi_dy, deta_dx, deta_dy);

    const std::vector<double> qw1 = quad1->GetWeight();
    const std::vector<double> qw2 = quad2->GetWeight();

    // Forward SF: contract first in xi, then in eta to obtain u and its derivatives.
    std::vector<PetscScalar> value_x(nx*ny, 0.0);
    std::vector<PetscScalar> derivative_x_stage(nx*ny, 0.0);
    for (int qx = 0; qx < nx; ++qx)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const PetscScalar value = Floc_in[j*nx+i];
                value_x[qx*ny+j] += B1[qx*nx+i] * value;
                derivative_x_stage[qx*ny+j] += dB1[qx*nx+i] * value;
            }
        }
    }

    std::vector<PetscScalar> scale_xi(nx*ny, 0.0);
    std::vector<PetscScalar> scale_eta(nx*ny, 0.0);
    std::vector<PetscScalar> scale_value(nx*ny, 0.0);
    for (int qy = 0; qy < ny; ++qy)
    {
        for (int qx = 0; qx < nx; ++qx)
        {
            const int qp = qy*nx+qx;
            PetscScalar value = 0.0;
            PetscScalar derivative_xi = 0.0;
            PetscScalar derivative_eta = 0.0;
            for (int j = 0; j < ny; ++j)
            {
                value += B2[qy*ny+j] * value_x[qx*ny+j];
                derivative_xi += B2[qy*ny+j] * derivative_x_stage[qx*ny+j];
                derivative_eta += dB2[qy*ny+j] * value_x[qx*ny+j];
            }

            value /= W[qp];
            derivative_xi = (derivative_xi - dW_dxi[qp]*value) / W[qp];
            derivative_eta = (derivative_eta - dW_deta[qp]*value) / W[qp];

            const PetscScalar derivative_x =
                dxi_dx[qp]*derivative_xi + deta_dx[qp]*derivative_eta;
            const PetscScalar derivative_y =
                dxi_dy[qp]*derivative_xi + deta_dy[qp]*derivative_eta;
            const PetscScalar coefficient_xi =
                dxi_dx[qp]*derivative_x + dxi_dy[qp]*derivative_y;
            const PetscScalar coefficient_eta =
                deta_dx[qp]*derivative_x + deta_dy[qp]*derivative_y;
            const double quadrature_scale = J[qp] * qw1[qx] * qw2[qy];

            scale_xi[qp] = coefficient_xi * quadrature_scale / W[qp];
            scale_eta[qp] = coefficient_eta * quadrature_scale / W[qp];
            scale_value[qp] =
                -(dW_dxi[qp]*coefficient_xi + dW_deta[qp]*coefficient_eta)
                * quadrature_scale / (W[qp]*W[qp]);
        }
    }

    // Backward SF: apply the transposed derivative/value bases in reverse order.
    std::vector<PetscScalar> backward_value(ny*nx, 0.0);
    std::vector<PetscScalar> backward_eta(ny*nx, 0.0);
    for (int qy = 0; qy < ny; ++qy)
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int qx = 0; qx < nx; ++qx)
            {
                const int qp = qy*nx+qx;
                backward_value[qy*nx+i] +=
                    dB1[qx*nx+i]*scale_xi[qp] + B1[qx*nx+i]*scale_value[qp];
                backward_eta[qy*nx+i] += B1[qx*nx+i]*scale_eta[qp];
            }
        }
    }

    ResetStiffnessLoadOut();
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            PetscScalar value = 0.0;
            for (int qy = 0; qy < ny; ++qy)
            {
                value += B2[qy*ny+j]*backward_value[qy*nx+i]
                    + dB2[qy*ny+j]*backward_eta[qy*nx+i];
            }
            Floc_out[j*nx+i] = -value;
        }
    }
}
