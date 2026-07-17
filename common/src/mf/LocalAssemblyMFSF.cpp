#include "LocalAssemblyMFSF.hpp"

void LocalAssemblyMFSF::AssemLocalLoad(ElementMFSF * const &elem,
    const std::vector<double> &eCP)
{
    std::vector<double> R{};
    std::vector<double> J{};
    elem->GenerateElement(quad1, quad2, eCP, R, J);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qw1 = quad1->GetWeight();
    const std::vector<double> qw2 = quad2->GetWeight();
    std::vector<double> qw{};
    for (int i = 0; i < nqp1; ++i)
    {
        for (int j = 0; j < nqp2; ++j)
        {
            qw.push_back(qw1[j] * qw2[i]);
        }
    }
    const int nqp = nqp1 * nqp2;
    const int n = elem->GetNumLocalBasis();

    ResetLoad();

    for (int ii = 0; ii < nqp; ++ii)
    {
        double x = 0.0;
        double y = 0.0;
        for (int jj = 0; jj < n; ++jj)
        {
            x += R[ii*n+jj] * eCP[2*jj];
            y += R[ii*n+jj] * eCP[2*jj+1];
        }

        for (int jj = 0; jj < n; ++jj)
        {
            Floc[jj] += R[ii*n+jj] * Getf(x, y) * qw[ii] * J[ii];
        }
    }
}

void LocalAssemblyMFSF::LocalMatMulMF(ElementMFSF * const &elem,
    const std::vector<double> &eCP)
{
    std::vector<double> B1, B2, dB1, dB2, W, J, dW_dx, dW_dy;
    elem->GenerateElement(quad1, quad2, eCP, B1, B2, dB1, dB2, W, J, dW_dx, dW_dy);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qw1 = quad1->GetWeight();
    const std::vector<double> qw2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    ResetStiffnessLoadOut();

    for (int i = 0; i < n; ++i)
    {
        int i1 = map[i*2];
        int i2 = map[i*2+1];
        for (int j = 0; j < n; ++j)
        {
            int j1 = map[j*2];
            int j2 = map[j*2+1];
            for (int ii = 0; ii < nqp1; ++ii)
            {
                double temp1 = 0.0;
                double temp2 = 0.0;
                double temp3 = 0.0;
                double temp4 = 0.0;
                double temp5 = 0.0;
                double temp6 = 0.0;
                double temp7 = 0.0;

                for (int jj = 0; jj < nqp2; ++jj)
                {
                    temp1 += qw2[jj] * B2[nqp2*jj+i2] * B2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj]);
                    temp2 += qw2[jj] * dB2[nqp2*jj+i2] * dB2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj]);
                    temp3 += qw2[jj] * B2[nqp2*jj+i2] * B2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj] * W[ii*nqp2+jj]) * dW_dx[ii*nqp2+jj];
                    temp4 += qw2[jj] * B2[nqp2*jj+i2] * B2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj] * W[ii*nqp2+jj]) * dW_dx[ii*nqp2+jj];
                    temp5 += qw2[jj] * B2[nqp2*jj+i2] * dB2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj] * W[ii*nqp2+jj]) * dW_dy[ii*nqp2+jj];
                    temp6 += qw2[jj] * dB2[nqp2*jj+i2] * B2[nqp2*jj+j2] * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj] * W[ii*nqp2+jj]) * dW_dy[ii*nqp2+jj];
                    temp7 += qw2[jj] * B2[nqp2*jj+i2] * B2[nqp2*jj+j2] * (dW_dx[ii*nqp2+jj] * dW_dx[ii*nqp2+jj] + dW_dy[ii*nqp2+jj] 
                            * dW_dy[ii*nqp2+jj]) * J[ii*nqp2+jj] / (W[ii*nqp2+jj] * W[ii*nqp2+jj]);
                }
                Kloc[i*n+j] -= qw1[ii] * (dB1[nqp1*ii+i1] * dB1[nqp1*ii+j1] * temp1
                    + B1[nqp1*ii+i1] * B1[nqp1*ii+j1] * temp2
                    - B1[nqp1*ii+i1] * dB1[nqp1*ii+j1] * temp3
                    - dB1[nqp1*ii+i1] * B1[nqp1*ii+j1] * temp4
                    - B1[nqp1*ii+i1] * B1[nqp1*ii+j1] * temp5
                    - B1[nqp1*ii+i1] * B1[nqp1*ii+j1] * temp6
                    - B1[nqp1*ii+i1] * B1[nqp1*ii+j1] * temp7);
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Floc_out[i] += Kloc[i*n+j] * Floc_in[j];
        }
    }
}