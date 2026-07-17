#include "LocalAssemblyMF.hpp"

void LocalAssemblyMF::AssemLocalLoad(ElementMF * const &elem,
    const std::vector<double> &eCP)
{
    elem->GenerateElement(quad1, quad2, eCP);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const int nqp = nqp1 * nqp2;
    const int n = elem->GetNumLocalBasis();

    ResetStiffnessLoad();

    for (int ii = 0; ii < nqp; ++ii)
    {
        double x = 0.0;
        double y = 0.0;
        for (int jj = 0; jj < n; ++jj)
        {
            x += elem->get_R(ii, jj) * eCP[2*jj];
            y += elem->get_R(ii, jj) * eCP[2*jj+1];
        }

        for (int jj = 0; jj < n; ++jj)
        {
            Floc[jj] += elem->get_R(ii, jj) * Getf(x, y) * elem->get_JxW(ii);
        }
    }
}

void LocalAssemblyMF::LocalMatMulMF(ElementMF * const &elem,
    const std::vector<double> &eCP)
{
    elem->GenerateElement(quad1, quad2, eCP);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const int nqp = nqp1 * nqp2;
    const int n = elem->GetNumLocalBasis();

    std::vector<double> temp_x(nqp, 0.0);
    std::vector<double> temp_y(nqp, 0.0);

    ResetLocalLoad();

    for (int ii = 0; ii < nqp; ++ii)
    {
        for (int jj = 0; jj < n; ++jj)
        {
            temp_x[ii] += elem->get_dR_dx(ii, jj) * Floc_in[jj];
            temp_y[ii] += elem->get_dR_dy(ii, jj) * Floc_in[jj];
        }
    }

    for (int ii = 0; ii < nqp; ++ii)
    {
        temp_x[ii] *= -elem->get_JxW(ii);
        temp_y[ii] *= -elem->get_JxW(ii);
    }

    for (int jj = 0; jj < nqp; ++jj)
    {
        for (int ii = 0; ii < n; ++ii)
        {
            Floc_out[ii] += elem->get_dR_dx(jj, ii) * temp_x[jj];
            Floc_out[ii] += elem->get_dR_dy(jj, ii) * temp_y[jj]; 
        }
    }
}