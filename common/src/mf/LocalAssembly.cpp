#include "LocalAssembly.hpp"

void LocalAssembly::AssemLocalStiffnessLoad(const Element * const &elem,
    const std::vector<double> &eCP)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> R{};
    std::vector<double> dR_dx{};
    std::vector<double> dR_dy{};
    double jacobian;
    double xx, yy;

    ResetStiffnessLoad();

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateElementSingleQP(qp1[ii], qp2[jj], eCP, R, dR_dx, dR_dy, xx, yy, jacobian);

            double J_W = w1[ii]*w2[jj]*jacobian;

            for (int kk = 0; kk < n; ++kk)
            {
                for (int ll = 0; ll < n; ++ll)
                {
                    Kloc[kk*n+ll] -= J_W * (dR_dx[kk]*dR_dx[ll] + dR_dy[kk]*dR_dy[ll]);
                }
                Floc[kk] += J_W * R[kk] * Getf(xx, yy);
            }
        }
    }
}

void LocalAssembly::AssemLocalStiffnessLoad(const ElementFEM * const &elem,
    const std::vector<double> &eCP)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> R{};
    std::vector<double> dR_dx{};
    std::vector<double> dR_dy{};
    double jacobian;
    double xx, yy;

    ResetStiffnessLoad();

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateElementSingleQP(qp1[ii], qp2[jj], eCP, R, dR_dx, dR_dy, xx, yy, jacobian);

            double J_W = w1[ii]*w2[jj]*jacobian;

            for (int kk = 0; kk < n; ++kk)
            {
                for (int ll = 0; ll < n; ++ll)
                {
                    Kloc[kk*n+ll] -= J_W * (dR_dx[kk]*dR_dx[ll] + dR_dy[kk]*dR_dy[ll]);
                }
                Floc[kk] += J_W * R[kk] * Getf(xx, yy);
            }
        }
    }
}