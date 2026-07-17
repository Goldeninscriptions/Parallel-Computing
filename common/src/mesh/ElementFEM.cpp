#include "ElementFEM.hpp"

void ElementFEM::GenerateElementSingleQP(const double &xi, const double &eta,
    const std::vector<double> &eCP, 
    std::vector<double> &R, std::vector<double> &dR_dx, std::vector<double> &dR_dy,
    double &x, double &y,
    double &jacobian) const
{
    std::vector<double> N1 = GenerateBasis1DSingleQP(xi, p);
    std::vector<double> dN1 = GenerateBasisDerivative1DSingleQP(xi, p);
    std::vector<double> N2 = GenerateBasis1DSingleQP(eta, q);
    std::vector<double> dN2 = GenerateBasisDerivative1DSingleQP(eta, q);

    R.clear();
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            R.push_back(N1[i] * N2[j]);
            dN_dxi.push_back(dN1[i] * N2[j]);
            dN_deta.push_back(N1[i] * dN2[j]);
        }
    }
    
    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    jacobian = 0.0;

    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    x = 0.0;
    y = 0.0;

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            x += eCP[(j*(p+1)+i)*2] * R[j*(p+1)+i];
            y += eCP[(j*(p+1)+i)*2+1] * R[j*(p+1)+i];

            dx_dxi += eCP[(j*(p+1)+i)*2] * dN_dxi[j*(p+1)+i];
            dx_deta += eCP[(j*(p+1)+i)*2] * dN_deta[j*(p+1)+i];
            dy_dxi += eCP[(j*(p+1)+i)*2+1] * dN_dxi[j*(p+1)+i];
            dy_deta += eCP[(j*(p+1)+i)*2+1] * dN_deta[j*(p+1)+i];
        }
    }

    jacobian = dx_dxi*dy_deta - dx_deta*dy_dxi;

    dxi_dx = dy_deta/jacobian;
    dxi_dy = -dy_dxi/jacobian;
    deta_dx = -dx_deta/jacobian;
    deta_dy = dx_dxi/jacobian;

    dR_dx.clear();
    dR_dy.clear();

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            dR_dx.push_back(dxi_dx*dN_dxi[j*(p+1)+i] + dxi_dy*dN_deta[j*(p+1)+i]);
            dR_dy.push_back(deta_dx*dN_dxi[j*(p+1)+i] + deta_dy*dN_deta[j*(p+1)+i]);
        }
    }
}

void ElementFEM::GenerateRefElementSingleQP(const double &xi, const double &eta,
    std::vector<double> &N,
    std::vector<double> &dN_dxi, std::vector<double> &dN_deta) const
{
    std::vector<double> N1 = GenerateBasis1DSingleQP(xi, p);
    std::vector<double> dN1 = GenerateBasisDerivative1DSingleQP(xi, p);
    std::vector<double> N2 = GenerateBasis1DSingleQP(eta, q);
    std::vector<double> dN2 = GenerateBasisDerivative1DSingleQP(eta, q);

    N.clear();
    dN_dxi.clear();
    dN_deta.clear();

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            N.push_back(N1[i] * N2[j]);
            dN_dxi.push_back(dN1[i] * N2[j]);
            dN_deta.push_back(N1[i] * dN2[j]);
        }
    }
}

std::vector<double> ElementFEM::GenerateBasis1DSingleQP(const double &xi, const int &pp) const
{
    switch (pp)
    {
        case 1:
            return {1.0-xi, xi};
        case 2:
            return {2*xi*xi - 3*xi + 1,
                    -4*xi*xi + 4*xi,
                    2*xi*xi - xi};
        case 3:
            return {-9*xi*xi*xi + 18*xi*xi - 9*xi + 1,
                    27*xi*xi*xi - 45*xi*xi + 18*xi,
                    -27*xi*xi*xi + 36*xi*xi - 9*xi,
                    9*xi*xi*xi - 9*xi*xi + xi};
        default:
            std::cerr << "Invalid polynomial degree" << std::endl;
            exit(1);
    }
}

std::vector<double> ElementFEM::GenerateBasisDerivative1DSingleQP(const double &xi, const int &pp) const
{
    switch (pp)
    {
        case 1:
            return {-1.0, 1.0};
        case 2:
            return {4*xi - 3,
                    -8*xi + 4,
                    4*xi - 1};
        case 3:
            return {-27*xi*xi + 36*xi - 9,
                    81*xi*xi - 90*xi + 18,
                    -81*xi*xi + 72*xi - 9,
                    27*xi*xi - 18*xi + 1};
        default:
            std::cerr << "Invalid polynomial degree" << std::endl;
            exit(1);
    }
}