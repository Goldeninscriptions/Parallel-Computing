#include "RefElement.hpp"

std::vector<double> RefElement::GenerateBasis1DSingleQP(const BernsteinBasis * const &bern,
    const std::vector<double> &extraction,
    const double &xi)
{
    const int p = bern->GetDegree();
    std::vector<double> basis{};
    const std::vector<double> BB = bern->GetBernsteinBasisSingleQP(xi);

    std::vector<double> temp(p+1, 0.0);
    for (int jj = 0; jj < p+1; ++jj)
        for (int kk = 0; kk < p+1; ++kk)
            temp[jj] += extraction[jj*(p+1)+kk] * BB[kk];
    basis.insert(basis.end(), temp.begin(), temp.end());

    return basis;
}

std::vector<double> RefElement::GenerateBasisDerivative1DSingleQP(const BernsteinBasis * const &bern,
    const std::vector<double> &extraction,
    const double &xi, const double &h)
{
    const int p = bern->GetDegree();
    std::vector<double> dbasis{};
    const std::vector<double> dBB = bern->GetBernsteinBasisDerivativeSingleQP(xi, h);

    std::vector<double> temp(p+1, 0.0);
    for (int jj = 0; jj < p+1; ++jj)
        for (int kk = 0; kk < p+1; ++kk)
            temp[jj] += extraction[jj*(p+1)+kk] * dBB[kk];
    dbasis.insert(dbasis.end(), temp.begin(), temp.end());

    return dbasis;
}