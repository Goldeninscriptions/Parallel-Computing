#include "AbscissaeGenerator.hpp"

std::vector<double> AbscissaeGenerator::GenerateAbscissae1D(
    const std::vector<double> &S, const int &p, const int &idx)
{
    return GrevilleAbscissae(S, p);
}

std::vector<double> AbscissaeGenerator::GenerateAbscissae2D(
    const std::vector<double> &S, const std::vector<double> &T,
    const int &p, const int &q)
{
    std::vector<double> abscissae{};
    std::vector<double> abscissae1 = GrevilleAbscissae(S, p);
    std::vector<double> abscissae2 = GrevilleAbscissae(T, q);
    for (auto &t : abscissae2)
    {
        for (auto &s : abscissae1)
        {
            abscissae.push_back(s);
            abscissae.push_back(t);
        }
    }
    return abscissae;
}

std::vector<double> AbscissaeGenerator::GrevilleAbscissae(
    const std::vector<double> &S, const int &p)
{
    const int n = S.size() - p - 1;
    std::vector<double> abscissae(n, 0.0);
    for (int ii = 0; ii < n; ++ii)
    {
        for (int jj = 1; jj <= p; ++jj)
        {
            abscissae[ii] += S[ii + jj];
        }
        abscissae[ii] /= p;
    }
    return abscissae;
}

// std::vector<double> AbscissaeGenerator::DemkoAbscissae(
//     const std::vector<double> &S, const int &p)
// {
//     const int n = p + 1;
//     return chbpnt(S, n);
// }