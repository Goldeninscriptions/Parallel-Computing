#include "AbscissaeGenerator.hpp"

#include <algorithm>
#include <cmath>

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

std::vector<double> AbscissaeGenerator::DemkoAbscissae(
    const std::vector<double> &S, const int &p)
{
    const int n = static_cast<int>(S.size()) - p - 1;
    std::vector<double> abscissae(n, 0.0);

    if (n <= 1)
    {
        if (n == 1) abscissae[0] = S.front();
        return abscissae;
    }

    const double a = S.front();
    const double b = S.back();
    const double mid = 0.5 * (a + b);
    const double half = 0.5 * (b - a);

    // Demko points are not implemented in the original codebase.
    // We approximate them with Chebyshev-Lobatto nodes on the physical interval,
    // which preserves endpoint inclusion and boundary clustering.
    for (int ii = 0; ii < n; ++ii)
    {
        const double theta = M_PI * static_cast<double>(n - 1 - ii) / static_cast<double>(n - 1);
        abscissae[ii] = mid + half * std::cos(theta);
    }

    std::sort(abscissae.begin(), abscissae.end());
    return abscissae;
}
