#include "IDGenerator.hpp"

std::vector<int> IDGenerator::GenerateID1D(const int &nFuncX)
{
    std::vector<int> ID(nFuncX, -1);

    for (int i = 1; i < nFuncX-1; ++i) ID[i] = i - 1;

    return ID;
}

std::vector<int> IDGenerator::GenerateID1D(const BSplineBasis * const &basis)
{
    const int p = basis->GetDegree();
    const std::vector<double> S = basis->GetKnotVector();
    const int nFunc = static_cast<int>(S.size()) - p - 1;

    std::vector<int> ID(nFunc, -1);

    for (int i = 1; i < nFunc-1; ++i) ID[i] = i - 1;

    return ID;
}

std::vector<int> IDGenerator::GenerateID2D(const int &nFuncX, const int &nFuncY)
{
    std::vector<int> ID(nFuncX*nFuncY, -1);

    int A = 0;
    for (int j = 1; j < nFuncY-1; ++j)
    {
        for (int i = 1; i < nFuncX-1; ++i)
        {
            ID[j*nFuncX+i] = A;
            A++;
        }
    }

    return ID;
}

std::vector<int> IDGenerator::GenerateID2D(const BSplineBasis * const &basis1,
    const BSplineBasis * const &basis2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;

    std::vector<int> ID(m*n, -1);

    int A = 0;
    for (int j = 1; j < n-1; ++j)
    {
        for (int i = 1; i < m-1; ++i)
        {
            ID[j*m+i] = A;
            A++;
        }
    }

    return ID;
}
