#include "NURBSExtractionGenerator.hpp"
#include <iostream>

std::vector<double> NURBSExtractionGenerator::GenerateExtraction1D(const BSplineBasis * const &basis)
{
    const std::vector<double> S = basis->GetKnotVector();
    const int p = basis->GetDegree();
    const int m = static_cast<int>(S.size());
    const int nFunc = m - p - 1;
    const int nElem = nFunc - p;

    const int a = p + 1;
    std::vector<double> C(a * a * nElem, 0.0);

    std::vector<double> I(a * a, 0.0);
    for (int i = 0; i < a; ++i) I[i * a + i] = 1.0;

    for (int i = 1; i <= a; ++i) C[idxC(i, i, 1, a)] = 1.0;

    int nb = 1;
    int A = p + 1;
    int B = A + 1;

    while (B < m)
    {
        for (int i = 1; i <= a; ++i) C[idxC(i, i, nb+1, a)] = 1.0;
        int iVal = B;

        while (B < m && S[B] == S[B-1]) B++;
        int mult = B - iVal + 1;

        if (mult < p)
        {
            double numer = S[B - 1] - S[A - 1];

            int r = p - mult;
            std::vector<double> alphas(r, 0.0);

            for (int j = p; j > mult; --j)
                alphas[j - mult - 1] = numer / (S[A + j - 1] - S[A - 1]);

            for (int j = 1; j <= r; ++j)
            {
                int save = r - j + 1;
                int s = mult + j;

                for (int k = p + 1; k > s; --k)
                {
                    double alpha = alphas[k - s - 1];

                    for (int row = 1; row <= a; ++row)
                    {
                        double oldVal_k = C[idxC(row, k, nb, a)];
                        double oldVal_km1 = C[idxC(row, k - 1, nb, a)];
                        C[idxC(row, k, nb, a)] = alpha * oldVal_k + (1.0 - alpha) * oldVal_km1;
                    }
                }

                if (B < m)
                {
                    for (int offset = 0; offset <= j; ++offset)
                    {
                        int rowLeft = save + offset;
                        int rowRight = p - j + 1 + offset;

                        C[idxC(rowLeft, save, nb+1, a)] = C[idxC(rowRight, p + 1, nb, a)];
                    }
                }
            }

            nb++;
            if (nb >= nElem) break;
            if (B < m)
            {
                A = B;
                B++;
            }
        }
    }

    return C;
}