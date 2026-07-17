#include "BSplineBasis.hpp"

int BSplineBasis::FindSpan(const double &u) const
{
    const int nFunc = static_cast<int>(S.size()) - p - 1;
    if (u == S[nFunc + 1]) return nFunc;
    int low = p;
    int high = nFunc + 1;
    int mid = (low + high) / 2;
    while (u < S[mid] || u >= S[mid + 1])
    {
        if (u < S[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

std::vector<double> BSplineBasis::BasisFuns(const double &u, const int &i) const
{
    std::vector<double> N(p + 1, 0.0);
    std::vector<double> left(p + 1, 0.0);
    std::vector<double> right(p + 1, 0.0);
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j)
    {
        left[j] = u - S[i + 1 - j];
        right[j] = S[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r)
        {
            double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    return N;
}

std::vector<double> BSplineBasis::DerBasisFuns(const double &u, const int &i, const int &n) const
{
    std::vector<double> ndu((p+1)*(p+1), 0.0);
    std::vector<double> ders((n+1)*(p+1), 0.0);
    std::vector<double> left(p + 1, 0.0);
    std::vector<double> right(p + 1, 0.0);
    ndu[0] = 1.0;
    for (int j = 1; j <= p; j++)
    {
        left[j] = u - S[i + 1 - j];
        right[j] = S[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; r++)
        {
            ndu[j*(p+1)+r] = right[r + 1] + left[j - r];
            double temp = ndu[r*(p+1)+j-1] / ndu[j*(p+1)+r];
            ndu[r*(p+1)+j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j*(p+1)+j] = saved;
    }
    for (int j = 0; j <= p; j++)
    {
        ders[j] = ndu[j*(p+1)+p];
    }
    for (int r = 0; r <= p; r++)
    {
        int s1 = 0;
        int s2 = 1;
        std::vector<double> a(2*(p+1), 0.0);
        a[0] = 1.0;
        for (int k = 1; k <= n; k++)
        {
            double d = 0.0;
            int rk = r - k;
            int pk = p - k;
            if (r >= k)
            {
                a[s2*(p+1)+0] = a[s1*(p+1)+0] / ndu[(pk+1)*(p+1)+rk];
                d = a[s2*(p+1)+0] * ndu[rk*(p+1)+pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; j++)
            {
                a[s2*(p+1)+j] = (a[s1*(p+1)+j] - a[s1*(p+1)+j-1]) / ndu[(pk+1)*(p+1)+rk+j];
                d += a[s2*(p+1)+j] * ndu[(rk+j)*(p+1)+pk];
            }
            if (r <= pk)
            {
                a[s2*(p+1)+k] = -a[s1*(p+1)+k-1] / ndu[(pk+1)*(p+1)+r];
                d += a[s2*(p+1)+k] * ndu[r*(p+1)+pk];
            }
            ders[(k*(p+1))+r] = d;
            int j = s1;
            s1 = s2;
            s2 = j;
        }
    }
    int r = p;
    for (int k = 1; k <= n; k++)
    {
        for (int j = 0; j <= p; j++)
        {
            ders[k*(p+1)+j] *= r;
        }
        r *= (p - k);
    }
    std::vector<double> dersOut(ders.begin()+(n*(p+1)), ders.end());
    return dersOut;
}