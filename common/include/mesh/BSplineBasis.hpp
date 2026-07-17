#ifndef BSPLINEBASIS_HPP
#define BSPLINEBASIS_HPP

// Ref : Les Piegl and Wayne Tiller, The NURBS Book, 2nd Edition, Springer, 1997
// FindSpan and BasisFuns are based on Algorithm A2.1 and A2.2
// DerBasisFuns is based on Algorithm A2.3

#include <vector>

class BSplineBasis
{
    public:
        BSplineBasis(const int &in_p, const std::vector<double> &in_S)
            : p(in_p), S(in_S) {}
        
        ~BSplineBasis() {}

        int GetDegree() const { return p; }

        std::vector<double> GetKnotVector() const { return S; }

        int FindSpan(const double &u) const;

        std::vector<double> BasisFuns(const double &u, const int &i) const;

        std::vector<double> DerBasisFuns(const double &u, const int &i, const int &n) const;

    private:
        const int p;
        const std::vector<double> S;
};

#endif