#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

#include <vector>
#include <iostream>

class BernsteinBasis {
    public:
        BernsteinBasis(const int &p) : p(p) {}
        ~BernsteinBasis(){}

        std::vector<double> GetBernsteinBasisSingleQP(const double &xi) const;

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const double &xi) const;

        int GetDegree() const { return p; }

    private:
        const int p;
};

#endif