#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <algorithm>
#include "RefElement.hpp"

class Element
{
    public:
        Element(const int &p, const int &q) : p(p), q(q) {}
        ~Element(){}

        void SetElement(const std::vector<double> &extraction1,
            const std::vector<double> &extraction2,
            const double &hx, const double &hy)
        {
            this->extraction1 = extraction1;
            this->extraction2 = extraction2;
            this->hx = hx;
            this->hy = hy;
        }

        int GetNumLocalBasis() const { return (p+1)*(q+1); }

        int GetNumLocalBasis1D(const int &dim) const
        {
            return (dim == 0) ? p+1 : q+1;
        }

        void GenerateElementSingleQP(const double &xi, const double &eta,
            const std::vector<double> &eCP, 
            std::vector<double> &R, std::vector<double> &dR_dx, std::vector<double> &dR_dy,
            double &x, double &y,
            double &jacobian) const;
    
    private:
        const int p;
        const int q;

        std::vector<double> extraction1{};
        std::vector<double> extraction2{};
        double hx;
        double hy;
};

#endif