#ifndef ELEMENTMFSF_HPP
#define ELEMENTMFSF_HPP

#include <algorithm>
#include "RefElement.hpp"
#include "QuadraturePoint.hpp"

class ElementMFSF
{
    public:
        ElementMFSF(const int &p, const int &q) : p(p), q(q), nLocBas((p+1)*(q+1)) {}
        ~ElementMFSF(){}

        void SetElement(const std::vector<double> &extraction1,
            const std::vector<double> &extraction2,
            const double &hx, const double &hy)
        {
            this->extraction1 = extraction1;
            this->extraction2 = extraction2;
            this->hx = hx;
            this->hy = hy;
        }

        int GetNumLocalBasis() const { return nLocBas; }

        int GetNumLocalBasis1D(const int &dim) const
        {
            return (dim == 0) ? p+1 : q+1;
        }

        void GenerateBSplineBasis1D(const double &xi,
            std::vector<double> &N1, std::vector<double> &N2) const;

        void GenerateBSplineBasis1D(const double &xi,
            std::vector<double> &N1, std::vector<double> &N2,
            std::vector<double> &dN1, std::vector<double> &dN2) const;
        
        void GenerateElementSingleQP(const std::vector<double> &eCP,
            const std::vector<double> &N1, const std::vector<double> &N2,
            const std::vector<double> &dN1, const std::vector<double> &dN2,
            std::vector<double> &R,
            double &jacobian) const;

        void GenerateElementSingleQP(const std::vector<double> &eCP, 
            const std::vector<double> &N1, const std::vector<double> &N2,
            const std::vector<double> &dN1, const std::vector<double> &dN2,
            double &w, double &jacobian,
            double &dw_dx, double &dw_dy) const;
        
        void GenerateElement(const QuadraturePoint * const &quad1,
            const QuadraturePoint * const &quad2,
            const std::vector<double> &eCP,
            std::vector<double> &R,
            std::vector<double> &J);
        
        void GenerateElement(const QuadraturePoint * const &quad1,
            const QuadraturePoint * const &quad2,
            const std::vector<double> &eCP,
            std::vector<double> &B1, std::vector<double> &B2,
            std::vector<double> &dB1, std::vector<double> &dB2,
            std::vector<double> &W, std::vector<double> &J,
            std::vector<double> &dW_dx, std::vector<double> &dW_dy) const;
    
    private:
        const int p;
        const int q;
        const int nLocBas;

        std::vector<double> extraction1{};
        std::vector<double> extraction2{};
        double hx;
        double hy;
};

#endif