#ifndef ELEMENTMF_HPP
#define ELEMENTMF_HPP

#include <algorithm>
#include "RefElement.hpp"
#include "QuadraturePoint.hpp"

class ElementMF
{
    public:
        ElementMF(const int &p, const int &q) : p(p), q(q), nLocBas((p+1)*(q+1)) {}
        ~ElementMF(){}

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

        void GenerateElementSingleQP(const double &xi, const double &eta,
            const std::vector<double> &eCP, 
            std::vector<double> &R, std::vector<double> &dR_dx, std::vector<double> &dR_dy,
            double &jacobian) const;
        
        void GenerateElement(const QuadraturePoint * const &quad1,
            const QuadraturePoint * const &quad2,
            const std::vector<double> &eCP);
        
        double get_R(const int &qpIdx, const int &basIdx)
        {
            return RR[qpIdx*nLocBas+basIdx];
        }

        double get_dR_dx(const int &qpIdx, const int &basIdx)
        {
            return dRR_dx[qpIdx*nLocBas+basIdx];
        }

        double get_dR_dy(const int &qpIdx, const int &basIdx)
        {
            return dRR_dy[qpIdx*nLocBas+basIdx];
        }

        double get_JxW(const int &qpIdx)
        {
            return JxW[qpIdx];
        }

        int get_p() const { return p; }
        int get_q() const { return q; }
    
    private:
        const int p;
        const int q;
        const int nLocBas;

        std::vector<double> RR{};
        std::vector<double> dRR_dx{};
        std::vector<double> dRR_dy{};

        std::vector<double> JxW{};

        std::vector<double> extraction1{};
        std::vector<double> extraction2{};
        double hx;
        double hy;
};

#endif