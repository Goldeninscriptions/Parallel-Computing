#ifndef QUADRATUREPOINT_HPP
#define QUADRATUREPOINT_HPP

#include <vector>
#include <iostream>

class QuadraturePoint
{
    public:
        QuadraturePoint(const int &nqp, const int &left, const int &right);
        ~QuadraturePoint() {}

        int GetNumQuadraturePoint() const { return nqp; }
        std::vector<double> GetQuadraturePoint() const { return qp; }
        std::vector<double> GetWeight() const { return w; }

    private:
        const int nqp;
        std::vector<double> qp{};
        std::vector<double> w{};
};

#endif