#ifndef ABSCISSAE_GENERATOR_HPP
#define ABSCISSAE_GENERATOR_HPP

#include <vector>

class AbscissaeGenerator
{
    public:
        std::vector<double> GenerateAbscissae1D(
            const std::vector<double> &S, const int &p, const int &idx);

        std::vector<double> GenerateAbscissae2D(const std::vector<double> &S,
            const std::vector<double> &T, const int &p, const int &q);
        
        std::vector<double> GrevilleAbscissae(
            const std::vector<double> &S, const int &p);

        // std::vector<double> DemkoAbscissae(
        //     const std::vector<double> &S, const int &p);
};
#endif