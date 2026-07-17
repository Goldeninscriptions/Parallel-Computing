#ifndef GLOBALASSEMBLY_HPP
#define GLOBALASSEMBLY_HPP

#include <petscmat.h>
#include "QuadraturePoint.hpp"
#include "ElementFEM.hpp"
#include "memory.cuh"
#include "assembly.cuh"

class GlobalAssembly
{
    public:
        Mat K;

        GlobalAssembly(const int &nLocBas,
            const int &nnz, const int &nlocalfunc,
            const int &nlocalelemx, const int &nlocalelemy,
            const std::vector<int> &rows,
            const std::vector<int> &cols);

        ~GlobalAssembly();
        
        void AssemStiffness(QuadraturePoint * const &quad1,
            QuadraturePoint * const &quad2,
            const std::vector<int> &IEN,
            const std::vector<int> &dir2coo,
            const std::vector<double> &CP,
            const std::vector<int> &elem2coo,
            ElementFEM * const &elem);
    
    private:
        const int nLocBas;
        const int nnz;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;
};

#endif