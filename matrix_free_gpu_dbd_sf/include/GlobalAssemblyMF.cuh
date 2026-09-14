#ifndef GLOBALASSEMBLYMF_HPP
#define GLOBALASSEMBLYMF_HPP

#include <petscmat.h>
#include "QuadraturePoint.hpp"
#include "ElementMF.hpp"
#include "memory.cuh"
#include "mult.cuh"
#include "InvLM.hpp"

class GlobalAssemblyMF
{
    public:
        Vec F;

        GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
            const int &nlocalelemx, const int &nlocalelemy);
            
        ~GlobalAssemblyMF()
        {
            VecDestroy(&F);
        }

        void AssemLoad(QuadraturePoint * const &quad1,
            QuadraturePoint * const &quad2,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir,
            const std::vector<int> &invlm_elemNum,
            const std::vector<int> &invlm_offset,
            const std::vector<int> &invlm_elemIdx,
            const std::vector<int> &invlm_baseIdx,
            const std::vector<int> &xelemIdx,
            const std::vector<int> &yelemIdx,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            ElementMF * const &elemmf,
            BernsteinBasis * const &bernstein);
        
        void MatMulMF(QuadraturePoint * const &quad1,
            QuadraturePoint * const &quad2,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir,
            const std::vector<int> &invlm_elemNum,
            const std::vector<int> &invlm_offset,
            const std::vector<int> &invlm_elemIdx,
            const std::vector<int> &invlm_baseIdx,
            const std::vector<int> &xelemIdx,
            const std::vector<int> &yelemIdx,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            ElementMF * const &elemmf,
            BernsteinBasis * const &bernstein,
            Vec x, Vec y);
    
    private:
        const int nLocBas;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;
};

#endif