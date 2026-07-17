#ifndef GLOBALASSEMBLYDR_HPP
#define GLOBALASSEMBLYDR_HPP

#include <petscmat.h>
#include "LocalAssembly.hpp"
#include "FileManager.hpp"

class GlobalAssemblyDR
{
    public:
        Mat K;

        GlobalAssemblyDR(FileManager * const &fm, const std::vector<int> &IEN,
            const std::vector<int> &ID, const std::vector<int> &Dir,
            LocalAssembly * const &locassem, const int &nLocBas,
            const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy);

        ~GlobalAssemblyDR();

        void NonZeroCoordinate(const Mat &K, int &nnz, std::vector<int> &rows, std::vector<int> &cols);

        void AssemNonZeroEstimate(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir);
    
    private:
        const int nLocBas;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;

        void DirichletBCK(const std::vector<int> &Dir);
};

#endif