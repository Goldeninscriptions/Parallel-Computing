#ifndef ELEM2COOGENERATOR_HPP
#define ELEM2COOGENERATOR_HPP

#include <vector>
#include <algorithm>
#include <iterator>
#include <omp.h>

class Elem2COOGenerator
{
    public:
        Elem2COOGenerator(const int &nLocBas, const int &nnz,
            const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy)
        : nLocBas(nLocBas), nnz(nnz), nlocalfunc(nlocalfunc),
            nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy) {}

        ~Elem2COOGenerator() {}

        void GenerateElem2COO(const std::vector<int> &IEN,
            const std::vector<int> &ID, const std::vector<int> &rows,
            const std::vector<int> &cols, std::vector<int> &elem2coo);

        void GenerateDir2COO(const std::vector<int> &ID,
            const std::vector<int> &rows, std::vector<int> &dir2coo);
    
    private:
        const int nLocBas;
        const int nnz;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;
};

#endif