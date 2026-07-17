#ifndef LOCALASSEMBLYMF_HPP
#define LOCALASSEMBLYMF_HPP

#include <petscmat.h>
#include "ElementMF.hpp"

class LocalAssemblyMF
{
    public:
        PetscScalar *Floc_in;
        PetscScalar *Floc_out;
        PetscScalar *Floc;

        LocalAssemblyMF(const int &p, const int &q)
            : n((p+1)*(q+1))
        {
            quad1 = new QuadraturePoint(p+1, 0, 1);
            quad2 = new QuadraturePoint(q+1, 0, 1);
            Floc = new PetscScalar[n];
            Floc_in = new PetscScalar[n];
            Floc_out = new PetscScalar[n];
        }

        ~LocalAssemblyMF()
        {
            delete quad1; quad1 = nullptr;
            delete quad2; quad2 = nullptr;
            delete[] Floc; Floc = nullptr;
            delete[] Floc_in; Floc_in = nullptr;
            delete[] Floc_out; Floc_out = nullptr;
        }

        void AssemLocalLoad(ElementMF * const &elem,
            const std::vector<double> &eCP);
        
        void LocalMatMulMF(ElementMF * const &elem,
            const std::vector<double> &eCP);

    private:
        const int n;
        QuadraturePoint * quad1;
        QuadraturePoint * quad2;

        double Getf(const double &xi, const double &eta)
        {
            return xi*(1-xi)*eta*(1-eta);
        }

        void ResetStiffnessLoad()
        {
            for (int i = 0; i < n; ++i)
                Floc[i] = 0.0;
        }

        void ResetLocalLoad()
        {
            for (int i = 0; i < n; ++i)
                Floc_out[i] = 0.0;
        }
};

#endif