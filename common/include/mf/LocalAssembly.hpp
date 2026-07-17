#ifndef LOCALASSEMBLY_HPP
#define LOCALASSEMBLY_HPP

#include <petscmat.h>
#include "Element.hpp"
#include "ElementFEM.hpp"
#include "QuadraturePoint.hpp"

class LocalAssembly
{
    public:
        PetscScalar *Kloc;
        PetscScalar *Floc;

        LocalAssembly(const int &p, const int &q)
            : n((p+1)*(q+1))
        {
            quad1 = new QuadraturePoint(p+1, 0, 1);
            quad2 = new QuadraturePoint(q+1, 0, 1);
            Kloc = new PetscScalar[n*n];
            Floc = new PetscScalar[n];
        }

        ~LocalAssembly()
        {
            delete quad1; quad1 = nullptr;
            delete quad2; quad2 = nullptr;
            delete[] Kloc; Kloc = nullptr;
            delete[] Floc; Floc = nullptr;
        }

        void AssemLocalStiffnessLoad(const Element * const &elem,
            const std::vector<double> &eCP);
        
        void AssemLocalStiffnessLoad(const ElementFEM * const &elem,
            const std::vector<double> &eCP);

        void AssemNonZero()
        {
            for (int i = 0; i < n*n; ++i)
                Kloc[i] = 1.0;
        }

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
            for (int i = 0; i < n*n; ++i)
                Kloc[i] = 0.0;
            for (int i = 0; i < n; ++i)
                Floc[i] = 0.0;
        }
};

#endif