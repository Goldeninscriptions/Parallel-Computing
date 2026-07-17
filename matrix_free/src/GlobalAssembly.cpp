#include "GlobalAssembly.hpp"
#include "petsc_compat.hpp"

GlobalAssembly::GlobalAssembly(const std::vector<int> &IEN, const std::vector<int> &ID,
    const std::vector<int> &Dir, LocalAssembly * const &locassem, const int &nLocBas,
    const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{
    const int dnz = 2*nLocBas+1;
    const int onz = dnz;
    MatCreateAIJ(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, PETSC_DETERMINE,
        PETSC_DETERMINE, dnz, PC_PETSC_NULLPTR, dnz, PC_PETSC_NULLPTR, &K);
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DECIDE);
    VecSetFromOptions(F);
    VecSet(F, 0.0);
    VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    AssemNonZeroEstimate(locassem, IEN, ID, Dir);

    std::vector<int> Dnz, Onz;
    NonZeroCount(K, Dnz, Onz);

    MatDestroy(&K);

    MatCreateAIJ(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, PETSC_DETERMINE,
        PETSC_DETERMINE, 0, &Dnz[0], 0, &Onz[0], &K);
    MatSetOption(K, MAT_SYMMETRIC, PETSC_TRUE);
}

GlobalAssembly::~GlobalAssembly()
{
    MatDestroy(&K);
    VecDestroy(&F);
}

void GlobalAssembly::NonZeroCount(const Mat &K, std::vector<int> &dnz, std::vector<int> &onz)
{
    PetscInt rstart, rend;
    PetscInt cstart, cend;
    MatGetOwnershipRange(K, &rstart, &rend);
    MatGetOwnershipRangeColumn(K, &cstart, &cend);

    const PetscInt nlocalfunc = rend - rstart;
    dnz.resize(nlocalfunc, 0);
    onz.resize(nlocalfunc, 0);

    for (PetscInt i = 0; i < nlocalfunc; ++i)
    {
        const PetscInt globalrow = rstart + i;

        PetscInt ncols;
        const PetscInt *cols;
        MatGetRow(K, globalrow, &ncols, &cols, PC_PETSC_NULLPTR);

        for (PetscInt j = 0; j < ncols; ++j)
        {
            const PetscInt globalcol = cols[j];
            if (globalcol >= cstart && globalcol < cend)
            {
                ++dnz[i];
            }
            else
            {
                ++onz[i];
            }
        }

        MatRestoreRow(K, globalrow, &ncols, &cols, nullptr);
    }
}

void GlobalAssembly::AssemNonZeroEstimate(LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir)
{
    PetscInt * eID = new PetscInt[nLocBas];

    const int nlocalelem = nlocalelemx * nlocalelemy;

    for (int i = 0; i < nlocalelem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            eID[j] = ID[IEN[i*nLocBas+j]];
        }

        locassem->AssemNonZero();

        MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);
    }

    delete[] eID; eID = nullptr;

    DirichletBCK(Dir);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

void GlobalAssembly::AssemStiffnessLoad(LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    Element * const &elem)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elem->GetNumLocalBasis1D(0);
    const int qq = elem->GetNumLocalBasis1D(1);
    std::vector<double> eNURBSExtraction1(pp*pp, 0.0);
    std::vector<double> eNURBSExtraction2(qq*qq, 0.0);

    for (int jj = 0; jj < nlocalelemy; ++jj)
    {
        for (int ii = 0; ii < nlocalelemx; ++ii)
        {
            int elemIndex = jj*nlocalelemx + ii;
            for (int j = 0; j < nLocBas; ++j)
            {
                eID[j] = ID[IEN[elemIndex*nLocBas+j]];
                eCP[2*j] = CP[2*IEN[elemIndex*nLocBas+j]];
                eCP[2*j+1] = CP[2*IEN[elemIndex*nLocBas+j]+1];
            }

            std::copy(NURBSExtraction1.begin() + ii * pp * pp, 
                NURBSExtraction1.begin() + (ii + 1) * pp * pp, 
                eNURBSExtraction1.begin());
            std::copy(NURBSExtraction2.begin() + jj * qq * qq,
                NURBSExtraction2.begin() + (jj + 1) * qq * qq,
                eNURBSExtraction2.begin());
        
            elem->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            locassem->AssemLocalStiffnessLoad(elem, eCP);

            MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void GlobalAssembly::AssemStiffnessLoad(LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    ElementFEM * const &elem)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);

    for (int jj = 0; jj < nlocalelemy; ++jj)
    {
        for (int ii = 0; ii < nlocalelemx; ++ii)
        {
            int elemIndex = jj*nlocalelemx + ii;
            for (int j = 0; j < nLocBas; ++j)
            {
                eID[j] = ID[IEN[elemIndex*nLocBas+j]];
                eCP[2*j] = CP[2*IEN[elemIndex*nLocBas+j]];
                eCP[2*j+1] = CP[2*IEN[elemIndex*nLocBas+j]+1];
            }

            locassem->AssemLocalStiffnessLoad(elem, eCP);

            MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void GlobalAssembly::DirichletBCK(const std::vector<int> &Dir)
{
    for (int i = 0; i < static_cast<int>(Dir.size()); ++i)
        MatSetValue(K, Dir[i], Dir[i], 1.0, ADD_VALUES);
}

void GlobalAssembly::DirichletBC(const std::vector<int> &Dir)
{
    for (int i = 0; i < static_cast<int>(Dir.size()); ++i)
    {
        MatSetValue(K, Dir[i], Dir[i], 1.0, ADD_VALUES);
        VecSetValue(F, Dir[i], 0.0, INSERT_VALUES);
    }
}
