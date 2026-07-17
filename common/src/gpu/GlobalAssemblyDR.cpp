#include "GlobalAssemblyDR.hpp"
#include "petsc_compat.hpp"

GlobalAssemblyDR::GlobalAssemblyDR(FileManager * const &fm, const std::vector<int> &IEN,
    const std::vector<int> &ID, const std::vector<int> &Dir,
    LocalAssembly * const &locassem, const int &nLocBas,
    const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{
    const int dnz = 2*nLocBas+1;
    const int onz = dnz;
    MatCreateAIJ(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, PETSC_DETERMINE,
        PETSC_DETERMINE, dnz, PC_PETSC_NULLPTR, dnz, PC_PETSC_NULLPTR, &K);
    MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    AssemNonZeroEstimate(locassem, IEN, ID, Dir);

    int nnz;
    std::vector<int> rows, cols;
    NonZeroCoordinate(K, nnz, rows, cols);

    PetscInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::string filename = fm->GetNonZeroCoordinateFilename("coordinate", rank);
    fm->WriteNonZeroCoordinate(filename, nnz, rows, cols);
}

GlobalAssemblyDR::~GlobalAssemblyDR()
{
    MatDestroy(&K);
}

void GlobalAssemblyDR::AssemNonZeroEstimate(LocalAssembly * const &locassem,
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

void GlobalAssemblyDR::NonZeroCoordinate(const Mat &K, int &nnz, std::vector<int> &rows, std::vector<int> &cols)
{
    PetscInt rstart, rend;
    MatGetOwnershipRange(K, &rstart, &rend);

    PetscInt nrows;
    const PetscInt *ia, *ja;
    PetscBool done;

    MatGetRowIJ(K, 0, PETSC_FALSE, PETSC_FALSE, &nrows, &ia, &ja, &done);

    nnz = ia[nrows] - ia[0];
    rows.clear();
    cols.clear();

    for (PetscInt i = 0; i < nrows; ++i)
    {
        for (PetscInt j = ia[i]; j < ia[i+1]; ++j)
        {
            rows.push_back(rstart + i);
            cols.push_back(ja[j]);
        }
    }
}

void GlobalAssemblyDR::DirichletBCK(const std::vector<int> &Dir)
{
    for (int ii = 0; ii < Dir.size(); ++ii)
    {
        int row = Dir[ii];
        MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
}
