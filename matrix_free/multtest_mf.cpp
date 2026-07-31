#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssemblyMF.hpp"

int main(int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscInt repeat = 20;
    PetscBool repeat_set = PETSC_FALSE;
    PetscOptionsGetInt(NULL, NULL, "-repeat", &repeat, &repeat_set);

    PetscMPIInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if (rank == 0)
    {
        std::cout << "p: " << p << std::endl;
        std::cout << "q: " << q << std::endl;
        std::cout << "Lx: " << Lx << std::endl;
        std::cout << "Ly: " << Ly << std::endl;
        std::cout << "nElemX: " << nElemX << std::endl;
        std::cout << "nElemY: " << nElemY << std::endl;
        std::cout << "part_num_1d: " << part_num_1d << std::endl;
        std::cout << "dim: " << dim << std::endl;
        std::cout << "base_name: " << base_name << std::endl;
        std::cout << "repeat: " << repeat << std::endl;
    }

    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;
    std::vector<int> ghostID;
    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> IEN;
    std::vector<int> Dir;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, ghostID, Dir, IEN,
        NURBSExtraction1, NURBSExtraction2);

    PetscLogDouble setup_start = 0.0;
    PetscLogDouble setup_end = 0.0;
    PetscTime(&setup_start);

    ElementMF * elem = new ElementMF(p, q);
    LocalAssemblyMF * locassem = new LocalAssemblyMF(p, q);
    GlobalAssemblyMF * globalassem = new GlobalAssemblyMF(elem->GetNumLocalBasis(),
        nlocalfunc, nlocalelemx, nlocalelemy, ghostID);

    globalassem->AssemLoad(locassem, IEN,
        ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elem);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscTime(&setup_end);

    Vec x;
    Vec y;
    VecDuplicate(globalassem->F, &x);
    VecDuplicate(globalassem->F, &y);
    VecCopy(globalassem->F, x);
    VecSet(y, 0.0);

    globalassem->MatMulMF(locassem,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2,
        elem, x, y);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscLogDouble mult_start = 0.0;
    PetscLogDouble mult_end = 0.0;
    PetscTime(&mult_start);
    for (PetscInt iter = 0; iter < repeat; ++iter)
    {
        globalassem->MatMulMF(locassem,
            IEN, ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2,
            elem, x, y);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscTime(&mult_end);

    PetscReal ynorm = 0.0;
    VecNorm(y, NORM_2, &ynorm);
    PetscPrintf(PETSC_COMM_WORLD, "Setup time: %.6f\n", setup_end - setup_start);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult total time: %.6f\n", mult_end - mult_start);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult average time: %.6f\n", (mult_end - mult_start) / repeat);
    PetscPrintf(PETSC_COMM_WORLD, "Output vector 2-norm: %.12e\n", ynorm);

    VecDestroy(&x);
    VecDestroy(&y);
    delete globalassem;
    delete locassem;
    delete elem;
    delete fm;

    PetscFinalize();
    return 0;
}
