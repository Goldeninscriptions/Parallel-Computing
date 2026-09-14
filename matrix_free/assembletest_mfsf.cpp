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

    FileManager *fm = new FileManager();
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

    ElementMFSF *elem = new ElementMFSF(p, q);
    LocalAssemblyMFSF *locassem = new LocalAssemblyMFSF(p, q);

    PetscLogDouble assemble_start = 0.0;
    PetscLogDouble assemble_end = 0.0;
    PetscTime(&assemble_start);

    PetscReal fnorm = 0.0;
    for (PetscInt iter = 0; iter < repeat; ++iter)
    {
        GlobalAssemblyMF *globalassem = new GlobalAssemblyMF(elem->GetNumLocalBasis(),
            nlocalfunc, nlocalelemx, nlocalelemy, ghostID);

        globalassem->AssemLoad(locassem, IEN,
            ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elem);

        if (iter == repeat - 1)
        {
            VecNorm(globalassem->F, NORM_2, &fnorm);
        }

        delete globalassem;
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscTime(&assemble_end);

    PetscPrintf(PETSC_COMM_WORLD, "Assemble total time: %.6f\n", assemble_end - assemble_start);
    PetscPrintf(PETSC_COMM_WORLD, "Assemble average time: %.6f\n", (assemble_end - assemble_start) / repeat);
    PetscPrintf(PETSC_COMM_WORLD, "Assembled load vector 2-norm: %.12e\n", fnorm);

    delete locassem;
    delete elem;
    delete fm;

    PetscFinalize();
    return 0;
}
