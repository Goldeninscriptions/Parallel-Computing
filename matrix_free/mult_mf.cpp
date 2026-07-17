#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssemblyMF.hpp"
#include "GlobalAssembly.hpp"
#include "MatrixFreeMode.hpp"

int main(int argc, char *argv[])
{
    MatrixFreeMode mode = MatrixFreeMode::MF;
    if (!ParseMatrixFreeMode(argc, argv, mode))
        return 0;

    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    PetscInitialize(&argc, &argv, NULL, NULL);

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
        std::cout << "matrix_free_mode: " << MatrixFreeModeName(mode) << std::endl;
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

    if (mode == MatrixFreeMode::MF)
    {
        ElementMF * elem = new ElementMF(p, q);
        LocalAssemblyMF * locassem = new LocalAssemblyMF(p, q);
        GlobalAssemblyMF * globalassem = new GlobalAssemblyMF(elem->GetNumLocalBasis(),
            nlocalfunc, nlocalelemx, nlocalelemy, ghostID);

        globalassem->AssemLoad(locassem, IEN,
            ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elem);

        MPI_Barrier(PETSC_COMM_WORLD);

        Vec u;
        VecDuplicate(globalassem->F, &u);
        VecSet(u, 0.0);

        globalassem->MatMulMF(locassem,
            IEN, ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2,
            elem, globalassem->F, u);

        VecDestroy(&u);
        delete globalassem;
        delete locassem;
        delete elem;
    }
    else
    {
        ElementMFSF * elem = new ElementMFSF(p, q);
        LocalAssemblyMFSF * locassem = new LocalAssemblyMFSF(p, q);
        GlobalAssemblyMF * globalassem = new GlobalAssemblyMF(elem->GetNumLocalBasis(),
            nlocalfunc, nlocalelemx, nlocalelemy, ghostID);

        globalassem->AssemLoad(locassem, IEN,
            ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elem);

        MPI_Barrier(PETSC_COMM_WORLD);

        Vec u;
        VecDuplicate(globalassem->F, &u);
        VecSet(u, 0.0);

        globalassem->MatMulMF(locassem,
            IEN, ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2,
            elem, globalassem->F, u);

        VecDestroy(&u);
        delete globalassem;
        delete locassem;
        delete elem;
    }

    delete fm;

    PetscFinalize();
    return 0;
}
