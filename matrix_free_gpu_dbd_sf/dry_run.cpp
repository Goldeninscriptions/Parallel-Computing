#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssemblyDR.hpp"

int main(int argc, char *argv[])
{
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
    }

    std::vector<double> CP_fem;
    std::vector<int> ID_fem;
    std::vector<int> Dir_fem;
    std::vector<int> IEN_fem;
    int nlocalfunc_fem;
    int nlocalelemx_fem;
    int nlocalelemy_fem;

    std::string base_name_fem = "part_fem";
    std::string filename_fem = fm->GetPartitionFilename(base_name_fem, rank);
    fm->ReadPartition(filename_fem, nlocalfunc_fem,
        nlocalelemx_fem, nlocalelemy_fem,
        CP_fem, ID_fem, Dir_fem, IEN_fem);

    ElementFEM * elem_fem = new ElementFEM(1, 1);
    LocalAssembly * locassem_fem = new LocalAssembly(1, 1);
    GlobalAssemblyDR * globalassem_fem = new GlobalAssemblyDR(fm,
        IEN_fem, ID_fem, Dir_fem, locassem_fem,
        4, nlocalfunc_fem, nlocalelemx_fem, nlocalelemy_fem);
    
    MatSetOption(globalassem_fem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling stiffness matrix and load vector...done\n");

    delete fm; fm = nullptr;
    delete elem_fem; elem_fem = nullptr;
    delete locassem_fem; locassem_fem = nullptr;
    delete globalassem_fem; globalassem_fem = nullptr;
    
    PetscFinalize();
    return 0;
}