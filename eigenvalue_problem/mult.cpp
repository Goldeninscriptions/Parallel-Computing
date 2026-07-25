#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssembly.hpp"

int main(int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_x, part_num_y, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_x, part_num_y, dim, base_name);

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
        std::cout << "part_num_x: " << part_num_x << std::endl;
        std::cout << "part_num_y: " << part_num_y << std::endl;
        std::cout << "dim: " << dim << std::endl;
        std::cout << "base_name: " << base_name << std::endl;
    }

    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> IEN;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;
    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, IEN, NURBSExtraction1, NURBSExtraction2);
    
    Element * elem = new Element(p, q);
    const int nLocBas = elem->GetNumLocalBasis();
    LocalAssembly * locassem = new LocalAssembly(p, q);
    GlobalAssembly * globalassem = new GlobalAssembly(IEN, ID, locassem,
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);
    
    MatSetOption(globalassem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    globalassem->AssemStiffnessLoad(locassem, IEN, ID, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elem);

    MPI_Barrier(PETSC_COMM_WORLD);

    Vec u;
    VecDuplicate(globalassem->F, &u);

    MatMult(globalassem->K, globalassem->F, u);

    delete globalassem;
    delete locassem;
    delete elem;
    delete fm;

    PetscFinalize();
    return 0;
}
