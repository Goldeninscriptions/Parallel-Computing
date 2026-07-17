#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssemblyMF.hpp"
#include "GlobalAssembly.hpp"

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

    ElementMFSF * elemmf = new ElementMFSF(p, q);
    const int nLocBas = elemmf->GetNumLocalBasis();
    LocalAssemblyMFSF * locassemmf = new LocalAssemblyMFSF(p, q);

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, ghostID, Dir, IEN,
        NURBSExtraction1, NURBSExtraction2);

    GlobalAssemblyMF * globalassem = new GlobalAssemblyMF(nLocBas, nlocalfunc,
        nlocalelemx, nlocalelemy, ghostID);

    globalassem->AssemLoad(locassemmf, IEN,
        ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf);

    MPI_Barrier(PETSC_COMM_WORLD);

    Vec u;
    VecDuplicate(globalassem->F, &u);
    VecSet(u, 0.0);

    globalassem->MatMulMF(locassemmf,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2,
        elemmf, globalassem->F, u);
    
    delete globalassem;
    delete locassemmf;
    delete elemmf;
    delete fm;

    PetscFinalize();
    return 0;
}
