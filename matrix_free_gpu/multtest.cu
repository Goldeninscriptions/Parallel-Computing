#include <petscmat.h>
#include <cuda_runtime.h>
#include "FileManager.hpp"
#include "GlobalAssemblyMF.cuh"
#include "IENGenerator.hpp"
#include "IDGenerator.hpp"

int main (int argc, char *argv[])
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

    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> ghostID;
    std::vector<int> Dir;
    std::vector<int> IEN;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;
    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;

    std::string filename = fm->GetPartitionFilename(base_name, 0);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, ghostID, Dir, IEN, NURBSExtraction1, NURBSExtraction2);

    ElementMF * elemmf = new ElementMF(p, q);
    int nLocBas = elemmf->GetNumLocalBasis();
    GlobalAssemblyMF * globalAssembly = new GlobalAssemblyMF(
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);

    QuadraturePoint * quad1 = new QuadraturePoint(p+1, 0, 1);
    QuadraturePoint * quad2 = new QuadraturePoint(q+1, 0, 1);

    BernsteinBasis * bernstein = new BernsteinBasis(p);

    PetscLogDouble setup_start = 0.0;
    PetscLogDouble setup_end = 0.0;
    PetscTime(&setup_start);
    globalAssembly->AssemLoad(quad1, quad2,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein);
    cudaDeviceSynchronize();
    PetscTime(&setup_end);
    
    Vec x;
    Vec y;
    VecDuplicate(globalAssembly->F, &x);
    VecDuplicate(globalAssembly->F, &y);
    VecCopy(globalAssembly->F, x);
    VecSet(y, 0.0);

    // Warm up once before repeated timing.
    globalAssembly->MatMulMF(quad1, quad2,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein,
        x, y);
    cudaDeviceSynchronize();

    PetscLogDouble mult_start = 0.0;
    PetscLogDouble mult_end = 0.0;
    PetscTime(&mult_start);
    for (PetscInt iter = 0; iter < repeat; ++iter)
    {
        globalAssembly->MatMulMF(quad1, quad2,
            IEN, ID, Dir, CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elemmf, bernstein,
            x, y);
    }
    cudaDeviceSynchronize();
    PetscTime(&mult_end);

    PetscReal ynorm = 0.0;
    VecNorm(y, NORM_2, &ynorm);
    PetscPrintf(PETSC_COMM_WORLD, "Setup time: %.6f\n", setup_end - setup_start);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult total time: %.6f\n", mult_end - mult_start);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult average time: %.6f\n", (mult_end - mult_start) / repeat);
    PetscPrintf(PETSC_COMM_WORLD, "Output vector 2-norm: %.12e\n", ynorm);
    
    delete globalAssembly;
    delete elemmf;
    delete quad1;
    delete quad2;

    delete fm;
    VecDestroy(&x);
    VecDestroy(&y);
    PetscFinalize();
    return 0;
}
