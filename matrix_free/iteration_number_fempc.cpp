#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssembly.hpp"

typedef struct {
    KSP innerksp;
} MyPCCtx;

PetscErrorCode MyPCApply(PC pc, Vec x, Vec y)
{
    MyPCCtx *ctx;
    PCShellGetContext(pc, (void**)&ctx);
    KSPSolve(ctx->innerksp, x, y);
    return 0;
}

PetscErrorCode MyPCDestroy(PC pc)
{
    MyPCCtx *ctx;
    PCShellGetContext(pc, (void**)&ctx);
    KSPDestroy(&ctx->innerksp);
    PetscCall(PetscFree(ctx));
    return 0;
}

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

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, ghostID, Dir, IEN, NURBSExtraction1, NURBSExtraction2);
    
    Element * elem = new Element(p, q);
    const int nLocBas = elem->GetNumLocalBasis();
    LocalAssembly * locassem = new LocalAssembly(p, q);
    GlobalAssembly * globalassem = new GlobalAssembly(IEN, ID, Dir, locassem,
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);
    
    MatSetOption(globalassem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    globalassem->AssemStiffnessLoad(locassem, IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elem);

    MPI_Barrier(PETSC_COMM_WORLD);

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
    GlobalAssembly * globalassem_fem = new GlobalAssembly(IEN_fem, ID_fem, Dir_fem, locassem_fem,
        4, nlocalfunc_fem, nlocalelemx_fem, nlocalelemy_fem);
    
    MatSetOption(globalassem_fem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    
    globalassem_fem->AssemStiffnessLoad(locassem_fem, IEN_fem, ID_fem, Dir_fem, CP_fem, elem_fem);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling stiffness matrix and load vector...done\n");

    // Vec u;
    // VecDuplicate(globalassem->F, &u);
    // VecSet(u, 0.0);
    // MatMult(globalassem->K, globalassem->F, u);

    // Vec u1;
    // VecDuplicate(u, &u1);
    // VecSet(u1, 0.0);
    // MatMult(globalassem->K, u, u1);

    // VecView(globalassem->F, PETSC_VIEWER_STDOUT_WORLD);
    // VecView(u, PETSC_VIEWER_STDOUT_WORLD);
    // VecView(u1, PETSC_VIEWER_STDOUT_WORLD);

    PetscLogDouble tstart, tend;
    PetscTime(&tstart);

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, globalassem->K, globalassem->K);
    // KSPSetComputeSingularValues(ksp, PETSC_TRUE);
    KSPSetFromOptions(ksp);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);

    MyPCCtx *ctx;
    PetscNew(&ctx);
    KSPCreate(PETSC_COMM_WORLD, &ctx->innerksp);
    KSPSetOperators(ctx->innerksp, globalassem_fem->K, globalassem_fem->K);
    KSPSetType(ctx->innerksp, KSPCG);

    PC innerpc;
    KSPGetPC(ctx->innerksp, &innerpc);
    PCSetType(innerpc, PCJACOBI);
    PCSetFromOptions(innerpc);

    PetscReal rtol = 1e-10;
    PetscReal abstol = 1e-10;
    PetscReal divtol = 1e4;
    PetscInt maxits = 10000;
    KSPSetTolerances(ctx->innerksp, rtol, abstol, divtol, maxits);
    KSPSetTolerances(ksp, rtol, abstol, divtol, maxits);

    PCShellSetContext(pc, ctx);
    PCShellSetApply(pc, MyPCApply);
    PCShellSetDestroy(pc, MyPCDestroy);

    Vec u;
    VecDuplicate(globalassem->F, &u);
    VecSet(u, 0.0);
    KSPSolve(ksp, globalassem->F, u);

    PetscTime(&tend);
    PetscLogDouble time = tend - tstart;
    PetscPrintf(PETSC_COMM_WORLD, "Time: %f\n", time);

    PetscInt num_iterations;
    KSPGetIterationNumber(ksp, &num_iterations);

    if (rank == 0)
    {
        std::cout << "Number of KSP iterations: " << num_iterations << std::endl;
    }

    VecView(u, PETSC_VIEWER_STDOUT_WORLD);

    delete fm; fm = nullptr;
    delete elem; elem = nullptr;
    delete locassem; locassem = nullptr;
    delete globalassem; globalassem = nullptr;
    delete elem_fem; elem_fem = nullptr;
    delete locassem_fem; locassem_fem = nullptr;
    delete globalassem_fem; globalassem_fem = nullptr;

    VecDestroy(&u);
    KSPDestroy(&ksp);
    
    PetscFinalize();
    return 0;
}
