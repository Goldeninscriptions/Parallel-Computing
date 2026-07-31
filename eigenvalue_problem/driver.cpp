#include <petscksp.h>
#include <slepceps.h>
#include <slepcmath.h>
#include <slepcsvd.h>
#include "FileManager.hpp"
#include "GlobalAssembly.hpp"

double GetMaximumEigenvalue(Mat &Kiga, Mat &Kfem, Vec &u0,
    const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, Kfem, Kfem);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    while (true)
    {
        MatMult(Kiga, u0, u);

        KSPSolve(ksp, u, u);

        KSPSolve(ksp, u, u);

        MatMult(Kiga, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);
    KSPDestroy(&ksp);

    return lambda_new;
}

double GetMaximumEigenvalue(Mat &K, Vec &u0, const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    while (true)
    {
        MatMult(K, u0, u);

        MatMult(K, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);

    return lambda_new;
}

double GetMinimumEigenvalue(Mat &K, Vec &u0, const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, K, K);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    while (true)
    {
        KSPSolve(ksp, u0, u);

        KSPSolve(ksp, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);
    KSPDestroy(&ksp);

    return lambda_new;
}

typedef struct
{
    Mat K;
    KSP ksp;
} UserCtx;

typedef struct
{
    KSP ksp;
} UserCtxK;

PetscReal ComputeExtremeSingularValue(Mat A, SVDWhich which, PetscReal tol, PetscInt max_it);
PetscReal ComputeExtremeEigenvalue(Mat A, EPSWhich which, PetscReal tol, PetscInt max_it);

enum class SpectralMethod
{
    SVD,
    EPS
};

void PrintUnavailable(const char *matrix_name)
{
    PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value of %s: unavailable\n", matrix_name);
    PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of %s: unavailable\n", matrix_name);
    PetscPrintf(PETSC_COMM_WORLD, "Condition number of %s: unavailable\n", matrix_name);
}

void AnalyzeWithEPS(Mat A_neg, const char *matrix_name, PetscReal tol, PetscInt max_it)
{
    PetscReal lmax = ComputeExtremeEigenvalue(A_neg, EPS_LARGEST_REAL, tol, max_it);
    PetscPrintf(PETSC_COMM_WORLD, "Maximum eigenvalue of -%s: %.15g\n", matrix_name, lmax);

    PetscReal lmin = ComputeExtremeEigenvalue(A_neg, EPS_SMALLEST_REAL, tol, max_it);
    if (lmin == PETSC_MAX_REAL)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Minimum eigenvalue of -%s: unavailable\n", matrix_name);
        PrintUnavailable(matrix_name);
        return;
    }

    PetscPrintf(PETSC_COMM_WORLD, "Minimum eigenvalue of -%s: %.15g\n", matrix_name, lmin);
    PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value of %s: %.15g\n", matrix_name, lmax);
    PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of %s: %.15g\n", matrix_name, lmin);
    PetscPrintf(PETSC_COMM_WORLD, "Condition number of %s: %.15g\n", matrix_name, lmax / lmin);
}

PetscReal ComputeExtremeSingularValue(Mat A, SVDWhich which, PetscReal tol, PetscInt max_it)
{
    SVD svd;
    SVDCreate(PETSC_COMM_WORLD, &svd);
    SVDSetOperators(svd, A, NULL);
    SVDSetProblemType(svd, SVD_STANDARD);
    SVDSetType(svd, SVDTRLANCZOS);
    SVDTRLanczosSetOneSide(svd, PETSC_FALSE);
    SVDSetTolerances(svd, tol, max_it);
    if (which == SVD_SMALLEST)
    {
        /*
         * TRLanczos requires ncv <= nsv + mpd. Since nsv=1 here, choose a
         * larger-but-valid subspace to help the smallest singular value solve
         * without triggering setup errors.
         */
        SVDSetDimensions(svd, 1, 25, 24);
    }
    else
    {
        SVDSetDimensions(svd, 1, 17, 16);
    }

    SVDSetFromOptions(svd);
    SVDSetWhichSingularTriplets(svd, which);
    SVDSolve(svd);

    SVDConvergedReason reason;
    SVDGetConvergedReason(svd, &reason);
    PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

    PetscReal sigma = 0.0;
    PetscInt nconv = 0;
    SVDGetConverged(svd, &nconv);
    if (nconv < 1)
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "No singular triplet converged for which=%d; skip singular value extraction.\n",
            (int)which);
        sigma = PETSC_MAX_REAL;
    }
    else
    {
        SVDGetSingularTriplet(svd, 0, &sigma, NULL, NULL);
    }

    SVDDestroy(&svd);
    return sigma;
}

PetscReal ComputeExtremeEigenvalue(Mat A, EPSWhich which, PetscReal tol, PetscInt max_it)
{
    EPS eps;
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, A, NULL);
    EPSSetProblemType(eps, EPS_HEP);
    EPSSetType(eps, EPSKRYLOVSCHUR);
    EPSSetTolerances(eps, tol, max_it);
    if (which == EPS_SMALLEST_REAL)
    {
        /*
         * Krylov-Schur requires ncv <= nev + mpd. Since nev=1 here, use a
         * larger-but-valid search space for the smallest eigenvalue solve.
         */
        EPSSetDimensions(eps, 1, 31, 30);
        EPSSetTarget(eps, 0.0);
        EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL);

        ST st;
        EPSGetST(eps, &st);
        STSetType(st, STSINVERT);
        STSetShift(st, 0.0);
    }
    else
    {
        EPSSetDimensions(eps, 1, 31, 30);
        EPSSetWhichEigenpairs(eps, which);
    }

    EPSSetFromOptions(eps);
    EPSSolve(eps);

    EPSConvergedReason reason;
    EPSGetConvergedReason(eps, &reason);
    PetscPrintf(PETSC_COMM_WORLD, "EPS convergence reason: %d\n", reason);

    PetscReal lambda = 0.0;
    PetscInt nconv = 0;
    EPSGetConverged(eps, &nconv);
    if (nconv < 1)
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "No eigenpair converged for which=%d; skip eigenvalue extraction.\n",
            (int)which);
        lambda = PETSC_MAX_REAL;
    }
    else
    {
        PetscScalar eig;
        EPSGetEigenvalue(eps, 0, &eig, NULL);
        lambda = PetscRealPart(eig);
    }

    EPSDestroy(&eps);
    return lambda;
}

void MyMatMult(Mat M, Vec x, Vec y)
{
    UserCtx *ctx;
    Vec t;
    VecDuplicate(x, &t);

    MatShellGetContext(M, &ctx);
    MatMult(ctx->K, x, t);
    KSPSolve(ctx->ksp, t, y);
    VecDestroy(&t);
}

void MyMatMultTranspose(Mat M, Vec x, Vec y)
{
    UserCtx *ctx;
    Vec t;
    VecDuplicate(x, &t);

    MatShellGetContext(M, &ctx);
    KSPSolve(ctx->ksp, x, t);
    MatMult(ctx->K, t, y);
    VecDestroy(&t);
}

void MyMatMultK(Mat M, Vec x, Vec y)
{
    UserCtxK *ctx;

    MatShellGetContext(M, &ctx);
    KSPSolve(ctx->ksp, x, y);
}

int main(int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_x, part_num_y, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_x, part_num_y, dim, base_name);

    SlepcInitialize(&argc, &argv, NULL, NULL);

    PetscMPIInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    char spectral_method_name[16] = "eps";
    PetscBool method_set = PETSC_FALSE;
    PetscOptionsGetString(NULL, NULL, "-spectral_method",
        spectral_method_name, sizeof(spectral_method_name), &method_set);

    SpectralMethod spectral_method = SpectralMethod::EPS;
    if (std::string(spectral_method_name) == "svd")
    {
        spectral_method = SpectralMethod::SVD;
    }
    else if (std::string(spectral_method_name) == "eps")
    {
        spectral_method = SpectralMethod::EPS;
    }
    else
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
            "Unsupported -spectral_method value. Use 'svd' or 'eps'.");
    }

    PetscBool svd_skip_kfull = PETSC_FALSE;
    PetscOptionsGetBool(NULL, NULL, "-svd_skip_kfull", &svd_skip_kfull, NULL);

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
        std::cout << "spectral_method: " << spectral_method_name << std::endl;
        std::cout << "svd_skip_kfull: " << (svd_skip_kfull ? "true" : "false") << std::endl;
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

    std::vector<double> CP_fem;
    std::vector<int> ID_fem;
    std::vector<int> IEN_fem;
    int nlocalfunc_fem;
    int nlocalelemx_fem;
    int nlocalelemy_fem;

    std::string base_name_fem = "part_fem";
    std::string filename_fem = fm->GetPartitionFilename(base_name_fem, rank);
    fm->ReadPartition(filename_fem, nlocalfunc_fem,
        nlocalelemx_fem, nlocalelemy_fem,
        CP_fem, ID_fem, IEN_fem);

    ElementFEM * elem_fem = new ElementFEM(1, 1);
    LocalAssembly * locassem_fem = new LocalAssembly(1, 1);
    GlobalAssembly * globalassem_fem = new GlobalAssembly(IEN_fem, ID_fem, locassem_fem,
        4, nlocalfunc_fem, nlocalelemx_fem, nlocalelemy_fem);
    
    MatSetOption(globalassem_fem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    
    globalassem_fem->AssemStiffnessLoad(locassem_fem, IEN_fem, ID_fem, CP_fem, elem_fem);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling stiffness matrix and load vector...done\n");

    PetscReal tol = 1e-10;
    PetscInt max_it = 50000;
    if (spectral_method == SpectralMethod::SVD)
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "Using SVD spectral analysis for Kfull, Kfem, and K.\n");

        PetscReal rtol = 1e-10, abstol = PETSC_DEFAULT, divtol = PETSC_DEFAULT;
        PetscInt maxits = PETSC_DEFAULT;

        UserCtx ctx;
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, globalassem_fem->K, globalassem_fem->K);
        KSPSetFromOptions(ksp);
        KSPSetTolerances(ksp, rtol, abstol, divtol, maxits);
        ctx.ksp = ksp;
        ctx.K = globalassem->K;

        SVDConvergedReason reason;

        if (!svd_skip_kfull)
        {
            PetscInt mm, nn;
            MatGetSize(globalassem->K, &mm, &nn);

            Mat M;
            MatCreateShell(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, mm, nn, &ctx, &M);
            MatShellSetOperation(M, MATOP_MULT, (void(*)(void))MyMatMult);
            MatShellSetOperation(M, MATOP_MULT_TRANSPOSE, (void(*)(void))MyMatMultTranspose);

            SVD svd;
            SVDCreate(PETSC_COMM_WORLD, &svd);
            SVDSetOperators(svd, M, NULL);
            SVDSetProblemType(svd, SVD_STANDARD);
            SVDSetType(svd, SVDTRLANCZOS);
            SVDSetTolerances(svd, tol, max_it);
            SVDSetFromOptions(svd);
            SVDSetDimensions(svd, 1, PETSC_DEFAULT, PETSC_DEFAULT);
            SVDSetWhichSingularTriplets(svd, SVD_LARGEST);
            SVDSolve(svd);

            SVDGetConvergedReason(svd, &reason);
            PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

            PetscReal smax = 0.0;
            SVDGetSingularTriplet(svd, 0, &smax, NULL, NULL);
            PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value of Kfull: %.15g\n", smax);

            SVDDestroy(&svd);

            KSP ksp_inv;
            KSPCreate(PETSC_COMM_WORLD, &ksp_inv);
            KSPSetOperators(ksp_inv, globalassem->K, globalassem->K);
            KSPSetFromOptions(ksp_inv);
            KSPSetTolerances(ksp_inv, rtol, abstol, divtol, maxits);

            UserCtx ctx_inv;
            ctx_inv.ksp = ksp_inv;
            ctx_inv.K = globalassem_fem->K;
            Mat M_inv;
            MatCreateShell(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, mm, nn, &ctx_inv, &M_inv);
            MatShellSetOperation(M_inv, MATOP_MULT, (void(*)(void))MyMatMult);
            MatShellSetOperation(M_inv, MATOP_MULT_TRANSPOSE, (void(*)(void))MyMatMultTranspose);

            SVD svd_inv;
            SVDCreate(PETSC_COMM_WORLD, &svd_inv);
            SVDSetOperators(svd_inv, M_inv, NULL);
            SVDSetProblemType(svd_inv, SVD_STANDARD);
            SVDSetType(svd_inv, SVDTRLANCZOS);
            SVDSetTolerances(svd_inv, tol, max_it);
            SVDSetFromOptions(svd_inv);
            SVDSetDimensions(svd_inv, 1, PETSC_DEFAULT, PETSC_DEFAULT);
            SVDSetWhichSingularTriplets(svd_inv, SVD_LARGEST);
            SVDSolve(svd_inv);

            SVDGetConvergedReason(svd_inv, &reason);
            PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

            PetscReal smax_inv = 0.0;
            SVDGetSingularTriplet(svd_inv, 0, &smax_inv, NULL, NULL);
            PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of Kfull: %.15g\n", 1.0 / smax_inv);
            PetscPrintf(PETSC_COMM_WORLD, "Condition number of Kfull: %.15g\n", smax * smax_inv);

            SVDDestroy(&svd_inv);
            MatDestroy(&M_inv);
            KSPDestroy(&ksp_inv);
            MatDestroy(&M);
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD,
                "Skipping Kfull in SVD mode because -svd_skip_kfull true was requested.\n");
        }

        PetscInt mm_fem, nn_fem;
        MatGetSize(globalassem_fem->K, &mm_fem, &nn_fem);

        SVD svd_fem;
        SVDCreate(PETSC_COMM_WORLD, &svd_fem);
        SVDSetOperators(svd_fem, globalassem_fem->K, NULL);
        SVDSetProblemType(svd_fem, SVD_STANDARD);
        SVDSetType(svd_fem, SVDTRLANCZOS);
        SVDSetTolerances(svd_fem, tol, max_it);
        SVDSetFromOptions(svd_fem);
        SVDSetDimensions(svd_fem, 1, PETSC_DEFAULT, PETSC_DEFAULT);
        SVDSetWhichSingularTriplets(svd_fem, SVD_LARGEST);
        SVDSolve(svd_fem);

        SVDGetConvergedReason(svd_fem, &reason);
        PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

        PetscReal smax_fem = 0.0;
        SVDGetSingularTriplet(svd_fem, 0, &smax_fem, NULL, NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value of Kfem: %.15g\n", smax_fem);
        SVDDestroy(&svd_fem);

        KSP ksp_fem_inv;
        KSPCreate(PETSC_COMM_WORLD, &ksp_fem_inv);
        KSPSetOperators(ksp_fem_inv, globalassem_fem->K, globalassem_fem->K);
        KSPSetFromOptions(ksp_fem_inv);
        KSPSetTolerances(ksp_fem_inv, rtol, abstol, divtol, maxits);

        UserCtxK ctx_fem_inv;
        ctx_fem_inv.ksp = ksp_fem_inv;
        Mat M_fem_inv;
        MatCreateShell(PETSC_COMM_WORLD, nlocalfunc_fem, nlocalfunc_fem, mm_fem, nn_fem, &ctx_fem_inv, &M_fem_inv);
        MatShellSetOperation(M_fem_inv, MATOP_MULT, (void(*)(void))MyMatMultK);
        MatShellSetOperation(M_fem_inv, MATOP_MULT_TRANSPOSE, (void(*)(void))MyMatMultK);

        SVD svd_fem_inv;
        SVDCreate(PETSC_COMM_WORLD, &svd_fem_inv);
        SVDSetOperators(svd_fem_inv, M_fem_inv, NULL);
        SVDSetProblemType(svd_fem_inv, SVD_STANDARD);
        SVDSetType(svd_fem_inv, SVDTRLANCZOS);
        SVDSetTolerances(svd_fem_inv, tol, max_it);
        SVDSetFromOptions(svd_fem_inv);
        SVDSetDimensions(svd_fem_inv, 1, PETSC_DEFAULT, PETSC_DEFAULT);
        SVDSetWhichSingularTriplets(svd_fem_inv, SVD_LARGEST);
        SVDSolve(svd_fem_inv);

        SVDGetConvergedReason(svd_fem_inv, &reason);
        PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

        PetscInt nconv_fem_inv = 0;
        SVDGetConverged(svd_fem_inv, &nconv_fem_inv);
        if (nconv_fem_inv < 1)
        {
            PetscPrintf(PETSC_COMM_WORLD,
                "No singular triplet converged for inverse Kfem shell; skip singular value extraction.\n");
            PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of Kfem: unavailable\n");
            PetscPrintf(PETSC_COMM_WORLD, "Condition number of Kfem: unavailable\n");
        }
        else
        {
            PetscReal smax_fem_inv = 0.0;
            SVDGetSingularTriplet(svd_fem_inv, 0, &smax_fem_inv, NULL, NULL);
            PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of Kfem: %.15g\n", 1.0 / smax_fem_inv);
            PetscPrintf(PETSC_COMM_WORLD, "Condition number of Kfem: %.15g\n", smax_fem * smax_fem_inv);
        }

        SVDDestroy(&svd_fem_inv);
        MatDestroy(&M_fem_inv);
        KSPDestroy(&ksp_fem_inv);

        PetscInt mm, nn;
        MatGetSize(globalassem->K, &mm, &nn);

        SVD svd_K;
        SVDCreate(PETSC_COMM_WORLD, &svd_K);
        SVDSetOperators(svd_K, globalassem->K, NULL);
        SVDSetProblemType(svd_K, SVD_STANDARD);
        SVDSetType(svd_K, SVDTRLANCZOS);
        SVDSetTolerances(svd_K, tol, max_it);
        SVDSetFromOptions(svd_K);
        SVDSetDimensions(svd_K, 1, PETSC_DEFAULT, PETSC_DEFAULT);
        SVDSetWhichSingularTriplets(svd_K, SVD_LARGEST);
        SVDSolve(svd_K);

        SVDGetConvergedReason(svd_K, &reason);
        PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

        PetscReal smax_K = 0.0;
        SVDGetSingularTriplet(svd_K, 0, &smax_K, NULL, NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value of K: %.15g\n", smax_K);
        SVDDestroy(&svd_K);

        KSP ksp_K_inv;
        KSPCreate(PETSC_COMM_WORLD, &ksp_K_inv);
        KSPSetOperators(ksp_K_inv, globalassem->K, globalassem->K);
        KSPSetFromOptions(ksp_K_inv);
        KSPSetTolerances(ksp_K_inv, rtol, abstol, divtol, maxits);

        UserCtxK ctx_K_inv;
        ctx_K_inv.ksp = ksp_K_inv;
        Mat M_K_inv;
        MatCreateShell(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, mm, nn, &ctx_K_inv, &M_K_inv);
        MatShellSetOperation(M_K_inv, MATOP_MULT, (void(*)(void))MyMatMultK);
        MatShellSetOperation(M_K_inv, MATOP_MULT_TRANSPOSE, (void(*)(void))MyMatMultK);

        SVD svd_K_inv;
        SVDCreate(PETSC_COMM_WORLD, &svd_K_inv);
        SVDSetOperators(svd_K_inv, M_K_inv, NULL);
        SVDSetProblemType(svd_K_inv, SVD_STANDARD);
        SVDSetType(svd_K_inv, SVDTRLANCZOS);
        SVDSetTolerances(svd_K_inv, tol, max_it);
        SVDSetFromOptions(svd_K_inv);
        SVDSetDimensions(svd_K_inv, 1, PETSC_DEFAULT, PETSC_DEFAULT);
        SVDSetWhichSingularTriplets(svd_K_inv, SVD_LARGEST);
        SVDSolve(svd_K_inv);

        SVDGetConvergedReason(svd_K_inv, &reason);
        PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %d\n", reason);

        PetscInt nconv_K_inv = 0;
        SVDGetConverged(svd_K_inv, &nconv_K_inv);
        if (nconv_K_inv < 1)
        {
            PetscPrintf(PETSC_COMM_WORLD,
                "No singular triplet converged for inverse K shell; skip singular value extraction.\n");
            PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of K: unavailable\n");
            PetscPrintf(PETSC_COMM_WORLD, "Condition number of K: unavailable\n");
        }
        else
        {
            PetscReal smax_K_inv = 0.0;
            SVDGetSingularTriplet(svd_K_inv, 0, &smax_K_inv, NULL, NULL);
            PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value of K: %.15g\n", 1.0 / smax_K_inv);
            PetscPrintf(PETSC_COMM_WORLD, "Condition number of K: %.15g\n", smax_K * smax_K_inv);
        }

        SVDDestroy(&svd_K_inv);
        MatDestroy(&M_K_inv);
        KSPDestroy(&ksp_K_inv);

        KSPDestroy(&ksp);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "Using EPS spectral analysis for explicit Kfem and K.\n");

        Mat K_neg = nullptr;
        Mat Kfem_neg = nullptr;
        MatDuplicate(globalassem->K, MAT_COPY_VALUES, &K_neg);
        MatDuplicate(globalassem_fem->K, MAT_COPY_VALUES, &Kfem_neg);
        MatScale(K_neg, -1.0);
        MatScale(Kfem_neg, -1.0);
        PetscPrintf(PETSC_COMM_WORLD,
            "Negated explicit matrices K and Kfem before EPS analysis.\n");

        AnalyzeWithEPS(Kfem_neg, "Kfem", tol, max_it);
        AnalyzeWithEPS(K_neg, "K", tol, max_it);

        PetscPrintf(PETSC_COMM_WORLD,
            "EPS mode does not analyze Kfull because the current EPS route is only implemented for explicit matrices.\n");

        MatDestroy(&K_neg);
        MatDestroy(&Kfem_neg);
    }
    
    delete fm; fm = nullptr;
    delete elem; elem = nullptr;
    delete locassem; locassem = nullptr;
    delete globalassem; globalassem = nullptr;
    delete elem_fem; elem_fem = nullptr;
    delete locassem_fem; locassem_fem = nullptr;
    delete globalassem_fem; globalassem_fem = nullptr;
    
    SlepcFinalize();
    return 0;
}
