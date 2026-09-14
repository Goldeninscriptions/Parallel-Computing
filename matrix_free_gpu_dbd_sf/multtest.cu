#include <petscmat.h>
#include <cuda_runtime.h>
#include <fstream>
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
    PetscInt setup_repeat = 20;
    PetscBool setup_repeat_set = PETSC_FALSE;
    PetscInt matmult_repeat = 20;
    PetscBool matmult_repeat_set = PETSC_FALSE;
    char dump_vec_path[PETSC_MAX_PATH_LEN] = "";
    PetscBool dump_vec_set = PETSC_FALSE;
    PetscOptionsGetInt(NULL, NULL, "-repeat", &repeat, &repeat_set);
    PetscOptionsGetInt(NULL, NULL, "-setup_repeat", &setup_repeat, &setup_repeat_set);
    PetscOptionsGetInt(NULL, NULL, "-matmult_repeat", &matmult_repeat, &matmult_repeat_set);
    PetscOptionsGetString(NULL, NULL, "-dump_vec", dump_vec_path, sizeof(dump_vec_path), &dump_vec_set);
    if (!setup_repeat_set)
        setup_repeat = repeat;
    if (!matmult_repeat_set)
        matmult_repeat = repeat;
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
    std::cout << "setup_repeat: " << setup_repeat << std::endl;
    std::cout << "matmult_repeat: " << matmult_repeat << std::endl;

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
        CP, ID, Dir, IEN, NURBSExtraction1, NURBSExtraction2);

    ElementMF * elemmf = new ElementMF(p, q);
    int nLocBas = elemmf->GetNumLocalBasis();
    GlobalAssemblyMF * globalAssembly = nullptr;

    InvLM * invlm = new InvLM((p+1)*(q+1), nlocalelemx*nlocalelemy, nlocalfunc, IEN);

    std::vector<int> invlm_elemNum = invlm->GetAllElemNum();
    std::vector<int> invlm_offset = invlm->GetAllOffset();
    std::vector<int> invlm_elemIdx = invlm->GetAllElemIdx();
    std::vector<int> invlm_baseIdx = invlm->GetAllBaseIdx();

     std::vector<int> xelemIdx{};
    std::vector<int> yelemIdx{};
    for (int ey = 0; ey < nlocalelemy; ++ey)
    {
        for (int ex = 0; ex < nlocalelemx; ++ex)
        {
            xelemIdx.push_back(ex);
            yelemIdx.push_back(ey);
        }
    }
    QuadraturePoint * quad1 = new QuadraturePoint(p+1, 0, 1);
    QuadraturePoint * quad2 = new QuadraturePoint(q+1, 0, 1);

    BernsteinBasis * bernstein = new BernsteinBasis(p);
    cudaEvent_t setup_event_start;
    cudaEvent_t setup_event_stop;
    cudaEvent_t matmult_event_start;
    cudaEvent_t matmult_event_stop;
    cudaEventCreate(&setup_event_start);
    cudaEventCreate(&setup_event_stop);
    cudaEventCreate(&matmult_event_start);
    cudaEventCreate(&matmult_event_stop);

    PetscLogDouble setup_start = 0.0;
    PetscLogDouble setup_end = 0.0;
    float setup_cuda_total_ms = 0.0f;
    PetscTime(&setup_start);
    for (PetscInt iter = 0; iter < setup_repeat; ++iter)
    {
        delete globalAssembly;
        globalAssembly = new GlobalAssemblyMF(
            nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);
        cudaEventRecord(setup_event_start);
        globalAssembly->AssemLoad(quad1, quad2,
            IEN, ID, Dir,
            invlm_elemNum, invlm_offset, invlm_elemIdx, invlm_baseIdx,
            xelemIdx, yelemIdx,
            CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elemmf, bernstein);
        cudaEventRecord(setup_event_stop);
        cudaEventSynchronize(setup_event_stop);
        float setup_cuda_iter_ms = 0.0f;
        cudaEventElapsedTime(&setup_cuda_iter_ms, setup_event_start, setup_event_stop);
        setup_cuda_total_ms += setup_cuda_iter_ms;
    }
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
        IEN, ID, Dir,
        invlm_elemNum, invlm_offset, invlm_elemIdx, invlm_baseIdx,
        xelemIdx, yelemIdx,
        CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein,
        x, y);
    cudaDeviceSynchronize();

    PetscLogDouble mult_start = 0.0;
    PetscLogDouble mult_end = 0.0;
    PetscLogDouble total_start = setup_start;
    PetscLogDouble total_end = 0.0;
    float matmult_cuda_total_ms = 0.0f;
    PetscTime(&mult_start);
    for (PetscInt iter = 0; iter < matmult_repeat; ++iter)
    {
        cudaEventRecord(matmult_event_start);
        globalAssembly->MatMulMF(quad1, quad2,
            IEN, ID, Dir,
            invlm_elemNum, invlm_offset, invlm_elemIdx, invlm_baseIdx,
            xelemIdx, yelemIdx,
            CP,
            NURBSExtraction1, NURBSExtraction2,
            elem_size1, elem_size2, elemmf, bernstein,
            x, y);
        cudaEventRecord(matmult_event_stop);
        cudaEventSynchronize(matmult_event_stop);
        float matmult_cuda_iter_ms = 0.0f;
        cudaEventElapsedTime(&matmult_cuda_iter_ms, matmult_event_start, matmult_event_stop);
        matmult_cuda_total_ms += matmult_cuda_iter_ms;
    }
    cudaDeviceSynchronize();
    PetscTime(&mult_end);
    total_end = mult_end;

    PetscReal ynorm = 0.0;
    VecNorm(y, NORM_2, &ynorm);
    PetscPrintf(PETSC_COMM_WORLD, "Setup total time: %.6f\n", setup_end - setup_start);
    PetscPrintf(PETSC_COMM_WORLD, "Setup average time: %.6f\n",
        (setup_end - setup_start) / setup_repeat);
    PetscPrintf(PETSC_COMM_WORLD, "Setup CUDA total time (ms): %.6f\n", setup_cuda_total_ms);
    PetscPrintf(PETSC_COMM_WORLD, "Setup CUDA average time (ms): %.6f\n",
        setup_cuda_total_ms / setup_repeat);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult total time: %.6f\n", mult_end - mult_start);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult average time: %.6f\n",
        (mult_end - mult_start) / matmult_repeat);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult CUDA total time (ms): %.6f\n", matmult_cuda_total_ms);
    PetscPrintf(PETSC_COMM_WORLD, "MatMult CUDA average time (ms): %.6f\n",
        matmult_cuda_total_ms / matmult_repeat);
    PetscPrintf(PETSC_COMM_WORLD, "Total time: %.6f\n", total_end - total_start);
    PetscPrintf(PETSC_COMM_WORLD, "Output vector 2-norm: %.12e\n", ynorm);
    if (dump_vec_set)
    {
        PetscViewer viewer = nullptr;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, dump_vec_path, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU);
        VecView(y, viewer);
        PetscViewerDestroy(&viewer);
        PetscPrintf(PETSC_COMM_WORLD, "Dumped output vector to %s\n", dump_vec_path);
    }

    delete globalAssembly;
    delete elemmf;
    delete quad1;
    delete quad2;

    delete fm;
    VecDestroy(&x);
    VecDestroy(&y);
    cudaEventDestroy(setup_event_start);
    cudaEventDestroy(setup_event_stop);
    cudaEventDestroy(matmult_event_start);
    cudaEventDestroy(matmult_event_stop);
    PetscFinalize();
    return 0;
}
