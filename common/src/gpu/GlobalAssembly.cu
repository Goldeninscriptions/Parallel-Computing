#include "GlobalAssembly.cuh"

GlobalAssembly::GlobalAssembly(const int &nLocBas,
    const int &nnz, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy,
    const std::vector<int> &rows,
    const std::vector<int> &cols)
: nLocBas(nlocalfunc), nnz(nnz), nlocalfunc(nlocalfunc),
  nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{
    PetscInt * d_rows, * d_cols;
    MallocDeviceMemory(&d_rows, nnz);
    MallocDeviceMemory(&d_cols, nnz);
    CopyToDevice(d_rows, rows.data(), nnz);
    CopyToDevice(d_cols, cols.data(), nnz);
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, nlocalfunc, nlocalfunc, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(K, MATAIJCUSPARSE);
    MatSetPreallocationCOO(K, nnz, d_rows, d_cols);
    FreeDeviceMemory(d_rows);
    FreeDeviceMemory(d_cols);
}

GlobalAssembly::~GlobalAssembly()
{
    MatDestroy(&K);
}

void GlobalAssembly::AssemStiffness(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &dir2coo,
    const std::vector<double> &CP,
    const std::vector<int> &elem2coo,
    ElementFEM * const &elem)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> sqp_N{};
    std::vector<double> sqp_dN_dxi{};
    std::vector<double> sqp_dN_deta{};

    std::vector<double> N{};
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    std::vector<double> weight{};

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateRefElementSingleQP(qp1[ii], qp2[jj], sqp_N, sqp_dN_dxi, sqp_dN_deta);
            N.insert(N.end(), sqp_N.begin(),sqp_N.end());
            dN_dxi.insert(dN_dxi.end(), sqp_dN_dxi.begin(), sqp_dN_dxi.end());
            dN_deta.insert(dN_deta.end(), sqp_dN_deta.begin(), sqp_dN_deta.end());
            weight.push_back(w1[ii] * w2[jj]);
        }
    }

    double * d_N, * d_dN_dxi, * d_dN_deta, * d_weight;
    MallocDeviceMemory(&d_N, N.size());
    MallocDeviceMemory(&d_dN_dxi, dN_dxi.size());
    MallocDeviceMemory(&d_dN_deta, dN_deta.size());
    MallocDeviceMemory(&d_weight, weight.size());
    CopyToDevice(d_N, N.data(), N.size());
    CopyToDevice(d_dN_dxi, dN_dxi.data(), dN_dxi.size());
    CopyToDevice(d_dN_deta, dN_deta.data(), dN_deta.size());
    CopyToDevice(d_weight, weight.data(), weight.size());

    int * d_IEN, * d_dir2coo, * d_elem2coo;
    double * d_CP, * d_val;
    MallocDeviceMemory(&d_IEN, IEN.size());
    MallocDeviceMemory(&d_dir2coo, dir2coo.size());
    MallocDeviceMemory(&d_CP, CP.size());
    MallocDeviceMemory(&d_elem2coo, elem2coo.size());
    MallocDeviceMemory(&d_val, nnz);
    CopyToDevice(d_IEN, IEN.data(), IEN.size());
    CopyToDevice(d_dir2coo, dir2coo.data(), dir2coo.size());
    CopyToDevice(d_CP, CP.data(), CP.size());
    CopyToDevice(d_elem2coo, elem2coo.data(), elem2coo.size());
    cudaMemset(d_val, 0, nnz * sizeof(double));

    AssembleStiffnessCUDA(nqp1, nqp2,
        nlocalelemx, nlocalelemy,
        d_N, d_dN_dxi, d_dN_deta,
        d_weight, d_IEN,
        d_CP, d_elem2coo, d_val);

    DirichletBCKCUDA(d_dir2coo, dir2coo.size(), d_val);

    MatSetValuesCOO(K, d_val, INSERT_VALUES);

    FreeDeviceMemory(d_N);
    FreeDeviceMemory(d_dN_dxi);
    FreeDeviceMemory(d_dN_deta);
    FreeDeviceMemory(d_weight);
    FreeDeviceMemory(d_IEN);
    FreeDeviceMemory(d_dir2coo);
    FreeDeviceMemory(d_CP);
    FreeDeviceMemory(d_elem2coo);
    FreeDeviceMemory(d_val);
}