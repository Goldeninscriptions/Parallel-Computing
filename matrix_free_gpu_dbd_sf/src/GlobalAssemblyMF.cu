#include "GlobalAssemblyMF.cuh"

GlobalAssemblyMF::GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{   
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DETERMINE);
    VecSetType(F, VECCUDA);
    VecSet(F, 0.0);
}

void GlobalAssemblyMF::AssemLoad(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<int> &invlm_elemNum,
    const std::vector<int> &invlm_offset,
    const std::vector<int> &invlm_elemIdx,
    const std::vector<int> &invlm_baseIdx,
    const std::vector<int> &xelemIdx,
    const std::vector<int> &yelemIdx,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();
    const int pp = elemmf->get_p();
    const int qq = elemmf->get_q();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};
    std::vector<double> temp{};

    for (int i = 0; i < nqp1; ++i)
    {
        temp = bernstein->GetBernsteinBasisSingleQP(qp1[i]);
        B1.insert(B1.end(), temp.begin(), temp.end());
        temp = bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]);
        dB1.insert(dB1.end(), temp.begin(), temp.end());
    }
    for (int j = 0; j < nqp2; ++j)
    {
        temp = bernstein->GetBernsteinBasisSingleQP(qp2[j]);
        B2.insert(B2.end(), temp.begin(), temp.end());
        temp = bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]);
        dB2.insert(dB2.end(), temp.begin(), temp.end());
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    MallocDeviceMemory(&d_B1, B1.size());
    MallocDeviceMemory(&d_B2, B2.size());
    MallocDeviceMemory(&d_dB1, dB1.size());
    MallocDeviceMemory(&d_dB2, dB2.size());
    CopyToDevice(d_B1, B1.data(), B1.size());
    CopyToDevice(d_B2, B2.data(), B2.size());
    CopyToDevice(d_dB1, dB1.data(), dB1.size());
    CopyToDevice(d_dB2, dB2.data(), dB2.size());

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    MallocDeviceMemory(&d_NURBSExtraction1, NURBSExtraction1.size());
    MallocDeviceMemory(&d_NURBSExtraction2, NURBSExtraction2.size());
    CopyToDevice(d_NURBSExtraction1, NURBSExtraction1.data(), NURBSExtraction1.size());
    CopyToDevice(d_NURBSExtraction2, NURBSExtraction2.data(), NURBSExtraction2.size());

    int * d_IEN, * d_ID;
    double * d_CP;
    MallocDeviceMemory(&d_IEN, IEN.size());
    MallocDeviceMemory(&d_ID, ID.size());
    MallocDeviceMemory(&d_CP, CP.size());
    CopyToDevice(d_IEN, IEN.data(), IEN.size());
    CopyToDevice(d_ID, ID.data(), ID.size());
    CopyToDevice(d_CP, CP.data(), CP.size());

    double * d_elem_size1, * d_elem_size2;
    MallocDeviceMemory(&d_elem_size1, elem_size1.size());
    MallocDeviceMemory(&d_elem_size2, elem_size2.size());
    CopyToDevice(d_elem_size1, elem_size1.data(), elem_size1.size());
    CopyToDevice(d_elem_size2, elem_size2.data(), elem_size2.size());

    double * qw1, * qw2;
    MallocDeviceMemory(&qw1, nqp1);
    MallocDeviceMemory(&qw2, nqp2);
    CopyToDevice(qw1, w1.data(), nqp1);
    CopyToDevice(qw2, w2.data(), nqp2);

    int * d_invlm_elemNum, * d_invlm_offset;
    int * d_invlm_elemIdx, * d_invlm_baseIdx;
    MallocDeviceMemory(&d_invlm_elemNum, invlm_elemNum.size());
    MallocDeviceMemory(&d_invlm_offset, invlm_offset.size());
    MallocDeviceMemory(&d_invlm_elemIdx, invlm_elemIdx.size());
    MallocDeviceMemory(&d_invlm_baseIdx, invlm_baseIdx.size());
    CopyToDevice(d_invlm_elemNum, invlm_elemNum.data(), invlm_elemNum.size());
    CopyToDevice(d_invlm_offset, invlm_offset.data(), invlm_offset.size());
    CopyToDevice(d_invlm_elemIdx, invlm_elemIdx.data(), invlm_elemIdx.size());
    CopyToDevice(d_invlm_baseIdx, invlm_baseIdx.data(), invlm_baseIdx.size());

    int nfunc = static_cast<int>(invlm_elemNum.size());

    int * d_xelemIdx, * d_yelemIdx;
    MallocDeviceMemory(&d_xelemIdx, xelemIdx.size());
    MallocDeviceMemory(&d_yelemIdx, yelemIdx.size());
    CopyToDevice(d_xelemIdx, xelemIdx.data(), xelemIdx.size());
    CopyToDevice(d_yelemIdx, yelemIdx.data(), yelemIdx.size());

    int * d_Dir = nullptr;
    if (!Dir.empty())
    {
        MallocDeviceMemory(&d_Dir, Dir.size());
        CopyToDevice(d_Dir, Dir.data(), Dir.size());
    }

    VecSet(F, 0.0);
    double *d_F_array;
    VecCUDAGetArray(F, &d_F_array);

    AssembleLoadCUDA(nfunc, pp, qq,
        nlocalelemx, nlocalelemy,
        d_B1, d_B2, d_dB1, d_dB2,
        d_NURBSExtraction1, d_NURBSExtraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_invlm_elemNum, d_invlm_offset,
        d_invlm_elemIdx, d_invlm_baseIdx,
        d_xelemIdx, d_yelemIdx, 
        d_F_array);
    
    if (d_Dir != nullptr)
        DirichletBCCUDA(d_Dir, static_cast<int>(Dir.size()), d_F_array, 0.0);

    VecCUDARestoreArray(F, &d_F_array);

    FreeDeviceMemory(d_B1);
    FreeDeviceMemory(d_B2);
    FreeDeviceMemory(d_dB1);
    FreeDeviceMemory(d_dB2);
    FreeDeviceMemory(d_NURBSExtraction1);
    FreeDeviceMemory(d_NURBSExtraction2);
    FreeDeviceMemory(d_IEN);
    FreeDeviceMemory(d_ID);
    FreeDeviceMemory(d_CP);
    FreeDeviceMemory(d_elem_size1);
    FreeDeviceMemory(d_elem_size2);
    FreeDeviceMemory(qw1);
    FreeDeviceMemory(qw2);
    FreeDeviceMemory(d_invlm_elemNum);
    FreeDeviceMemory(d_invlm_offset);
    FreeDeviceMemory(d_invlm_elemIdx);
    FreeDeviceMemory(d_invlm_baseIdx);
    FreeDeviceMemory(d_xelemIdx);
    FreeDeviceMemory(d_yelemIdx);
    FreeDeviceMemory(d_Dir);
}

void GlobalAssemblyMF::MatMulMF(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<int> &invlm_elemNum,
    const std::vector<int> &invlm_offset,
    const std::vector<int> &invlm_elemIdx,
    const std::vector<int> &invlm_baseIdx,
    const std::vector<int> &xelemIdx,
    const std::vector<int> &yelemIdx,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein,
    Vec x, Vec y)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();
    const int pp = elemmf->get_p();
    const int qq = elemmf->get_q();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};
    std::vector<double> temp{};

    for (int i = 0; i < nqp1; ++i)
    {
        temp = bernstein->GetBernsteinBasisSingleQP(qp1[i]);
        B1.insert(B1.end(), temp.begin(), temp.end());
        temp = bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]);
        dB1.insert(dB1.end(), temp.begin(), temp.end());
    }
    for (int j = 0; j < nqp2; ++j)
    {
        temp = bernstein->GetBernsteinBasisSingleQP(qp2[j]);
        B2.insert(B2.end(), temp.begin(), temp.end());
        temp = bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]);
        dB2.insert(dB2.end(), temp.begin(), temp.end());
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    MallocDeviceMemory(&d_B1, B1.size());
    MallocDeviceMemory(&d_B2, B2.size());
    MallocDeviceMemory(&d_dB1, dB1.size());
    MallocDeviceMemory(&d_dB2, dB2.size());
    CopyToDevice(d_B1, B1.data(), B1.size());
    CopyToDevice(d_B2, B2.data(), B2.size());
    CopyToDevice(d_dB1, dB1.data(), dB1.size());
    CopyToDevice(d_dB2, dB2.data(), dB2.size());

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    MallocDeviceMemory(&d_NURBSExtraction1, NURBSExtraction1.size());
    MallocDeviceMemory(&d_NURBSExtraction2, NURBSExtraction2.size());
    CopyToDevice(d_NURBSExtraction1, NURBSExtraction1.data(),
        NURBSExtraction1.size());
    CopyToDevice(d_NURBSExtraction2, NURBSExtraction2.data(),
        NURBSExtraction2.size());
    
    int * d_IEN, * d_ID;
    double * d_CP;
    MallocDeviceMemory(&d_IEN, IEN.size());
    MallocDeviceMemory(&d_ID, ID.size());
    MallocDeviceMemory(&d_CP, CP.size());
    CopyToDevice(d_IEN, IEN.data(), IEN.size());
    CopyToDevice(d_ID, ID.data(), ID.size());
    CopyToDevice(d_CP, CP.data(), CP.size());

    double * d_elem_size1, * d_elem_size2;
    MallocDeviceMemory(&d_elem_size1, elem_size1.size());
    MallocDeviceMemory(&d_elem_size2, elem_size2.size());
    CopyToDevice(d_elem_size1, elem_size1.data(), elem_size1.size());
    CopyToDevice(d_elem_size2, elem_size2.data(), elem_size2.size());

    double * qw1, * qw2;
    MallocDeviceMemory(&qw1, nqp1);
    MallocDeviceMemory(&qw2, nqp2);
    CopyToDevice(qw1, w1.data(), nqp1);
    CopyToDevice(qw2, w2.data(), nqp2);

    int * d_invlm_elemNum, * d_invlm_offset;
    int * d_invlm_elemIdx, * d_invlm_baseIdx;
    MallocDeviceMemory(&d_invlm_elemNum, invlm_elemNum.size());
    MallocDeviceMemory(&d_invlm_offset, invlm_offset.size());
    MallocDeviceMemory(&d_invlm_elemIdx, invlm_elemIdx.size());
    MallocDeviceMemory(&d_invlm_baseIdx, invlm_baseIdx.size());
    CopyToDevice(d_invlm_elemNum, invlm_elemNum.data(), invlm_elemNum.size());
    CopyToDevice(d_invlm_offset, invlm_offset.data(), invlm_offset.size());
    CopyToDevice(d_invlm_elemIdx, invlm_elemIdx.data(), invlm_elemIdx.size());
    CopyToDevice(d_invlm_baseIdx, invlm_baseIdx.data(), invlm_baseIdx.size());

    int nfunc = static_cast<int>(invlm_elemNum.size());

    int * d_xelemIdx, * d_yelemIdx;
    MallocDeviceMemory(&d_xelemIdx, xelemIdx.size());
    MallocDeviceMemory(&d_yelemIdx, yelemIdx.size());
    CopyToDevice(d_xelemIdx, xelemIdx.data(), xelemIdx.size());
    CopyToDevice(d_yelemIdx, yelemIdx.data(), yelemIdx.size());

    int * d_Dir = nullptr;
    if (!Dir.empty())
    {
        MallocDeviceMemory(&d_Dir, Dir.size());
        CopyToDevice(d_Dir, Dir.data(), Dir.size());
    }

    VecSet(y, 0.0);
    const double *d_x_array;
    double *d_y_array;
    VecCUDAGetArrayRead(x, &d_x_array);
    VecCUDAGetArray(y, &d_y_array);

    MatrixFreeMatMultCUDA(nfunc, pp, qq,
        nlocalelemx, nlocalelemy,
        d_B1, d_B2, d_dB1, d_dB2,
        d_NURBSExtraction1, d_NURBSExtraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_invlm_elemNum, d_invlm_offset,
        d_invlm_elemIdx, d_invlm_baseIdx, 
        d_xelemIdx, d_yelemIdx,
        d_x_array, d_y_array);

    if (d_Dir != nullptr)
        DirichletBCCUDA(d_Dir, static_cast<int>(Dir.size()), d_y_array, 0.0);
    
    VecCUDARestoreArrayRead(x, &d_x_array);
    VecCUDARestoreArray(y, &d_y_array);

    FreeDeviceMemory(d_B1);
    FreeDeviceMemory(d_B2);
    FreeDeviceMemory(d_dB1);
    FreeDeviceMemory(d_dB2);
    FreeDeviceMemory(d_NURBSExtraction1);
    FreeDeviceMemory(d_NURBSExtraction2);
    FreeDeviceMemory(d_IEN);
    FreeDeviceMemory(d_ID);
    FreeDeviceMemory(d_CP);
    FreeDeviceMemory(d_elem_size1);
    FreeDeviceMemory(d_elem_size2);
    FreeDeviceMemory(qw1);
    FreeDeviceMemory(qw2);
    FreeDeviceMemory(d_invlm_elemNum);
    FreeDeviceMemory(d_invlm_offset);
    FreeDeviceMemory(d_invlm_elemIdx);
    FreeDeviceMemory(d_invlm_baseIdx);
    FreeDeviceMemory(d_xelemIdx);
    FreeDeviceMemory(d_yelemIdx);
    FreeDeviceMemory(d_Dir);
}
