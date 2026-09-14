#include "FileManager.hpp"
#include "InvLM.hpp"

int main(int argc, char **argv)
{
    FileManager * fm = new FileManager();

    std::string file_info = "info.txt";

    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    fm->ReadPreprocessInfo(file_info,
        p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    std::cout << "Preprocessing information of matrix-free:" << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "Lx: " << Lx << std::endl;
    std::cout << "Ly: " << Ly << std::endl;
    std::cout << "nElemX: " << nElemX << std::endl;
    std::cout << "nElemY: " << nElemY << std::endl;
    std::cout << "part_num_1d: " << part_num_1d << std::endl;
    std::cout << "dim: " << dim << std::endl;
    std::cout << "base_name: " << base_name << std::endl;

    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;
    std::vector<int> ID{};
    std::vector<int> Dir{};
    std::vector<int> IEN{};
    std::vector<double> CP{};
    std::vector<double> elem_size1{};
    std::vector<double> elem_size2{};
    std::vector<double> NURBSExtraction1{};
    std::vector<double> NURBSExtraction2{};

    std::string filename = fm->GetPartitionFilename(base_name, 0);
    fm->ReadPartition(filename,
        nlocalfunc,
        nlocalelemx,
        nlocalelemy,
        elem_size1,
        elem_size2,
        CP,
        ID,
        Dir,
        IEN,
        NURBSExtraction1,
        NURBSExtraction2);
    
    InvLM * invlm = new InvLM(
        (p + 1) * (q + 1),
        nlocalelemx * nlocalelemy,
        nlocalfunc,
        IEN);
    
    std::vector<int> invlm_elemnum = invlm->GetAllElemNum();
    std::vector<int> invlm_offset = invlm->GetAllOffset();
    std::cout << "InvLM elemNum:" << std::endl;
    for (size_t i = 0; i < invlm_elemnum.size(); ++i)
    {
        std::cout << "Function " << i << ": " << invlm_elemnum[i] << std::endl;
    }
    std::cout << "InvLM offset:" << std::endl;
    for (size_t i = 0; i < invlm_offset.size(); ++i)
    {
        std::cout << "Function " << i << ": " << invlm_offset[i] << std::endl;
    }
    std::vector<int> invlm_elemidx = invlm->GetAllElemIdx();
    std::vector<int> invlm_baseidx = invlm->GetAllBaseIdx();
    std::cout << "InvLM elemIdx and baseIdx:" << std::endl;
    for (int i = 0; i < nlocalfunc; ++i)
    {
        std::cout << "Function " << i << ":" << std::endl;
        std::vector<int> elemidx = invlm->GetElemIdx(i);
        std::vector<int> baseidx = invlm->GetBaseIdx(i);
        for (size_t j = 0; j < elemidx.size(); ++j)
        {
            std::cout << "  Element " << elemidx[j]
                      << ", Local Basis " << baseidx[j] << std::endl;
        }
    }

    return 0;
}