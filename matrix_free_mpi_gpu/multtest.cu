#include <petscmat.h>
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
    std::cout << "p: " << p << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "Lx: " << Lx << std::endl;
    std::cout << "Ly: " << Ly << std::endl;
    std::cout << "nElemX: " << nElemX << std::endl;
    std::cout << "nElemY: " << nElemY << std::endl;
    std::cout << "part_num_1d: " << part_num_1d << std::endl;
    std::cout << "dim: " << dim << std::endl;
    std::cout << "base_name: " << base_name << std::endl;

    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> localID;
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
        CP, ID, localID, ghostID, Dir, IEN, NURBSExtraction1, NURBSExtraction2);

    ElementMF * elemmf = new ElementMF(p, q);
    int nLocBas = elemmf->GetNumLocalBasis();
    GlobalAssemblyMF * globalAssembly = new GlobalAssemblyMF(
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);

    QuadraturePoint * quad1 = new QuadraturePoint(p+1, 0, 1);
    QuadraturePoint * quad2 = new QuadraturePoint(q+1, 0, 1);

    BernsteinBasis * bernstein = new BernsteinBasis(p);

    globalAssembly->AssemLoad(quad1, quad2,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein);
    VecView(globalAssembly->F, PETSC_VIEWER_STDOUT_WORLD);
    
    Vec x;
    VecDuplicate(globalAssembly->F, &x);
    VecCopy(globalAssembly->F, x);

    Vec y;
    VecDuplicate(globalAssembly->F, &y);

    globalAssembly->MatMulMF(quad1, quad2,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein,
        x, y);
    
    VecView(y, PETSC_VIEWER_STDOUT_WORLD);

    VecCopy(y,x);
    globalAssembly->MatMulMF(quad1, quad2,
        IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elemmf, bernstein,
        x, y);
    VecView(y, PETSC_VIEWER_STDOUT_WORLD);

    VecDestroy(&x);
    VecDestroy(&y);
    delete globalAssembly;
    delete elemmf;
    delete quad1;
    delete quad2;

    delete fm;
    PetscFinalize();
    return 0;
}
