#include <petscmat.h>
#include "FileManager.hpp"
#include "memory.cuh"

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

    int nfunc = 0;
    for (int ii = 0; ii < part_num_1d * part_num_1d; ++ii)
    {
        int nlocalfunc = 0;
        std::string base_name_fem = "part_fem";
        std::string filename_fem = fm->GetPartitionFilename(base_name_fem, ii);
        fm->ReadPartition(filename_fem, nlocalfunc);
        nfunc += nlocalfunc;
    }

    std::string filename_dr = "coordinate";

    std::vector<int> rows{};
    std::vector<int> cols{};
    int nnz = 0;
    for (int ii = 0; ii < part_num_1d * part_num_1d; ++ii)
    {
        std::vector<int> rows_temp{};
        std::vector<int> cols_temp{};
        int temp_nnz = 0;
        std::string filename = fm->GetNonZeroCoordinateFilename(filename_dr, ii);
        fm->ReadNonZeroCoordinate(filename, temp_nnz, rows_temp, cols_temp);
        rows.insert(rows.end(), rows_temp.begin(), rows_temp.end());
        cols.insert(cols.end(), cols_temp.begin(), cols_temp.end());
        nnz += temp_nnz;
    }

    std::string map_name = "new_to_old_mapping.txt";
    std::vector<int> new_to_old;
    fm->ReadNewToOldMapping(map_name, new_to_old);

    for (int i = 0; i < rows.size(); ++i)
    {
        rows[i] = new_to_old[rows[i]];
        cols[i] = new_to_old[cols[i]];
    }

    Mat K;
    PetscInt * d_rows, * d_cols;
    MallocDeviceMemory(&d_rows, nnz);
    MallocDeviceMemory(&d_cols, nnz);
    CopyToDevice(d_rows, rows.data(), nnz);
    CopyToDevice(d_cols, cols.data(), nnz);
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetType(K, MATAIJCUSPARSE);
    MatSetSizes(K, nfunc, nfunc, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetPreallocationCOO(K, nnz, d_rows, d_cols);
    FreeDeviceMemory(d_rows);
    FreeDeviceMemory(d_cols);

    PetscScalar * d_val;
    MallocDeviceMemory(&d_val, nnz);
    MatSetValuesCOO(K, d_val, INSERT_VALUES);
    FreeDeviceMemory(d_val);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

    PetscPrintf(PETSC_COMM_WORLD, "Matrix K assembled with %d non-zero entries.\n", nnz);

    MatView(K, PETSC_VIEWER_STDOUT_WORLD);

    delete fm;

    PetscFinalize();
    return 0;
}