#include <petscmat.h>
#include "FileManager.hpp"
#include "memory.cuh"
#include "GlobalAssembly.cuh"
#include "Elem2COOGenerator.hpp"
#include "AbscissaeGenerator.hpp"
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

    double hx = Lx / (nElemX + 2);
    double hy = Ly / (nElemY + 2);
    
    int nFuncX = p + nElemX;
    int nFuncY = q + nElemY;

    std::vector<double> S;
    std::vector<double> T;

    for (int i = 0; i < p - 1; ++i) S.push_back(0.0);
    for (int i = 1; i < nElemX+2; ++i) S.push_back(i * hx);
    for (int i = 0; i < p - 1; ++i) S.push_back(Lx);

    for (int i = 0; i < q - 1; ++i) T.push_back(0.0);
    for (int i = 1; i < nElemY+2; ++i) T.push_back(i * hy);
    for (int i = 0; i < q - 1; ++i) T.push_back(Ly);

    AbscissaeGenerator * absgen = new AbscissaeGenerator();
    std::vector<double> CP = absgen->GenerateAbscissae2D(S, T, p-2, q-2);

    nElemX = nFuncX - 1;
    nElemY = nFuncY - 1;

    IENGenerator * ien = new IENGenerator();
    std::vector<int> IEN = ien->GenerateIEN2D(nElemX, nElemY);

    IDGenerator * idgen = new IDGenerator();
    std::vector<int> ID = idgen->GenerateID2D(nFuncX, nFuncY);

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

    std::vector<int> old_rows = rows;
    std::vector<int> old_cols = cols;

    for (int i = 0; i < rows.size(); ++i)
    {
        old_rows[i] = new_to_old[rows[i]];
        old_cols[i] = new_to_old[cols[i]];
    }

    std::vector<std::pair<int, int>> row_col_pairs;
    row_col_pairs.reserve(old_rows.size());
    for (size_t i = 0; i < old_rows.size(); ++i) {
        row_col_pairs.emplace_back(old_rows[i], old_cols[i]);
    }
    std::sort(row_col_pairs.begin(), row_col_pairs.end());
    for (size_t i = 0; i < row_col_pairs.size(); ++i) {
        old_rows[i] = row_col_pairs[i].first;
        old_cols[i] = row_col_pairs[i].second;
    }

    int temp_row_index = 0;
    std::cout << "Row 0: ";
    for (int j = 0; j < nnz; ++j)
    {
        if (old_rows[j] != temp_row_index)
        {
            std::cout << std::endl;
            std::cout << "Row " << old_rows[j] << ": " << old_cols[j] << " ";
            temp_row_index = old_rows[j];
        }
        else
        {
            std::cout << old_cols[j] << " ";
        }
    }

    delete fm;
    delete absgen;
    delete ien;
    delete idgen;

    PetscFinalize();
    return 0;
}