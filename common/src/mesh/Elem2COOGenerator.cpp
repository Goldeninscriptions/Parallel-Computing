#include "Elem2COOGenerator.hpp"

void Elem2COOGenerator::GenerateElem2COO(const std::vector<int> &IEN,
    const std::vector<int> &ID, const std::vector<int> &rows,
    const std::vector<int> &cols, std::vector<int> &elem2coo)
{
    elem2coo.clear();
    elem2coo.resize(IEN.size() * nLocBas);

    std::vector<int> row_counts(ID.size(), 0);
    std::vector<int> row_offsets(ID.size(), 0);

    for (size_t i = 0; i < rows.size(); ++i) {
        row_counts[rows[i]]++;
    }
    for (size_t i = 1; i < row_counts.size(); ++i) {
        row_offsets[i] = row_offsets[i - 1] + row_counts[i - 1];
    }

    #pragma omp parallel for
    for (int elem = 0; elem < nlocalelemx * nlocalelemy; ++elem)
    {
        for (int i = 0; i < nLocBas; ++i)
        {
            int II = ID[IEN[elem * nLocBas + i]];
            int base_idx = elem * nLocBas * nLocBas + i * nLocBas;
            if(II >= 0)
            {
                int Index = row_offsets[II];
                int counts = row_counts[II];
                for (int j = 0; j < nLocBas; ++j)
                {
                    int JJ = ID[IEN[elem * nLocBas + j]];
                    if(JJ >= 0)
                    {
                        auto it2 = std::find(cols.begin() + Index, cols.begin() + Index + counts, JJ);
                        int colIndex = std::distance(cols.begin() + Index, it2);
                        elem2coo[base_idx + j] = Index + colIndex;
                    }
                    else
                    {
                        elem2coo[base_idx + j] = -1; // If JJ is negative, we do not add it to the COO format
                    }
                }
            }
            else
            {
                for (int j = 0; j < nLocBas; ++j)
                {
                    elem2coo[base_idx + j] = -1; // If II is negative, we do not add it to the COO format
                }
            }
        }
    }
}

void Elem2COOGenerator::GenerateDir2COO(const std::vector<int> &ID,
    const std::vector<int> &rows, std::vector<int> &dir2coo)
{
    dir2coo.clear();
    for (size_t i = 0; i < ID.size(); ++i)
    {
        if(ID[i] == -1)
        {
            auto it = std::find(rows.begin(), rows.end(), i);
            dir2coo.push_back(std::distance(rows.begin(), it));
        }
    }
}