#include "Partition.hpp"

void Partition::GeneratePartition(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID,
    const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;
    const int nElemX = m - p;
    const int nElemY = n - q;
    const int nLocBas = (p + 1) * (q + 1);
    const int extSizeX = (p + 1) * (p + 1);
    const int extSizeY = (q + 1) * (q + 1);
    const double hx = (S.back() - S.front()) / nElemX;
    const double hy = (T.back() - T.front()) / nElemY;

    FileManager * fm = new FileManager();

    const int part_size_x = (m % part_num_x == 0) ? m / part_num_x : m / part_num_x + 1;
    const int part_size_y = (n % part_num_y == 0) ? n / part_num_y : n / part_num_y + 1;
    
    std::vector<int> num_local_funcs_x(part_num_x, 0);
    std::vector<int> num_local_funcs_y(part_num_y, 0);

    for (int i = 0; i < m; i += part_size_x) num_local_funcs_x[i/part_size_x] = part_size_x;
    num_local_funcs_x[part_num_x - 1] = m - (part_num_x - 1) * part_size_x;

    for (int i = 0; i < n; i += part_size_y) num_local_funcs_y[i/part_size_y] = part_size_y;
    num_local_funcs_y[part_num_y - 1] = n - (part_num_y - 1) * part_size_y;

    std::vector<int> node_start_idx_x(part_num_x, 0);
    std::vector<int> node_start_idx_y(part_num_y, 0);
    std::vector<int> node_end_idx_x(part_num_x, 0);
    std::vector<int> node_end_idx_y(part_num_y, 0);

    for (int i = 1; i < part_num_x; ++i)
    {
        node_start_idx_x[i] = node_start_idx_x[i - 1] + num_local_funcs_x[i - 1];
        node_end_idx_x[i - 1] = node_start_idx_x[i] - 1;
    }
    node_end_idx_x[part_num_x - 1] = m - 1;

    for (int i = 1; i < part_num_y; ++i)
    {
        node_start_idx_y[i] = node_start_idx_y[i - 1] + num_local_funcs_y[i - 1];
        node_end_idx_y[i - 1] = node_start_idx_y[i] - 1;
    }
    node_end_idx_y[part_num_y - 1] = n - 1;

    std::vector<int> elem_start_idx_x(part_num_x, 0);
    std::vector<int> elem_start_idx_y(part_num_y, 0);
    std::vector<int> elem_end_idx_x(part_num_x, 0);
    std::vector<int> elem_end_idx_y(part_num_y, 0);

    for (int i = 1; i < part_num_x; ++i)
    {
        elem_start_idx_x[i] = node_start_idx_x[i] - p/2;
        elem_end_idx_x[i - 1] = elem_start_idx_x[i] - 1;
    }
    elem_end_idx_x[part_num_x - 1] = nElemX - 1;
    
    for (int i = 1; i < part_num_y; ++i)
    {
        elem_start_idx_y[i] = node_start_idx_y[i] - q/2;
        elem_end_idx_y[i - 1] = elem_start_idx_y[i] - 1;
    }
    elem_end_idx_y[part_num_y - 1] = nElemY - 1;

    std::vector<int> old_to_new(m*n, 0);
    int old_to_new_count = 0;
    for (int jj = 0; jj < part_num_y; ++jj)
    {
        for (int ii = 0; ii < part_num_x; ++ii)
        {
            for (int localj = 0; localj < num_local_funcs_y[jj]; ++localj)
            {
                for (int locali = 0; locali < num_local_funcs_x[ii]; ++locali)
                {
                    const int globali = ii * part_size_x + locali;
                    const int globalj = jj * part_size_y + localj;
                    const int global = globalj * m + globali;
                    old_to_new[global] = old_to_new_count;
                    old_to_new_count++;
                }
            }
        }
    }

    std::vector<double> newCP = CP;
    for (int ii = 0; ii < CP.size() / dim; ++ii)
    {
        for (int jj = 0; jj <  dim; ++jj)
        {
            newCP[old_to_new[ii] * dim + jj] = CP[ii * dim + jj];
        }
    }

    std::vector<int> newID = ID;
    for (int ii = 0; ii < ID.size(); ++ii)
    {
        if (ID[ii] == -1)
            newID[old_to_new[ii]] = -1;
        else
            newID[old_to_new[ii]] = old_to_new[ii];
    }

    std::vector<int> newIEN = IEN;
    for (int ii = 0; ii < IEN.size(); ++ii)
    {
        newIEN[ii] = old_to_new[IEN[ii]];
    }

    int i,j;
    #pragma omp parallel for collapse(2) private(i, j)
    for (j = 0; j < part_num_y; ++j)
    {
        for (i = 0; i < part_num_x; ++i)
        {
            int count = j * part_num_y + i;
            std::cout << "Generating partition " << count << "..." << std::endl;

            std::vector<double> elem_size1(elem_end_idx_x[i] - elem_start_idx_x[i] + 1, hx);
            std::vector<double> elem_size2(elem_end_idx_y[j] - elem_start_idx_y[j] + 1, hy);

            std::vector<int> local_to_global_total{};
            for (int localj = elem_start_idx_y[j]; localj <= elem_end_idx_y[j]; ++localj)
            {
                for (int locali = elem_start_idx_x[i]; locali <= elem_end_idx_x[i]; ++locali)
                {
                    int index = localj * nElemX + locali;
                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        local_to_global_total.push_back(newIEN[index * nLocBas + ii]);
                    }
                }
            }

            std::sort(local_to_global_total.begin(), local_to_global_total.end());
            local_to_global_total.erase(std::unique(local_to_global_total.begin(), local_to_global_total.end()), local_to_global_total.end());

            std::vector<int> local_to_global{};
            for (int localj = 0; localj < num_local_funcs_y[j]; ++localj)
            {
                for (int locali = 0; locali < num_local_funcs_x[i]; ++locali)
                {
                    const int globali = i * part_size_x + locali;
                    const int globalj = j * part_size_y + localj;
                    const int global = globalj * m + globali;
                    local_to_global.push_back(old_to_new[global]);
                }
            }

            std::vector<int> local_to_global_ghost{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                if (std::find(local_to_global.begin(), local_to_global.end(), local_to_global_total[ii]) == local_to_global.end())
                {
                    local_to_global_ghost.push_back(local_to_global_total[ii]);
                }
            }

            local_to_global_total.clear();
            local_to_global_total.insert(local_to_global_total.end(), local_to_global.begin(), local_to_global.end());
            local_to_global_total.insert(local_to_global_total.end(), local_to_global_ghost.begin(), local_to_global_ghost.end());

            std::vector<double> localCP{};
            std::vector<int> localID{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    localCP.push_back(newCP[local_to_global_total[ii] * dim + jj]);
                }
                localID.push_back(newID[local_to_global_total[ii]]);
            }

            std::vector<int> localDir{};
            for (int ii = 0; ii < local_to_global.size(); ++ii)
            {
                if (localID[ii] == -1)
                    localDir.push_back(local_to_global[ii]);
            }

            std::vector<int> localIEN{};
            for (int localj = elem_start_idx_y[j]; localj <= elem_end_idx_y[j]; ++localj)
            {
                for (int locali = elem_start_idx_x[i]; locali <= elem_end_idx_x[i]; ++locali)
                {
                    int index = localj * nElemX + locali;
                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        localIEN.push_back(std::distance(local_to_global_total.begin(),
                            std::find(local_to_global_total.begin(), local_to_global_total.end(), newIEN[index * nLocBas + ii])));
                    }
                }
            }

            std::vector<double> localNURBSExtraction1{};
            for (int ii = elem_start_idx_x[i]; ii <= elem_end_idx_x[i]; ++ii)
            {
                size_t start = ii * extSizeX;
                size_t end = (ii + 1) * extSizeX;
                    
                std::copy(NURBSExtraction1.begin() + start, NURBSExtraction1.begin() + end,
                    std::back_inserter(localNURBSExtraction1));
            }

            std::vector<double> localNURBSExtraction2{};
            for (int jj = elem_start_idx_y[j]; jj <= elem_end_idx_y[j]; ++jj)
            {
                size_t start = jj * extSizeY;
                size_t end = (jj + 1) * extSizeY;
                    
                std::copy(NURBSExtraction2.begin() + start, NURBSExtraction2.begin() + end,
                    std::back_inserter(localNURBSExtraction2));
            }

            const int nlocalfunc = local_to_global.size();
            const int nlocalelemx = elem_end_idx_x[i] - elem_start_idx_x[i] + 1;
            const int nlocalelemy = elem_end_idx_y[j] - elem_start_idx_y[j] + 1;

            std::string filename = fm->GetPartitionFilename(base_name, count);
            fm->WritePartition(filename, nlocalfunc,
                nlocalelemx, nlocalelemy,
                elem_size1, elem_size2, localCP, localID, local_to_global_ghost,
                localDir, localIEN, localNURBSExtraction1, localNURBSExtraction2);

            std::cout << "Partition " << count << " generated." << std::endl;
        }
    }
    delete fm;
}

void Partition::GeneratePartitionSerial(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID,
    const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;
    const int nElemX = m - p;
    const int nElemY = n - q;
    const double hx = (S.back() - S.front()) / nElemX;
    const double hy = (T.back() - T.front()) / nElemY;

    std::vector<double> elem_size1(nElemX, hx);
    std::vector<double> elem_size2(nElemY, hy);

    FileManager * fm = new FileManager();

    std::vector<int> newID = ID;
    std::vector<int> Dir{};
    for (int ii = 0; ii < ID.size(); ++ii)
    {
        if (ID[ii] == -1)
        {
            newID[ii] = -1;
            Dir.push_back(ii);
        }
        else
            newID[ii] = ii;
    }

    int count = 0;
    std::cout << "Generating partition " << count << "..." << std::endl;

    std::string filename = fm->GetPartitionFilename(base_name, 0);
    fm->WritePartition(filename, m * n, nElemX, nElemY,
                elem_size1, elem_size2, CP, newID,
                Dir, IEN, NURBSExtraction1, NURBSExtraction2);

    std::cout << "Partition " << count << " generated." << std::endl;
    delete fm;
}

void Partition::GeneratePartition(const int &nElemX, const int &nElemY,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID)
{
    FileManager * fm = new FileManager();

    const int m = nElemX + 1;
    const int n = nElemY + 1;
    const int nLocBas = 4;

    const int part_size_x = (m % part_num_x == 0) ? m / part_num_x : m / part_num_x + 1;
    const int part_size_y = (n % part_num_y == 0) ? n / part_num_y : n / part_num_y + 1;
    
    std::vector<int> num_local_funcs_x(part_num_x, 0);
    std::vector<int> num_local_funcs_y(part_num_y, 0);

    for (int i = 0; i < m; i += part_size_x) num_local_funcs_x[i/part_size_x] = part_size_x;
    num_local_funcs_x[part_num_x - 1] = m - (part_num_x - 1) * part_size_x;

    for (int i = 0; i < n; i += part_size_y) num_local_funcs_y[i/part_size_y] = part_size_y;
    num_local_funcs_y[part_num_y - 1] = n - (part_num_y - 1) * part_size_y;

    std::vector<int> node_start_idx_x(part_num_x, 0);
    std::vector<int> node_start_idx_y(part_num_y, 0);
    std::vector<int> node_end_idx_x(part_num_x, 0);
    std::vector<int> node_end_idx_y(part_num_y, 0);

    for (int i = 1; i < part_num_x; ++i)
    {
        node_start_idx_x[i] = node_start_idx_x[i - 1] + num_local_funcs_x[i - 1];
        node_end_idx_x[i - 1] = node_start_idx_x[i] - 1;
    }
    node_end_idx_x[part_num_x - 1] = m - 1;

    for (int i = 1; i < part_num_y; ++i)
    {
        node_start_idx_y[i] = node_start_idx_y[i - 1] + num_local_funcs_y[i - 1];
        node_end_idx_y[i - 1] = node_start_idx_y[i] - 1;
    }
    node_end_idx_y[part_num_y - 1] = n - 1;

    std::vector<int> elem_start_idx_x(part_num_x, 0);
    std::vector<int> elem_start_idx_y(part_num_y, 0);
    std::vector<int> elem_end_idx_x(part_num_x, 0);
    std::vector<int> elem_end_idx_y(part_num_y, 0);

    for (int i = 1; i < part_num_x; ++i)
    {
        elem_start_idx_x[i] = node_start_idx_x[i];
        elem_end_idx_x[i - 1] = elem_start_idx_x[i] - 1;
    }
    elem_end_idx_x[part_num_x - 1] = nElemX - 1;
    
    for (int i = 1; i < part_num_y; ++i)
    {
        elem_start_idx_y[i] = node_start_idx_y[i];
        elem_end_idx_y[i - 1] = elem_start_idx_y[i] - 1;
    }
    elem_end_idx_y[part_num_y - 1] = nElemY - 1;

    std::vector<int> old_to_new(m*n, 0);
    std::vector<int> new_to_old(m*n, 0);
    int old_to_new_count = 0;
    for (int jj = 0; jj < part_num_y; ++jj)
    {
        for (int ii = 0; ii < part_num_x; ++ii)
        {
            for (int localj = 0; localj < num_local_funcs_y[jj]; ++localj)
            {
                for (int locali = 0; locali < num_local_funcs_x[ii]; ++locali)
                {
                    const int globali = ii * part_size_x + locali;
                    const int globalj = jj * part_size_y + localj;
                    const int global = globalj * m + globali;
                    old_to_new[global] = old_to_new_count;
                    new_to_old[old_to_new_count] = global;
                    old_to_new_count++;
                }
            }
        }
    }

    std::string map_name = "new_to_old_mapping.txt";
    fm->WriteNewToOldMapping(map_name, new_to_old);

    std::vector<double> newCP = CP;
    for (int ii = 0; ii < CP.size() / dim; ++ii)
    {
        for (int jj = 0; jj <  dim; ++jj)
        {
            newCP[old_to_new[ii] * dim + jj] = CP[ii * dim + jj];
        }
    }

    std::vector<int> newID = ID;
    for (int ii = 0; ii < ID.size(); ++ii)
    {
        if (ID[ii] == -1)
            newID[old_to_new[ii]] = -1;
        else
            newID[old_to_new[ii]] = old_to_new[ii];
    }

    std::vector<int> newIEN = IEN;
    for (int ii = 0; ii < IEN.size(); ++ii)
    {
        newIEN[ii] = old_to_new[IEN[ii]];
    }

    int i,j;
    #pragma omp parallel for collapse(2) private(i, j)
    for (j = 0; j < part_num_y; ++j)
    {
        for (i = 0; i < part_num_x; ++i)
        {
            int count = j * part_num_x + i;
            std::cout << "Generating partition " << count << "..." << std::endl;

            std::vector<int> local_to_global_total{};
            for (int localj = elem_start_idx_y[j]; localj <= elem_end_idx_y[j]; ++localj)
            {
                for (int locali = elem_start_idx_x[i]; locali <= elem_end_idx_x[i]; ++locali)
                {
                    int index = localj * nElemX + locali;
                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        local_to_global_total.push_back(newIEN[index * nLocBas + ii]);
                    }
                }
            }

            std::sort(local_to_global_total.begin(), local_to_global_total.end());
            local_to_global_total.erase(std::unique(local_to_global_total.begin(), local_to_global_total.end()), local_to_global_total.end());

            std::vector<int> local_to_global{};
            for (int localj = 0; localj < num_local_funcs_y[j]; ++localj)
            {
                for (int locali = 0; locali < num_local_funcs_x[i]; ++locali)
                {
                    const int globali = i * part_size_x + locali;
                    const int globalj = j * part_size_y + localj;
                    const int global = globalj * m + globali;
                    local_to_global.push_back(old_to_new[global]);
                }
            }

            std::vector<int> local_to_global_ghost{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                if (std::find(local_to_global.begin(), local_to_global.end(), local_to_global_total[ii]) == local_to_global.end())
                {
                    local_to_global_ghost.push_back(local_to_global_total[ii]);
                }
            }

            local_to_global_total.clear();
            local_to_global_total.insert(local_to_global_total.end(), local_to_global.begin(), local_to_global.end());
            local_to_global_total.insert(local_to_global_total.end(), local_to_global_ghost.begin(), local_to_global_ghost.end());

            std::vector<double> localCP{};
            std::vector<int> localID{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    localCP.push_back(newCP[local_to_global_total[ii] * dim + jj]);
                }
                localID.push_back(newID[local_to_global_total[ii]]);
            }

            std::vector<int> localDir{};
            for (int ii = 0; ii < local_to_global.size(); ++ii)
            {
                if (localID[ii] == -1)
                    localDir.push_back(local_to_global[ii]);
            }

            std::vector<int> localIEN{};
            for (int localj = elem_start_idx_y[j]; localj <= elem_end_idx_y[j]; ++localj)
            {
                for (int locali = elem_start_idx_x[i]; locali <= elem_end_idx_x[i]; ++locali)
                {
                    int index = localj * nElemX + locali;
                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        localIEN.push_back(std::distance(local_to_global_total.begin(),
                            std::find(local_to_global_total.begin(), local_to_global_total.end(), newIEN[index * nLocBas + ii])));
                    }
                }
            }

            const int nlocalfunc = local_to_global.size();
            const int nlocalelemx = elem_end_idx_x[i] - elem_start_idx_x[i] + 1;
            const int nlocalelemy = elem_end_idx_y[j] - elem_start_idx_y[j] + 1;

            std::string filename = fm->GetPartitionFilename(base_name, count);
            fm->WritePartition(filename, nlocalfunc,
                nlocalelemx, nlocalelemy, localCP, localID, localDir, localIEN);

            std::cout << "Partition " << count << " generated." << std::endl;
        }
    }
    delete fm;
}

void Partition::GeneratePartitionSerial(const int &nElemX, const int &nElemY,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID)
{
    FileManager * fm = new FileManager();

    const int m = nElemX + 1;
    const int n = nElemY + 1;
    const int nLocBas = 4;

    std::vector<int> newID = ID;
    std::vector<int> Dir{};
    for (int ii = 0; ii < ID.size(); ++ii)
    {
        if (ID[ii] == -1)
        {
            newID[ii] = -1;
            Dir.push_back(ii);
        }
        else
            newID[ii] = ii;
    }

    int count = 0;
    std::cout << "Generating partition " << count << "..." << std::endl;

    std::string filename = fm->GetPartitionFilename(base_name, 0);

    fm->WritePartition(filename, m * n, nElemX, nElemY, CP, newID, Dir, IEN);
    std::cout << "Partition " << count << " generated." << std::endl;
    delete fm;
}

