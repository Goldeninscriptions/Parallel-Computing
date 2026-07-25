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
    const int num_parts = part_num_x * part_num_y;
    const int total_owned = std::count_if(ID.begin(), ID.end(), [](int id) { return id != -1; });

    FileManager * fm = new FileManager();

    std::vector<int> num_local_funcs_x(part_num_x, 0);
    std::vector<int> num_local_funcs_y(part_num_y, 0);
    const int base_size_x = m / part_num_x;
    const int rem_size_x = m % part_num_x;
    const int base_size_y = n / part_num_y;
    const int rem_size_y = n % part_num_y;

    for (int i = 0; i < part_num_x; ++i)
    {
        num_local_funcs_x[i] = base_size_x + (i < rem_size_x ? 1 : 0);
    }

    for (int i = 0; i < part_num_y; ++i)
    {
        num_local_funcs_y[i] = base_size_y + (i < rem_size_y ? 1 : 0);
    }

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
        elem_start_idx_x[i] = std::max(0, node_start_idx_x[i] - p);
        elem_end_idx_x[i - 1] = elem_start_idx_x[i] - 1;
    }
    elem_end_idx_x[part_num_x - 1] = nElemX - 1;
    
    for (int i = 1; i < part_num_y; ++i)
    {
        elem_start_idx_y[i] = std::max(0, node_start_idx_y[i] - q);
        elem_end_idx_y[i - 1] = elem_start_idx_y[i] - 1;
    }
    elem_end_idx_y[part_num_y - 1] = nElemY - 1;

    int i,j;
    #pragma omp parallel for collapse(2) private(i, j)
    for (j = 0; j < part_num_y; ++j)
    {
        for (i = 0; i < part_num_x; ++i)
        {
            int count = j * part_num_x + i;
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
                        local_to_global_total.push_back(IEN[index * nLocBas + ii]);
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
                    const int global = (node_start_idx_y[j] + localj) * m + (node_start_idx_x[i] + locali);
                    local_to_global.push_back(global);
                }
            }

            std::vector<int> tempID{};
            std::vector<double> localCP{};
            std::vector<int> localID{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    localCP.push_back(CP[local_to_global_total[ii] * dim + jj]);
                }
                if (std::find(local_to_global.begin(), local_to_global.end(), local_to_global_total[ii]) != local_to_global.end())
                    tempID.push_back(ID[local_to_global_total[ii]]);
                else
                    tempID.push_back(-1);
                localID.push_back(ID[local_to_global_total[ii]]);
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
                            std::find(local_to_global_total.begin(), local_to_global_total.end(), IEN[index * nLocBas + ii])));
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

            const int nlocalfunc = total_owned / num_parts + (count < total_owned % num_parts ? 1 : 0);
            const int nlocalelemx = elem_end_idx_x[i] - elem_start_idx_x[i] + 1;
            const int nlocalelemy = elem_end_idx_y[j] - elem_start_idx_y[j] + 1;

            std::string filename = fm->GetPartitionFilename(base_name, count);
            fm->WritePartition(filename, nlocalfunc,
                nlocalelemx, nlocalelemy,
                elem_size1, elem_size2, localCP, localID,
                localIEN, localNURBSExtraction1, localNURBSExtraction2);

            std::cout << "Partition " << count << " generated." << std::endl;
        }
    }
}

void Partition::GeneratePartition(const int &nElemX, const int &nElemY,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID)
{
    FileManager * fm = new FileManager();

    const int m = nElemX + 1;
    const int n = nElemY + 1;
    const int nLocBas = 4;
    const int num_parts = part_num_x * part_num_y;
    const int total_owned = std::count_if(ID.begin(), ID.end(), [](int id) { return id != -1; });

    std::vector<int> num_local_funcs_x(part_num_x, 0);
    std::vector<int> num_local_funcs_y(part_num_y, 0);
    const int base_size_x = m / part_num_x;
    const int rem_size_x = m % part_num_x;
    const int base_size_y = n / part_num_y;
    const int rem_size_y = n % part_num_y;

    for (int i = 0; i < part_num_x; ++i)
    {
        num_local_funcs_x[i] = base_size_x + (i < rem_size_x ? 1 : 0);
    }

    for (int i = 0; i < part_num_y; ++i)
    {
        num_local_funcs_y[i] = base_size_y + (i < rem_size_y ? 1 : 0);
    }

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
                        local_to_global_total.push_back(IEN[index * nLocBas + ii]);
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
                    const int global = (node_start_idx_y[j] + localj) * m + (node_start_idx_x[i] + locali);
                    local_to_global.push_back(global);
                }
            }

            std::vector<int> tempID{};
            std::vector<double> localCP{};
            std::vector<int> localID{};
            for (int ii = 0; ii < local_to_global_total.size(); ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    localCP.push_back(CP[local_to_global_total[ii] * dim + jj]);
                }
                if (std::find(local_to_global.begin(), local_to_global.end(), local_to_global_total[ii]) != local_to_global.end())
                    tempID.push_back(ID[local_to_global_total[ii]]);
                else
                    tempID.push_back(-1);
                localID.push_back(ID[local_to_global_total[ii]]);
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
                            std::find(local_to_global_total.begin(), local_to_global_total.end(), IEN[index * nLocBas + ii])));
                    }
                }
            }

            const int nlocalfunc = total_owned / num_parts + (count < total_owned % num_parts ? 1 : 0);
            const int nlocalelemx = elem_end_idx_x[i] - elem_start_idx_x[i] + 1;
            const int nlocalelemy = elem_end_idx_y[j] - elem_start_idx_y[j] + 1;

            std::string filename = fm->GetPartitionFilename(base_name, count);
            fm->WritePartition(filename, nlocalfunc,
                nlocalelemx, nlocalelemy, localCP, localID, localIEN);

            std::cout << "Partition " << count << " generated." << std::endl;
        }
    }
}
