#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "FileManager.hpp"
#include <algorithm>
#include <iterator>
#include "BSplineBasis.hpp"
#include "omp.h"

class Partition
{
    public:
        Partition(const int &in_part_num_x, const int &in_part_num_y,
            const int &in_dim, const std::string &in_base_name) :
            part_num_x(in_part_num_x), part_num_y(in_part_num_y),
            dim(in_dim), base_name(in_base_name) {}
        
        ~Partition(){}

        void GeneratePartition(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
            const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID,
            const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2);
        
        void GeneratePartitionSerial(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
            const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID,
            const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2);

        void GeneratePartition(const int &nElemX, const int &nElemY,
            const std::vector<double> &CP, const std::vector<int> &ID, const std::vector<int> &IEN);

        void GeneratePartitionSerial(const int &nElemX, const int &nElemY,
            const std::vector<double> &CP, const std::vector<int> &ID, const std::vector<int> &IEN);

    private:
        const int part_num_x, part_num_y;
        const int dim;
        const std::string base_name;
};

#endif