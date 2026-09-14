#ifndef FILEMANAGER_HPP
#define FILEMANAGER_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class FileManager
{
public:
    void WritePartition(const std::string &filename,
        const int &nlocalfunc,
        const int &nlocalelemx,
        const int &nlocalelemy,
        const std::vector<double> &elem_size1,
        const std::vector<double> &elem_size2,
        const std::vector<double> &CP,
        const std::vector<int> &ID,
        const std::vector<int> &ghostID,
        const std::vector<int> &Dir,
        const std::vector<int> &IEN,
        const std::vector<double> &NURBSExtraction1,
        const std::vector<double> &NURBSExtraction2) const;

    void WritePartition(const std::string &filename,
        const int &nlocalfunc,
        const int &nlocalelemx,
        const int &nlocalelemy,
        const std::vector<double> &elem_size1,
        const std::vector<double> &elem_size2,
        const std::vector<double> &CP,
        const std::vector<int> &ID,
        const std::vector<int> &Dir,
        const std::vector<int> &IEN,
        const std::vector<double> &NURBSExtraction1,
        const std::vector<double> &NURBSExtraction2) const;
    
    void WritePartition(const std::string &filename,
        const int &nlocalfunc,
        const int &nlocalelemx,
        const int &nlocalelemy,
        const std::vector<double> &CP,
        const std::vector<int> &ID,
        const std::vector<int> &Dir,
        const std::vector<int> &IEN) const;
    
    void ReadPartition(const std::string &filename,
        int &nlocalfunc,
        int &nlocalelemx,
        int &nlocalelemy,
        std::vector<double> &elem_size1,
        std::vector<double> &elem_size2,
        std::vector<double> &CP,
        std::vector<int> &ID,
        std::vector<int> &ghostID,
        std::vector<int> &Dir,
        std::vector<int> &IEN,
        std::vector<double> &NURBSExtraction1,
        std::vector<double> &NURBSExtraction2) const;

    void ReadPartition(const std::string &filename,
        int &nlocalfunc,
        int &nlocalelemx,
        int &nlocalelemy,
        std::vector<double> &elem_size1,
        std::vector<double> &elem_size2,
        std::vector<double> &CP,
        std::vector<int> &ID,
        std::vector<int> &Dir,
        std::vector<int> &IEN,
        std::vector<double> &NURBSExtraction1,
        std::vector<double> &NURBSExtraction2) const;
    
    void ReadPartition(const std::string &filename,
        int &nlocalfunc,
        int &nlocalelemx,
        int &nlocalelemy,
        std::vector<double> &CP,
        std::vector<int> &ID,
        std::vector<int> &Dir,
        std::vector<int> &IEN) const;

    void ReadPartition(const std::string &filename,
        int &nlocalfunc) const;

    void WritePreprocessInfo(const std::string &filename,
        const int &p, const int &q,
        const double &Lx, const double &Ly,
        const int &nElemX, const int &nElemY,
        const int &part_num_1d, const int &dim,
        const std::string &base_name) const;

    void ReadPreprocessInfo(const std::string &filename,
        int &p, int &q, double &Lx, double &Ly,
        int &nElemX, int &nElemY,
        int &part_num_1d, int &dim,
        std::string &base_name) const;
    
    void WriteNonZeroCoordinate(const std::string &filename, const int &nnz,
        const std::vector<int> &rows, const std::vector<int> &cols) const;

    void ReadNonZeroCoordinate(const std::string &filename, int &nnz,
        std::vector<int> &rows, std::vector<int> &cols) const;

    void WriteNewToOldMapping(const std::string &filename, const std::vector<int> &new_to_old) const;

    void ReadNewToOldMapping(const std::string &filename, std::vector<int> &new_to_old) const;

    std::string GetPartitionFilename(const std::string &base_name, const int &rank) const;

    std::string GetNonZeroCoordinateFilename(const std::string &base_name, const int &rank) const;
};

#endif