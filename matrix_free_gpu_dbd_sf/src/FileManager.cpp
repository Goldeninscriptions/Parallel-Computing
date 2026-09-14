#include "FileManager.hpp"

void FileManager::WritePartition(const std::string &filename,
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
    const std::vector<double> &NURBSExtraction2) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "nlocalfunc: " << nlocalfunc << std::endl;

    file << "nlocalelemx: " << nlocalelemx << std::endl;

    file << "nlocalelemy: " << nlocalelemy << std::endl;

    file << "ID" << std::endl;
    for (int ii = 0; ii < ID.size(); ++ii) file << ID[ii] << " ";
    file << std::endl;

    file << "ghostID" << std::endl;
    for (int ii = 0; ii < ghostID.size(); ++ii) file << ghostID[ii] << " ";
    file << std::endl;

    file << "Dir" << std::endl;
    for (int ii = 0; ii < Dir.size(); ++ii) file << Dir[ii] << " ";
    file << std::endl;

    file << "IEN" << std::endl;
    for (int ii = 0; ii < IEN.size(); ++ii) file << IEN[ii] << " ";
    file << std::endl;

    file << std::setprecision(16);

    file << "ElemSize1" << std::endl;
    for (int ii = 0; ii < elem_size1.size(); ++ii) file << elem_size1[ii] << " ";
    file << std::endl;

    file << "ElemSize2" << std::endl;
    for (int ii = 0; ii < elem_size2.size(); ++ii) file << elem_size2[ii] << " ";
    file << std::endl;

    file << "CP" << std::endl;
    for (int ii = 0; ii < CP.size(); ++ii) file << CP[ii] << " ";
    file << std::endl;

    file << "NURBSExtraction1" << std::endl;
    for (int ii = 0; ii < NURBSExtraction1.size(); ++ii) file << NURBSExtraction1[ii] << " ";
    file << std::endl;

    file << "NURBSExtraction2" << std::endl;
    for (int ii = 0; ii < NURBSExtraction2.size(); ++ii) file << NURBSExtraction2[ii] << " ";
    file << std::endl;

    file.close();
}

void FileManager::WritePartition(const std::string &filename,
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
    const std::vector<double> &NURBSExtraction2) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "nlocalfunc: " << nlocalfunc << std::endl;

    file << "nlocalelemx: " << nlocalelemx << std::endl;

    file << "nlocalelemy: " << nlocalelemy << std::endl;

    file << "ID" << std::endl;
    for (int ii = 0; ii < ID.size(); ++ii) file << ID[ii] << " ";
    file << std::endl;

    file << "Dir" << std::endl;
    for (int ii = 0; ii < Dir.size(); ++ii) file << Dir[ii] << " ";
    file << std::endl;

    file << "IEN" << std::endl;
    for (int ii = 0; ii < IEN.size(); ++ii) file << IEN[ii] << " ";
    file << std::endl;

    file << std::setprecision(16);

    file << "ElemSize1" << std::endl;
    for (int ii = 0; ii < elem_size1.size(); ++ii) file << elem_size1[ii] << " ";
    file << std::endl;

    file << "ElemSize2" << std::endl;
    for (int ii = 0; ii < elem_size2.size(); ++ii) file << elem_size2[ii] << " ";
    file << std::endl;

    file << "CP" << std::endl;
    for (int ii = 0; ii < CP.size(); ++ii) file << CP[ii] << " ";
    file << std::endl;

    file << "NURBSExtraction1" << std::endl;
    for (int ii = 0; ii < NURBSExtraction1.size(); ++ii) file << NURBSExtraction1[ii] << " ";
    file << std::endl;

    file << "NURBSExtraction2" << std::endl;
    for (int ii = 0; ii < NURBSExtraction2.size(); ++ii) file << NURBSExtraction2[ii] << " ";
    file << std::endl;

    file.close();
}

void FileManager::WritePartition(const std::string &filename,
    const int &nlocalfunc,
    const int &nlocalelemx,
    const int &nlocalelemy,
    const std::vector<double> &CP,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<int> &IEN) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "nlocalfunc: " << nlocalfunc << std::endl;

    file << "nlocalelemx: " << nlocalelemx << std::endl;

    file << "nlocalelemy: " << nlocalelemy << std::endl;

    file << "ID" << std::endl;
    for (int ii = 0; ii < ID.size(); ++ii) file << ID[ii] << " ";
    file << std::endl;

    file << "Dir" << std::endl;
    for (int ii = 0; ii < Dir.size(); ++ii) file << Dir[ii] << " ";
    file << std::endl;

    file << "IEN" << std::endl;
    for (int ii = 0; ii < IEN.size(); ++ii) file << IEN[ii] << " ";
    file << std::endl;

    file << std::setprecision(16);

    file << "CP" << std::endl;
    for (int ii = 0; ii < CP.size(); ++ii) file << CP[ii] << " ";
    file << std::endl;

    file.close();
}

void FileManager::ReadPartition(const std::string &filename,
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
    std::vector<double> &NURBSExtraction2) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    std::getline(file, line);
    if (line.find("nlocalfunc: ") != std::string::npos)
    {
        nlocalfunc = std::stoi(line.substr(12));
    }
    std::getline(file, line);
    if (line.find("nlocalelemx: ") != std::string::npos)
    {
        nlocalelemx = std::stoi(line.substr(13));
    }
    std::getline(file, line);
    if (line.find("nlocalelemy: ") != std::string::npos)
    {
        nlocalelemy = std::stoi(line.substr(13));
    }
    while (std::getline(file, line))
    {
        if (line == "ElemSize1")
        {
            std::string elem_size_str;
            std::getline(file, elem_size_str);
            std::istringstream elem_size_ss(elem_size_str);
            elem_size1.clear();
            double d;
            while (elem_size_ss >> d)
            {
                elem_size1.push_back(d);
                if (elem_size_ss.peek() == ' ')
                    elem_size_ss.ignore();
            }
        }
        else if (line == "ElemSize2")
        {
            std::string elem_size_str;
            std::getline(file, elem_size_str);
            std::istringstream elem_size_ss(elem_size_str);
            elem_size2.clear();
            double d;
            while (elem_size_ss >> d)
            {
                elem_size2.push_back(d);
                if (elem_size_ss.peek() == ' ')
                    elem_size_ss.ignore();
            }
        }
        else if (line == "CP")
        {
            std::string CP_str;
            std::getline(file, CP_str);
            std::istringstream CP_ss(CP_str);
            CP.clear();
            double d;
            while (CP_ss >> d)
            {
                CP.push_back(d);
                if (CP_ss.peek() == ' ')
                    CP_ss.ignore();
            }
        }
        else if (line == "ID")
        {
            std::string ID_str;
            std::getline(file, ID_str);
            std::istringstream ID_ss(ID_str);
            ID.clear();
            int i;
            while (ID_ss >> i)
            {
                ID.push_back(i);
                if (ID_ss.peek() == ' ')
                    ID_ss.ignore();
            }
        }
        else if (line == "ghostID")
        {
            std::string ghostID_str;
            std::getline(file, ghostID_str);
            std::istringstream ghostID_ss(ghostID_str);
            ghostID.clear();
            int i;
            while (ghostID_ss >> i)
            {
                ghostID.push_back(i);
                if (ghostID_ss.peek() == ' ')
                    ghostID_ss.ignore();
            }
        }
        else if (line == "Dir")
        {
            std::string Dir_str;
            std::getline(file, Dir_str);
            std::istringstream Dir_ss(Dir_str);
            Dir.clear();
            int i;
            while (Dir_ss >> i)
            {
                Dir.push_back(i);
                if (Dir_ss.peek() == ' ')
                    Dir_ss.ignore();
            }
        }
        else if (line == "IEN")
        {
            std::string IEN_str;
            std::getline(file, IEN_str);
            std::istringstream IEN_ss(IEN_str);
            IEN.clear();
            int i;
            while (IEN_ss >> i)
            {
                IEN.push_back(i);
                if (IEN_ss.peek() == ' ')
                    IEN_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction1")
        {
            std::string NURBSExtraction1_str;
            std::getline(file, NURBSExtraction1_str);
            std::istringstream NURBSExtraction1_ss(NURBSExtraction1_str);
            NURBSExtraction1.clear();
            double d;
            while (NURBSExtraction1_ss >> d)
            {
                NURBSExtraction1.push_back(d);
                if (NURBSExtraction1_ss.peek() == ' ')
                    NURBSExtraction1_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction2")
        {
            std::string NURBSExtraction2_str;
            std::getline(file, NURBSExtraction2_str);
            std::istringstream NURBSExtraction2_ss(NURBSExtraction2_str);
            NURBSExtraction2.clear();
            double d;
            while (NURBSExtraction2_ss >> d)
            {
                NURBSExtraction2.push_back(d);
                if (NURBSExtraction2_ss.peek() == ' ')
                    NURBSExtraction2_ss.ignore();
            }
        }
    }
}

void FileManager::ReadPartition(const std::string &filename,
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
    std::vector<double> &NURBSExtraction2) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    std::getline(file, line);
    if (line.find("nlocalfunc: ") != std::string::npos)
    {
        nlocalfunc = std::stoi(line.substr(12));
    }
    std::getline(file, line);
    if (line.find("nlocalelemx: ") != std::string::npos)
    {
        nlocalelemx = std::stoi(line.substr(13));
    }
    std::getline(file, line);
    if (line.find("nlocalelemy: ") != std::string::npos)
    {
        nlocalelemy = std::stoi(line.substr(13));
    }
    while (std::getline(file, line))
    {
        if (line == "ElemSize1")
        {
            std::string elem_size_str;
            std::getline(file, elem_size_str);
            std::istringstream elem_size_ss(elem_size_str);
            elem_size1.clear();
            double d;
            while (elem_size_ss >> d)
            {
                elem_size1.push_back(d);
                if (elem_size_ss.peek() == ' ')
                    elem_size_ss.ignore();
            }
        }
        else if (line == "ElemSize2")
        {
            std::string elem_size_str;
            std::getline(file, elem_size_str);
            std::istringstream elem_size_ss(elem_size_str);
            elem_size2.clear();
            double d;
            while (elem_size_ss >> d)
            {
                elem_size2.push_back(d);
                if (elem_size_ss.peek() == ' ')
                    elem_size_ss.ignore();
            }
        }
        else if (line == "CP")
        {
            std::string CP_str;
            std::getline(file, CP_str);
            std::istringstream CP_ss(CP_str);
            CP.clear();
            double d;
            while (CP_ss >> d)
            {
                CP.push_back(d);
                if (CP_ss.peek() == ' ')
                    CP_ss.ignore();
            }
        }
        else if (line == "ID")
        {
            std::string ID_str;
            std::getline(file, ID_str);
            std::istringstream ID_ss(ID_str);
            ID.clear();
            int i;
            while (ID_ss >> i)
            {
                ID.push_back(i);
                if (ID_ss.peek() == ' ')
                    ID_ss.ignore();
            }
        }
        else if (line == "Dir")
        {
            std::string Dir_str;
            std::getline(file, Dir_str);
            std::istringstream Dir_ss(Dir_str);
            Dir.clear();
            int i;
            while (Dir_ss >> i)
            {
                Dir.push_back(i);
                if (Dir_ss.peek() == ' ')
                    Dir_ss.ignore();
            }
        }
        else if (line == "IEN")
        {
            std::string IEN_str;
            std::getline(file, IEN_str);
            std::istringstream IEN_ss(IEN_str);
            IEN.clear();
            int i;
            while (IEN_ss >> i)
            {
                IEN.push_back(i);
                if (IEN_ss.peek() == ' ')
                    IEN_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction1")
        {
            std::string NURBSExtraction1_str;
            std::getline(file, NURBSExtraction1_str);
            std::istringstream NURBSExtraction1_ss(NURBSExtraction1_str);
            NURBSExtraction1.clear();
            double d;
            while (NURBSExtraction1_ss >> d)
            {
                NURBSExtraction1.push_back(d);
                if (NURBSExtraction1_ss.peek() == ' ')
                    NURBSExtraction1_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction2")
        {
            std::string NURBSExtraction2_str;
            std::getline(file, NURBSExtraction2_str);
            std::istringstream NURBSExtraction2_ss(NURBSExtraction2_str);
            NURBSExtraction2.clear();
            double d;
            while (NURBSExtraction2_ss >> d)
            {
                NURBSExtraction2.push_back(d);
                if (NURBSExtraction2_ss.peek() == ' ')
                    NURBSExtraction2_ss.ignore();
            }
        }
    }
}

void FileManager::ReadPartition(const std::string &filename,
    int &nlocalfunc,
    int &nlocalelemx,
    int &nlocalelemy,
    std::vector<double> &CP,
    std::vector<int> &ID,
    std::vector<int> &Dir,
    std::vector<int> &IEN) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    std::getline(file, line);
    if (line.find("nlocalfunc: ") != std::string::npos)
    {
        nlocalfunc = std::stoi(line.substr(12));
    }
    std::getline(file, line);
    if (line.find("nlocalelemx: ") != std::string::npos)
    {
        nlocalelemx = std::stoi(line.substr(13));
    }
    std::getline(file, line);
    if (line.find("nlocalelemy: ") != std::string::npos)
    {
        nlocalelemy = std::stoi(line.substr(13));
    }
    while (std::getline(file, line))
    {
        if (line == "CP")
        {
            std::string CP_str;
            std::getline(file, CP_str);
            std::istringstream CP_ss(CP_str);
            CP.clear();
            double d;
            while (CP_ss >> d)
            {
                CP.push_back(d);
                if (CP_ss.peek() == ' ')
                    CP_ss.ignore();
            }
        }
        else if (line == "ID")
        {
            std::string ID_str;
            std::getline(file, ID_str);
            std::istringstream ID_ss(ID_str);
            ID.clear();
            int i;
            while (ID_ss >> i)
            {
                ID.push_back(i);
                if (ID_ss.peek() == ' ')
                    ID_ss.ignore();
            }
        }
        else if (line == "Dir")
        {
            std::string Dir_str;
            std::getline(file, Dir_str);
            std::istringstream Dir_ss(Dir_str);
            Dir.clear();
            int i;
            while (Dir_ss >> i)
            {
                Dir.push_back(i);
                if (Dir_ss.peek() == ' ')
                    Dir_ss.ignore();
            }
        }
        else if (line == "IEN")
        {
            std::string IEN_str;
            std::getline(file, IEN_str);
            std::istringstream IEN_ss(IEN_str);
            IEN.clear();
            int i;
            while (IEN_ss >> i)
            {
                IEN.push_back(i);
                if (IEN_ss.peek() == ' ')
                    IEN_ss.ignore();
            }
        }
    }
}

void FileManager::ReadPartition(const std::string &filename,
    int &nlocalfunc) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    std::getline(file, line);
    if (line.find("nlocalfunc: ") != std::string::npos)
    {
        nlocalfunc = std::stoi(line.substr(12));
    }
}

void FileManager::WritePreprocessInfo(const std::string &filename, const int &p, const int &q, const double &Lx, const double &Ly,
    const int &nElemX, const int &nElemY, const int &part_num_1d, const int &dim, const std::string &base_name) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "p: " << p << std::endl;
    file << "q: " << q << std::endl;
    file << "Lx: " << Lx << std::endl;
    file << "Ly: " << Ly << std::endl;
    file << "nElemX: " << nElemX << std::endl;
    file << "nElemY: " << nElemY << std::endl;
    file << "part_num_1d: " << part_num_1d << std::endl;
    file << "dim: " << dim << std::endl;
    file << "base_name: " << base_name << std::endl;

    file.close();
}

void FileManager::ReadPreprocessInfo(const std::string &filename, int &p, int &q, double &Lx, double &Ly,
    int &nElemX, int &nElemY, int &part_num_1d, int &dim, std::string &base_name) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("p: ") != std::string::npos)
        {
            p = std::stoi(line.substr(3));
        }
        else if (line.find("q: ") != std::string::npos)
        {
            q = std::stoi(line.substr(3));
        }
        else if (line.find("Lx: ") != std::string::npos)
        {
            Lx = std::stod(line.substr(4));
        }
        else if (line.find("Ly: ") != std::string::npos)
        {
            Ly = std::stod(line.substr(4));
        }
        else if (line.find("nElemX: ") != std::string::npos)
        {
            nElemX = std::stoi(line.substr(8));
        }
        else if (line.find("nElemY: ") != std::string::npos)
        {
            nElemY = std::stoi(line.substr(8));
        }
        else if (line.find("part_num_1d: ") != std::string::npos)
        {
            part_num_1d = std::stoi(line.substr(13));
        }
        else if (line.find("dim: ") != std::string::npos)
        {
            dim = std::stoi(line.substr(5));
        }
        else if (line.find("base_name: ") != std::string::npos)
        {
            base_name = line.substr(11);
        }
    }
}

void FileManager::WriteNonZeroCoordinate(const std::string &filename, const int &nnz,
    const std::vector<int> &rows, const std::vector<int> &cols) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "nnz: " << nnz << std::endl;

    file << "rows: ";
    for (const auto &row : rows)
    {
        file << row << " ";
    }
    file << std::endl;

    file << "cols: ";
    for (const auto &col : cols)
    {
        file << col << " ";
    }
    file << std::endl;

    file.close();
}

void FileManager::ReadNonZeroCoordinate(const std::string &filename, int &nnz,
    std::vector<int> &rows, std::vector<int> &cols) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("nnz: ") != std::string::npos)
        {
            nnz = std::stoi(line.substr(5));
        }
        else if (line.find("rows: ") != std::string::npos)
        {
            std::istringstream iss(line.substr(6));
            int row;
            while (iss >> row)
            {
                rows.push_back(row);
            }
        }
        else if (line.find("cols: ") != std::string::npos)
        {
            std::istringstream iss(line.substr(6));
            int col;
            while (iss >> col)
            {
                cols.push_back(col);
            }
        }
    }
}

void FileManager::WriteNewToOldMapping(const std::string &filename, const std::vector<int> &new_to_old) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    for (const auto &value : new_to_old)
    {
        file << value << " ";
    }
    file << std::endl;

    file.close();
}

void FileManager::ReadNewToOldMapping(const std::string &filename, std::vector<int> &new_to_old) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    int value;
    while (iss >> value)
    {
        new_to_old.push_back(value);
    }
    file.close();
}

std::string FileManager::GetPartitionFilename(const std::string &base_name, const int &rank) const
{
    return base_name + "_" + std::to_string(rank) + ".txt";
}

std::string FileManager::GetNonZeroCoordinateFilename(const std::string &base_name, const int &rank) const
{
    return base_name + "_nonzero_" + std::to_string(rank) + ".txt";
}