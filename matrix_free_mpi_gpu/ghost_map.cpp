#include "FileManager.hpp"
#include <algorithm>

int main(int argc, char **argv)
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    if (part_num_1d > 1)
    {
        std::cerr << "Error: part_num_1d should be greater than 1 to generate ghost map." << std::endl;
        exit(1);
    }

    std::string ghost_map_filebasename = "ghost_map";

    std::vector<int> ownlocalID;
    std::vector<int> ownghostID;
    std::vector<int> otherlocalID;
    std::vector<int> otherghostID;

    const int total_ranks = part_num_1d * part_num_1d;
    for (int rank = 0; rank < total_ranks; ++rank)
    {
        int roffset = 0;
        int soffset = 0;
        std::vector<int> recieveIndices;
        std::vector<int> recieveRanks;
        std::vector<int> recieveOffsets;
        std::vector<int> sendIndices;
        std::vector<int> sendRanks;
        std::vector<int> sendOffsets;
        std::string filename = fm->GetPartitionFilename(base_name, rank);
        fm->ReadPartition(filename, ownlocalID, ownghostID);
        for (int other_rank = 0; other_rank < total_ranks; ++other_rank)
        {
            if (other_rank != rank)
            {
                std::string other_filename = fm->GetPartitionFilename(base_name, other_rank);
                fm->ReadPartition(other_filename, otherlocalID, otherghostID);
                for (const auto &id : ownlocalID)
                {
                    if (std::find(otherghostID.begin(), otherghostID.end(), id) != otherghostID.end())
                    {
                        sendIndices.push_back(id);
                        sendRanks.push_back(other_rank);
                        soffset++;
                    }
                }
                for (const auto &id : ownghostID)
                {
                    if (std::find(otherlocalID.begin(), otherlocalID.end(), id) != otherlocalID.end())
                    {
                        recieveIndices.push_back(id);
                        recieveRanks.push_back(other_rank);
                        roffset++;
                    }
                }
                recieveOffsets.push_back(roffset);
                sendOffsets.push_back(soffset);
            }
        }
        recieveRanks.erase(std::unique(recieveRanks.begin(), recieveRanks.end()), recieveRanks.end());
        recieveOffsets.erase(std::unique(recieveOffsets.begin(), recieveOffsets.end()), recieveOffsets.end());
        sendRanks.erase(std::unique(sendRanks.begin(), sendRanks.end()), sendRanks.end());
        sendOffsets.erase(std::unique(sendOffsets.begin(), sendOffsets.end()), sendOffsets.end());

        std::string ghost_map_filename = ghost_map_filebasename + "_" + std::to_string(rank) + ".txt";
        fm->WriteGhostMap(ghost_map_filename, recieveRanks, recieveOffsets,
            recieveIndices, sendRanks, sendOffsets, sendIndices);
    }
}
