#include "InvLM.hpp"

InvLM::InvLM(const int &nLocBas,
    const int &nLocElem,
    const int &nlocfunc,
    const std::vector<int> IEN)
{
    elemNum.resize(nlocfunc, 0);
    offset.resize(nlocfunc, 0);
    std::vector<int> size_count(nlocfunc, 0);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = IEN[i*nLocBas + j];
            elemNum[global_func] += 1;
        }
    }

    int size = 0;
    for (int i = 0; i < nlocfunc; ++i)
    {
        offset[i] = size;
        size += elemNum[i];
    }

    elemIdx.resize(size, -1);
    baseIdx.resize(size, -1);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = IEN[i*nLocBas + j];
            if ( size_count[global_func] == elemNum[global_func] )
            {
                break;
            }
            elemIdx[size_count[global_func] + offset[global_func]] = i;
            baseIdx[size_count[global_func] + offset[global_func]] = j;
            size_count[global_func] += 1;
        }
    }
}