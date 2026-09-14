#ifndef INVLM_HPP
#define INVLM_HPP

#include <vector>

class InvLM
{
public:
    InvLM(const int &nLocBas,
    const int &nLocElem,
    const int &nlocfunc,
    const std::vector<int> IEN);

    int GetElemNum(const int &i) const
    {
        return elemNum[i];
    }

    int GetOffset(const int &i) const
    {
        return offset[i];
    }

    std::vector<int> GetElemIdx(const int &i) const
    {
        std::vector<int> result(elemNum[i], -1);
        for (int j = 0; j < elemNum[i]; ++j)
        {
            result[j] = elemIdx[offset[i] + j];
        }
        return result;
    }

    std::vector<int> GetBaseIdx(const int &i) const
    {
        std::vector<int> result(elemNum[i], -1);
        for (int j = 0; j < elemNum[i]; ++j)
        {
            result[j] = baseIdx[offset[i] + j];
        }
        return result;
    }

    std::vector<int> GetAllElemIdx() const
    {
        return elemIdx;
    }

    std::vector<int> GetAllBaseIdx() const
    {
        return baseIdx;
    }

    std::vector<int> GetAllElemNum() const
    {
        return elemNum;
    }

    std::vector<int> GetAllOffset() const
    {
        return offset;
    }

private:
    std::vector<int> elemNum{};
    std::vector<int> offset{};
    std::vector<int> elemIdx{};
    std::vector<int> baseIdx{};
};

#endif