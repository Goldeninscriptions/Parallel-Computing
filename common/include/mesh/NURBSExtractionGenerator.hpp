#ifndef NURBSEXTRACTIONGENERATOR_HPP
#define NURBSEXTRACTIONGENERATOR_HPP

// Ref : Borden et al. 2011 IJNME, 87:15-47
// GenerateExtraction1D : Algorithm 1

#include "BSplineBasis.hpp"

class NURBSExtractionGenerator
{
    public:
        std::vector<double> GenerateExtraction1D(const BSplineBasis * const &basis);

        int idxC(const int &i, const int &j, const int &k, const int &a) const
        {
            return (j-1) + (i-1) * a + (k-1) * a * a;
        }
};

#endif