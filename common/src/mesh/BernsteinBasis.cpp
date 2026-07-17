#include "BernsteinBasis.hpp"

std::vector<double> BernsteinBasis::GetBernsteinBasisSingleQP(const double &xi) const
{
    std::vector<double> B(p + 1, 0.0);
    switch (p) {
        case 1:
            B[0] = 1 - xi;
            B[1] = xi;
            break;
        case 2:
            B[0] = (1 - xi) * (1 - xi);
            B[1] = 2 * xi * (1 - xi);
            B[2] = xi * xi;
            break;
        case 3:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 3 * xi * (1 - xi) * (1 - xi);
            B[2] = 3 * xi * xi * (1 - xi);
            B[3] = xi * xi * xi;
            break;
        case 4:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 4 * xi * (1 - xi) * (1 - xi) * (1 - xi);
            B[2] = 6 * xi * xi * (1 - xi) * (1 - xi);
            B[3] = 4 * xi * xi * xi * (1 - xi);
            B[4] = xi * xi * xi * xi;
            break;
        case 5:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 5 * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[2] = 10 * xi * xi * (1 - xi) * (1 - xi) * (1 - xi);
            B[3] = 10 * xi * xi * xi * (1 - xi) * (1 - xi);
            B[4] = 5 * xi * xi * xi * xi * (1 - xi);
            B[5] = xi * xi * xi * xi * xi;
            break;
        case 6:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 6 * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[2] = 15 * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[3] = 20 * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi);
            B[4] = 15 * xi * xi * xi * xi * (1 - xi) * (1 - xi);
            B[5] = 6 * xi * xi * xi * xi * xi * (1 - xi);
            B[6] = xi * xi * xi * xi * xi * xi;
            break;
        case 7:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 7 * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[2] = 21 * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[3] = 35 * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[4] = 35 * xi * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi);
            B[5] = 21 * xi * xi * xi * xi * xi * (1 - xi) * (1 - xi);
            B[6] = 7 * xi * xi * xi * xi * xi * xi * (1 - xi);
            B[7] = xi * xi * xi * xi * xi * xi * xi;
            break;
        case 8:
            B[0] = (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[1] = 8 * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[2] = 28 * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[3] = 56 * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[4] = 70 * xi * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi) * (1 - xi);
            B[5] = 56 * xi * xi * xi * xi * xi * (1 - xi) * (1 - xi) * (1 - xi);
            B[6] = 28 * xi * xi * xi * xi * xi * xi * (1 - xi) * (1 - xi);
            B[7] = 8 * xi * xi * xi * xi * xi * xi * xi * (1 - xi);
            B[8] = xi * xi * xi * xi * xi * xi * xi * xi;
            break;
        default:
            std::cerr << "Error: Invalid polynomial degree" << std::endl;
            exit(1);
    }
    return B;
}

std::vector<double> BernsteinBasis::GetBernsteinBasisDerivativeSingleQP(const double &xi) const
{
    std::vector<double> dB(p + 1, 0.0);
    switch (p) {
        case 1:
            dB[0] = -1;
            dB[1] = 1;
            break;
        case 2:
            dB[0] = 2 * (xi - 1);
            dB[1] = 2 - 4 * xi;
            dB[2] = 2 * xi;
            break;
        case 3:
            dB[0] = -3 * (1 - xi) * (1 - xi);
            dB[1] = 9 * xi * xi - 12 * xi + 3;
            dB[2] = 6 * xi - 9 * xi * xi;
            dB[3] = 3 * xi * xi;
            break;
        case 4:
            dB[0] = 4 * (xi - 1) * (xi - 1) * (xi - 1);
            dB[1] = -16 * xi * xi * xi + 36 * xi * xi - 24 * xi + 4;
            dB[2] = 24 * xi * xi * xi - 36 * xi * xi + 12 * xi;
            dB[3] = 12 * xi * xi - 16 * xi * xi * xi;
            dB[4] = 4 * xi * xi * xi;
            break;
        case 5:
            dB[0] = -5 * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[1] = 5 * (5 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[2] = -10 * xi * (5 * xi - 2) * (xi - 1) * (xi - 1);
            dB[3] = 10 * xi * xi * (5 * xi * xi - 8 * xi + 3);
            dB[4] = -5 * xi * xi * xi * (5 * xi - 4);
            dB[5] = 5 * xi * xi * xi * xi;
            break;
        case 6:
            dB[0] = 6 * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[1] = -6 * (6 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[2] = 30 * xi * (3 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[3] = -60 * xi * xi * (2 * xi - 1) * (xi - 1) * (xi - 1);
            dB[4] = 30 * xi * xi * xi * (3 * xi * xi - 5 * xi + 2);
            dB[5] = -6 * xi * xi * xi * xi * (6 * xi - 5);
            dB[6] = 6 * xi * xi * xi * xi * xi;
            break;
        case 7:
            dB[0] = -7 * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[1] = 7 * (7 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[2] = -21 * xi * (7 * xi - 2) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[3] = 35 * xi * xi * (7 * xi - 3) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[4] = -35 * xi * xi * xi * (7 * xi - 4) * (xi - 1) * (xi - 1);
            dB[5] = 21 * xi * xi * xi * xi * (7 * xi * xi - 12 * xi + 5);
            dB[6] = -7 * xi * xi * xi * xi * xi * (7 * xi - 6);
            dB[7] = 7 * xi * xi * xi * xi * xi * xi;
            break;
        case 8:
            dB[0] = 8 * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[1] = -8 * (8 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[2] = 56 * xi * (4 * xi - 1) * (xi - 1) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[3] = -56 * xi * xi * (8 * xi - 3) * (xi - 1) * (xi - 1) * (xi - 1);
            dB[4] = 280 * xi * xi * xi * (2 * xi - 1) * (xi - 1) * (xi - 1);
            dB[5] = -56 * xi * xi * xi * xi * (8 * xi - 5) * (xi - 1);
            dB[6] = 56 * xi * xi * xi * xi * xi * (4 * xi * xi - 7 * xi + 3);
            dB[7] = -8 * xi * xi * xi * xi * xi * xi * (8 * xi - 7);
            dB[8] = 8 * xi * xi * xi * xi * xi * xi * xi;
            break;
        default:
            std::cerr << "Error: Invalid polynomial degree" << std::endl;
            exit(1);
    }
    return dB;
}