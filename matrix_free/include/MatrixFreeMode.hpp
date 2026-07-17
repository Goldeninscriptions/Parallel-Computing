#ifndef MATRIXFREEMODE_HPP
#define MATRIXFREEMODE_HPP

#include <iostream>
#include <string>

enum class MatrixFreeMode
{
    MF,
    MFSF
};

inline const char * MatrixFreeModeName(const MatrixFreeMode mode)
{
    return (mode == MatrixFreeMode::MF) ? "MF" : "MFSF";
}

inline bool ParseMatrixFreeMode(int argc, char *argv[], MatrixFreeMode &mode)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "--mf" || arg == "--mf-mode=mf" || arg == "--matrix-free-mode=mf")
        {
            mode = MatrixFreeMode::MF;
        }
        else if (arg == "--mfsf" || arg == "--mf-mode=mfsf" || arg == "--matrix-free-mode=mfsf")
        {
            mode = MatrixFreeMode::MFSF;
        }
        else if (arg == "--help" || arg == "-h")
        {
            std::cout << "Options:\n"
                      << "  --mf                  Use standard matrix-free mode\n"
                      << "  --mfsf                Use MFSF matrix-free mode\n"
                      << "  --mf-mode=mf|mfsf     Select matrix-free mode\n";
            return false;
        }
    }

    return true;
}

#endif
