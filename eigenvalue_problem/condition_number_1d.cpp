#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "BSplineBasis.hpp"
#include "QuadraturePoint.hpp"

namespace
{
std::vector<double> BuildOpenUniformKnotVector(int degree, int nElem, double length)
{
    std::vector<double> knots;
    knots.reserve(nElem + 2 * degree + 1);

    for (int i = 0; i <= degree; ++i)
        knots.push_back(0.0);

    const double h = length / static_cast<double>(nElem);
    for (int i = 1; i < nElem; ++i)
        knots.push_back(i * h);

    for (int i = 0; i <= degree; ++i)
        knots.push_back(length);

    return knots;
}

Eigen::MatrixXd AssembleStiffness1D(const BSplineBasis &basis, int nElem, double length)
{
    const int degree = basis.GetDegree();
    const int nFunc = static_cast<int>(basis.GetKnotVector().size()) - degree - 1;
    const double h = length / static_cast<double>(nElem);

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nFunc, nFunc);
    QuadraturePoint quad(degree + 1, -1, 1);
    const std::vector<double> qpts = quad.GetQuadraturePoint();
    const std::vector<double> weights = quad.GetWeight();

    for (int elem = 0; elem < nElem; ++elem)
    {
        const double x_left = elem * h;
        const double x_right = (elem + 1) * h;
        const double jacobian = 0.5 * (x_right - x_left);

        for (int iq = 0; iq < quad.GetNumQuadraturePoint(); ++iq)
        {
            const double x = 0.5 * ((x_right - x_left) * qpts[iq] + x_right + x_left);
            int span = basis.FindSpan(x);
            if (span == nFunc)
                span = nFunc - 1;

            const std::vector<double> dN = basis.DerBasisFuns(x, span, 1);

            for (int a = 0; a <= degree; ++a)
            {
                const int A = span - degree + a;
                for (int b = 0; b <= degree; ++b)
                {
                    const int B = span - degree + b;
                    K(A, B) += weights[iq] * jacobian * dN[a] * dN[b];
                }
            }
        }
    }

    return K;
}

Eigen::MatrixXd ApplyBoundaryCondition(const Eigen::MatrixXd &K, const std::string &bc)
{
    const int n = static_cast<int>(K.rows());
    if (bc == "dirichlet-dirichlet")
    {
        return K.block(1, 1, n - 2, n - 2);
    }
    if (bc == "dirichlet-neumann")
    {
        return K.block(1, 1, n - 1, n - 1);
    }

    std::cerr << "Unsupported bc_type: " << bc << std::endl;
    std::exit(1);
}

void ParseArgs(int argc, char **argv, int &degree, int &nElem, double &length, std::string &bc)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--degree" && i + 1 < argc)
        {
            degree = std::stoi(argv[++i]);
        }
        else if (arg == "--nelem" && i + 1 < argc)
        {
            nElem = std::stoi(argv[++i]);
        }
        else if (arg == "--length" && i + 1 < argc)
        {
            length = std::stod(argv[++i]);
        }
        else if (arg == "--bc" && i + 1 < argc)
        {
            bc = argv[++i];
        }
        else
        {
            std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
            std::exit(1);
        }
    }
}
}

int main(int argc, char **argv)
{
    int degree = 3;
    int nElem = 16;
    double length = 1.0;
    std::string bc = "dirichlet-dirichlet";

    ParseArgs(argc, argv, degree, nElem, length, bc);

    if (degree < 1 || nElem < 1 || length <= 0.0)
    {
        std::cerr << "Invalid input parameters." << std::endl;
        return 1;
    }

    const std::vector<double> knots = BuildOpenUniformKnotVector(degree, nElem, length);
    BSplineBasis basis(degree, knots);
    const Eigen::MatrixXd K_full = AssembleStiffness1D(basis, nElem, length);
    const Eigen::MatrixXd K = ApplyBoundaryCondition(K_full, bc);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(K);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Eigen decomposition failed." << std::endl;
        return 1;
    }

    const Eigen::VectorXd eigvals = solver.eigenvalues();
    const double tol = std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max<Eigen::Index>(1, K.rows()));

    double lambda_min = -1.0;
    for (int i = 0; i < eigvals.size(); ++i)
    {
        if (eigvals[i] > tol)
        {
            lambda_min = eigvals[i];
            break;
        }
    }

    if (lambda_min <= 0.0)
    {
        std::cerr << "Failed to locate a positive minimum eigenvalue." << std::endl;
        return 1;
    }

    const double lambda_max = eigvals[eigvals.size() - 1];
    const double condition_number = lambda_max / lambda_min;

    std::cout << "degree: " << degree << std::endl;
    std::cout << "nElem: " << nElem << std::endl;
    std::cout << "length: " << length << std::endl;
    std::cout << "bc_type: " << bc << std::endl;
    std::cout << "num_basis: " << K_full.rows() << std::endl;
    std::cout << "system_size: " << K.rows() << std::endl;
    std::cout << "lambda_min: " << lambda_min << std::endl;
    std::cout << "lambda_max: " << lambda_max << std::endl;
    std::cout << "condition_number: " << condition_number << std::endl;

    return 0;
}
