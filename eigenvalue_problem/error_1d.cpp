#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "BSplineBasis.hpp"
#include "QuadraturePoint.hpp"

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;

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

double ExactSolution(double x)
{
    return std::sin(2.0 * kPi * x) + std::exp(x) - 1.0;
}

double ExactDerivative(double x)
{
    return 2.0 * kPi * std::cos(2.0 * kPi * x) + std::exp(x);
}

double SourceTerm(double x)
{
    return 4.0 * kPi * kPi * std::sin(2.0 * kPi * x) - std::exp(x);
}

double LeftDirichletValue()
{
    return ExactSolution(0.0);
}

double RightDirichletValue()
{
    return ExactSolution(1.0);
}

double RightNeumannValue()
{
    return ExactDerivative(1.0);
}

int ClampSpan(const BSplineBasis &basis, double x)
{
    const int nFunc = static_cast<int>(basis.GetKnotVector().size()) - basis.GetDegree() - 1;
    int span = basis.FindSpan(x);
    if (span == nFunc)
        span = nFunc - 1;
    return span;
}

Eigen::MatrixXd AssembleStiffness1D(const BSplineBasis &basis, int nElem, double length)
{
    const int degree = basis.GetDegree();
    const int nFunc = static_cast<int>(basis.GetKnotVector().size()) - degree - 1;
    const double h = length / static_cast<double>(nElem);

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nFunc, nFunc);
    QuadraturePoint quad(std::min(10, degree + 3), -1, 1);
    const std::vector<double> qpts = quad.GetQuadraturePoint();
    const std::vector<double> weights = quad.GetWeight();

    for (int elem = 0; elem < nElem; ++elem)
    {
        const double xLeft = elem * h;
        const double xRight = (elem + 1) * h;
        const double jacobian = 0.5 * (xRight - xLeft);

        for (int iq = 0; iq < quad.GetNumQuadraturePoint(); ++iq)
        {
            const double x = 0.5 * ((xRight - xLeft) * qpts[iq] + xRight + xLeft);
            const int span = ClampSpan(basis, x);
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

Eigen::VectorXd AssembleLoad1D(const BSplineBasis &basis, int nElem, double length, const std::string &bc)
{
    const int degree = basis.GetDegree();
    const int nFunc = static_cast<int>(basis.GetKnotVector().size()) - degree - 1;
    const double h = length / static_cast<double>(nElem);

    Eigen::VectorXd F = Eigen::VectorXd::Zero(nFunc);
    QuadraturePoint quad(std::min(10, degree + 5), -1, 1);
    const std::vector<double> qpts = quad.GetQuadraturePoint();
    const std::vector<double> weights = quad.GetWeight();

    for (int elem = 0; elem < nElem; ++elem)
    {
        const double xLeft = elem * h;
        const double xRight = (elem + 1) * h;
        const double jacobian = 0.5 * (xRight - xLeft);

        for (int iq = 0; iq < quad.GetNumQuadraturePoint(); ++iq)
        {
            const double x = 0.5 * ((xRight - xLeft) * qpts[iq] + xRight + xLeft);
            const int span = ClampSpan(basis, x);
            const std::vector<double> N = basis.BasisFuns(x, span);
            const double rhsValue = SourceTerm(x);

            for (int a = 0; a <= degree; ++a)
            {
                const int A = span - degree + a;
                F[A] += weights[iq] * jacobian * rhsValue * N[a];
            }
        }
    }

    if (bc == "dirichlet-neumann")
    {
        F[nFunc - 1] += RightNeumannValue();
        return F;
    }

    if (bc == "dirichlet-dirichlet")
        return F;

    std::cerr << "Unsupported bc_type: " << bc << std::endl;
    std::exit(1);
}

Eigen::VectorXd SolveSystem(
    const Eigen::MatrixXd &KFull,
    const Eigen::VectorXd &FFull,
    const std::string &bc)
{
    const int n = static_cast<int>(KFull.rows());
    const double gLeft = LeftDirichletValue();
    Eigen::VectorXd U = Eigen::VectorXd::Zero(n);
    U[0] = gLeft;

    if (bc == "dirichlet-dirichlet")
    {
        const double gRight = RightDirichletValue();
        U[n - 1] = gRight;

        Eigen::VectorXd rhs = FFull.segment(1, n - 2);
        rhs -= KFull.block(1, 0, n - 2, 1) * gLeft;
        rhs -= KFull.block(1, n - 1, n - 2, 1) * gRight;

        const Eigen::MatrixXd K = KFull.block(1, 1, n - 2, n - 2);
        const Eigen::VectorXd inner = K.ldlt().solve(rhs);
        U.segment(1, n - 2) = inner;
        return U;
    }

    if (bc == "dirichlet-neumann")
    {
        Eigen::VectorXd rhs = FFull.segment(1, n - 1);
        rhs -= KFull.block(1, 0, n - 1, 1) * gLeft;

        const Eigen::MatrixXd K = KFull.block(1, 1, n - 1, n - 1);
        const Eigen::VectorXd inner = K.ldlt().solve(rhs);
        U.segment(1, n - 1) = inner;
        return U;
    }

    std::cerr << "Unsupported bc_type: " << bc << std::endl;
    std::exit(1);
}

void ComputeErrors(
    const BSplineBasis &basis,
    const Eigen::VectorXd &coeff,
    int nElem,
    double length,
    double &l2Error,
    double &h1Error)
{
    const int degree = basis.GetDegree();
    const double h = length / static_cast<double>(nElem);
    double l2Accum = 0.0;
    double h1Accum = 0.0;

    QuadraturePoint quad(std::min(10, degree + 5), -1, 1);
    const std::vector<double> qpts = quad.GetQuadraturePoint();
    const std::vector<double> weights = quad.GetWeight();

    for (int elem = 0; elem < nElem; ++elem)
    {
        const double xLeft = elem * h;
        const double xRight = (elem + 1) * h;
        const double jacobian = 0.5 * (xRight - xLeft);

        for (int iq = 0; iq < quad.GetNumQuadraturePoint(); ++iq)
        {
            const double x = 0.5 * ((xRight - xLeft) * qpts[iq] + xRight + xLeft);
            const int span = ClampSpan(basis, x);
            const std::vector<double> N = basis.BasisFuns(x, span);
            const std::vector<double> dN = basis.DerBasisFuns(x, span, 1);

            double uh = 0.0;
            double duh = 0.0;
            for (int a = 0; a <= degree; ++a)
            {
                const int A = span - degree + a;
                uh += N[a] * coeff[A];
                duh += dN[a] * coeff[A];
            }

            const double u = ExactSolution(x);
            const double du = ExactDerivative(x);
            l2Accum += weights[iq] * jacobian * (u - uh) * (u - uh);
            h1Accum += weights[iq] * jacobian * (du - duh) * (du - duh);
        }
    }

    l2Error = std::sqrt(l2Accum);
    h1Error = std::sqrt(h1Accum);
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
    std::string bc = "dirichlet-neumann";

    ParseArgs(argc, argv, degree, nElem, length, bc);

    if (degree < 1 || nElem < 1 || length <= 0.0)
    {
        std::cerr << "Invalid input parameters." << std::endl;
        return 1;
    }

    const std::vector<double> knots = BuildOpenUniformKnotVector(degree, nElem, length);
    const BSplineBasis basis(degree, knots);
    const Eigen::MatrixXd K = AssembleStiffness1D(basis, nElem, length);
    const Eigen::VectorXd F = AssembleLoad1D(basis, nElem, length, bc);
    const Eigen::VectorXd U = SolveSystem(K, F, bc);

    double l2Error = 0.0;
    double h1Error = 0.0;
    ComputeErrors(basis, U, nElem, length, l2Error, h1Error);

    std::cout << "degree: " << degree << std::endl;
    std::cout << "nElem: " << nElem << std::endl;
    std::cout << "length: " << length << std::endl;
    std::cout << "bc_type: " << bc << std::endl;
    std::cout << "num_basis: " << U.size() << std::endl;
    std::cout << "l2_error: " << l2Error << std::endl;
    std::cout << "h1_error: " << h1Error << std::endl;

    return 0;
}
