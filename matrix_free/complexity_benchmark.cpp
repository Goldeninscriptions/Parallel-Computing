#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <petscsys.h>

#include "Element.hpp"
#include "ElementMFSF.hpp"
#include "FileManager.hpp"
#include "LocalAssemblyMFSF.hpp"
#include "QuadraturePoint.hpp"

namespace
{
struct ElementData
{
    std::vector<double> control_points;
    std::vector<double> extraction_x;
    std::vector<double> extraction_y;
    std::vector<PetscScalar> input;
    double hx;
    double hy;
};

void AssembleNaiveStiffness(const Element &element,
    const QuadraturePoint &quad_x, const QuadraturePoint &quad_y,
    const std::vector<double> &control_points,
    std::vector<PetscScalar> &stiffness)
{
    const std::vector<double> points_x = quad_x.GetQuadraturePoint();
    const std::vector<double> points_y = quad_y.GetQuadraturePoint();
    const std::vector<double> weights_x = quad_x.GetWeight();
    const std::vector<double> weights_y = quad_y.GetWeight();
    const int n_basis = element.GetNumLocalBasis();

    std::fill(stiffness.begin(), stiffness.end(), 0.0);
    for (int qy = 0; qy < static_cast<int>(points_y.size()); ++qy)
    {
        for (int qx = 0; qx < static_cast<int>(points_x.size()); ++qx)
        {
            std::vector<double> basis;
            std::vector<double> basis_dx;
            std::vector<double> basis_dy;
            double x = 0.0;
            double y = 0.0;
            double jacobian = 0.0;
            element.GenerateElementSingleQP(points_x[qx], points_y[qy],
                control_points, basis, basis_dx, basis_dy, x, y, jacobian);
            const double scale = weights_x[qx] * weights_y[qy] * jacobian;

            for (int i = 0; i < n_basis; ++i)
            {
                for (int j = 0; j < n_basis; ++j)
                {
                    stiffness[i*n_basis+j] -= scale
                        * (basis_dx[i]*basis_dx[j] + basis_dy[i]*basis_dy[j]);
                }
            }
        }
    }
}

double RelativeError(const std::vector<PetscScalar> &reference,
    const PetscScalar *candidate)
{
    double difference_norm_squared = 0.0;
    double reference_norm_squared = 0.0;
    for (int i = 0; i < static_cast<int>(reference.size()); ++i)
    {
        const double difference = PetscRealPart(reference[i] - candidate[i]);
        const double value = PetscRealPart(reference[i]);
        difference_norm_squared += difference*difference;
        reference_norm_squared += value*value;
    }
    return std::sqrt(difference_norm_squared / reference_norm_squared);
}
}

int main(int argc, char *argv[])
{
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    PetscMPIInt rank = 0;
    PetscMPIInt size = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size != 1)
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "complexity_benchmark requires exactly one MPI rank.\n");
        PetscFinalize();
        return 1;
    }

    PetscInt repeat_naive = 1;
    PetscInt repeat_mfsf = 10;
    PetscOptionsGetInt(nullptr, nullptr, "-repeat_naive", &repeat_naive, nullptr);
    PetscOptionsGetInt(nullptr, nullptr, "-repeat_mfsf", &repeat_mfsf, nullptr);
    if (repeat_naive < 1 || repeat_mfsf < 1)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Repeat counts must be positive.\n");
        PetscFinalize();
        return 1;
    }

    int p = 0;
    int q = 0;
    int n_elem_x = 0;
    int n_elem_y = 0;
    int part_num_1d = 0;
    int dim = 0;
    double length_x = 0.0;
    double length_y = 0.0;
    std::string base_name;
    FileManager file_manager;
    file_manager.ReadPreprocessInfo("info.txt", p, q, length_x, length_y,
        n_elem_x, n_elem_y, part_num_1d, dim, base_name);
    if (p != q)
    {
        PetscPrintf(PETSC_COMM_WORLD,
            "This experiment requires equal polynomial degrees p and q.\n");
        PetscFinalize();
        return 1;
    }

    int n_local_functions = 0;
    int n_local_elem_x = 0;
    int n_local_elem_y = 0;
    std::vector<int> ghost_ids;
    std::vector<double> control_points;
    std::vector<int> ids;
    std::vector<int> ien;
    std::vector<int> dirichlet;
    std::vector<double> elem_size_x;
    std::vector<double> elem_size_y;
    std::vector<double> extraction_x;
    std::vector<double> extraction_y;
    file_manager.ReadPartition(file_manager.GetPartitionFilename(base_name, rank),
        n_local_functions, n_local_elem_x, n_local_elem_y,
        elem_size_x, elem_size_y, control_points, ids, ghost_ids, dirichlet, ien,
        extraction_x, extraction_y);

    const int nx = p+1;
    const int ny = q+1;
    const int n_basis = nx*ny;
    const int n_elements = n_local_elem_x*n_local_elem_y;
    std::vector<ElementData> elements;
    elements.reserve(n_elements);
    for (int ey = 0; ey < n_local_elem_y; ++ey)
    {
        for (int ex = 0; ex < n_local_elem_x; ++ex)
        {
            const int element_index = ey*n_local_elem_x+ex;
            ElementData data;
            data.control_points.resize(2*n_basis);
            data.extraction_x.assign(
                extraction_x.begin()+ex*nx*nx,
                extraction_x.begin()+(ex+1)*nx*nx);
            data.extraction_y.assign(
                extraction_y.begin()+ey*ny*ny,
                extraction_y.begin()+(ey+1)*ny*ny);
            data.input.resize(n_basis);
            data.hx = elem_size_x[ex];
            data.hy = elem_size_y[ey];
            for (int i = 0; i < n_basis; ++i)
            {
                const int global_basis = ien[element_index*n_basis+i];
                data.control_points[2*i] = control_points[2*global_basis];
                data.control_points[2*i+1] = control_points[2*global_basis+1];
                data.input[i] = std::sin(0.17*(i+1)) + 0.25*std::cos(0.11*(i+1));
            }
            elements.push_back(std::move(data));
        }
    }

    Element naive_element(p, q);
    ElementMFSF mfsf_element(p, q);
    LocalAssemblyMFSF mfsf_assembly(p, q);
    QuadraturePoint quad_x(nx, 0, 1);
    QuadraturePoint quad_y(ny, 0, 1);
    std::vector<PetscScalar> stiffness(n_basis*n_basis, 0.0);

    const ElementData &validation_data = elements.front();
    naive_element.SetElement(validation_data.extraction_x,
        validation_data.extraction_y, validation_data.hx, validation_data.hy);
    AssembleNaiveStiffness(naive_element, quad_x, quad_y,
        validation_data.control_points, stiffness);
    std::vector<PetscScalar> reference(n_basis, 0.0);
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j)
            reference[i] += stiffness[i*n_basis+j] * validation_data.input[j];

    mfsf_element.SetElement(validation_data.extraction_x,
        validation_data.extraction_y, validation_data.hx, validation_data.hy);
    std::copy(validation_data.input.begin(), validation_data.input.end(),
        mfsf_assembly.Floc_in);
    mfsf_assembly.LocalMatMulMF(&mfsf_element, validation_data.control_points);
    const double relative_error = RelativeError(reference, mfsf_assembly.Floc_out);

    auto run_naive = [&]()
    {
        double checksum = 0.0;
        for (const ElementData &data : elements)
        {
            naive_element.SetElement(data.extraction_x, data.extraction_y,
                data.hx, data.hy);
            AssembleNaiveStiffness(naive_element, quad_x, quad_y,
                data.control_points, stiffness);
            checksum += PetscRealPart(stiffness.front());
        }
        return checksum;
    };

    auto run_mfsf = [&]()
    {
        double checksum = 0.0;
        for (const ElementData &data : elements)
        {
            mfsf_element.SetElement(data.extraction_x, data.extraction_y,
                data.hx, data.hy);
            std::copy(data.input.begin(), data.input.end(), mfsf_assembly.Floc_in);
            mfsf_assembly.LocalMatMulMF(&mfsf_element, data.control_points);
            checksum += PetscRealPart(mfsf_assembly.Floc_out[0]);
        }
        return checksum;
    };

    run_naive();
    run_mfsf();

    PetscLogDouble start = 0.0;
    PetscLogDouble end = 0.0;
    double naive_checksum = 0.0;
    PetscTime(&start);
    for (PetscInt repeat = 0; repeat < repeat_naive; ++repeat)
        naive_checksum += run_naive();
    PetscTime(&end);
    const double naive_seconds = (end-start)/repeat_naive;

    double mfsf_checksum = 0.0;
    PetscTime(&start);
    for (PetscInt repeat = 0; repeat < repeat_mfsf; ++repeat)
        mfsf_checksum += run_mfsf();
    PetscTime(&end);
    const double mfsf_seconds = (end-start)/repeat_mfsf;

    std::cout << std::setprecision(15);
    std::cout << "p: " << p << '\n'
              << "n_basis_1d: " << nx << '\n'
              << "n_elements: " << n_elements << '\n'
              << "naive average seconds per mesh sweep: " << naive_seconds << '\n'
              << "MFSF average seconds per mesh sweep: " << mfsf_seconds << '\n'
              << "MFSF/naive time ratio: " << mfsf_seconds/naive_seconds << '\n'
              << "relative operator error: " << relative_error << '\n'
              << "naive checksum: " << naive_checksum << '\n'
              << "MFSF checksum: " << mfsf_checksum << '\n';
    std::cout << "BENCHMARK_CSV,"
              << p << ',' << nx << ',' << n_elem_x << ',' << n_elem_y << ','
              << n_elements << ',' << repeat_naive << ',' << repeat_mfsf << ','
              << naive_seconds << ',' << naive_seconds/n_elements << ','
              << mfsf_seconds << ',' << mfsf_seconds/n_elements << ','
              << mfsf_seconds/naive_seconds << ',' << relative_error << '\n';

    PetscFinalize();
    return 0;
}
