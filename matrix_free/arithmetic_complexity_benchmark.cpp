#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <petscsys.h>

#include "Element.hpp"
#include "ElementMFSF.hpp"
#include "FileManager.hpp"
#include "QuadraturePoint.hpp"

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
            "arithmetic_complexity_benchmark requires one MPI rank.\n");
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
        PetscPrintf(PETSC_COMM_WORLD, "This benchmark requires p=q.\n");
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

    const int n = p+1;
    const int n_basis = n*n;
    const int n_quadrature = n*n;
    const int n_elements = n_local_elem_x*n_local_elem_y;
    const int representative_x = n_local_elem_x/2;
    const int representative_y = n_local_elem_y/2;
    const int representative_element = representative_y*n_local_elem_x+representative_x;

    std::vector<double> element_control_points(2*n_basis, 0.0);
    for (int i = 0; i < n_basis; ++i)
    {
        const int global_basis = ien[representative_element*n_basis+i];
        element_control_points[2*i] = control_points[2*global_basis];
        element_control_points[2*i+1] = control_points[2*global_basis+1];
    }
    std::vector<double> element_extraction_x(
        extraction_x.begin()+representative_x*n*n,
        extraction_x.begin()+(representative_x+1)*n*n);
    std::vector<double> element_extraction_y(
        extraction_y.begin()+representative_y*n*n,
        extraction_y.begin()+(representative_y+1)*n*n);

    QuadraturePoint quad_x(n, 0, 1);
    QuadraturePoint quad_y(n, 0, 1);
    const std::vector<double> points_x = quad_x.GetQuadraturePoint();
    const std::vector<double> points_y = quad_y.GetQuadraturePoint();
    const std::vector<double> weights_x = quad_x.GetWeight();
    const std::vector<double> weights_y = quad_y.GetWeight();

    // Precompute the complete basis-pair product used by naive assembly.
    Element naive_element(p, q);
    naive_element.SetElement(element_extraction_x, element_extraction_y,
        elem_size_x[representative_x], elem_size_y[representative_y]);
    // Pair-major storage keeps each quadrature reduction contiguous in memory.
    std::vector<double> basis_pair(n_basis*n_basis*n_quadrature, 0.0);
    for (int qy = 0; qy < n; ++qy)
    {
        for (int qx = 0; qx < n; ++qx)
        {
            const int qp = qy*n+qx;
            std::vector<double> basis;
            std::vector<double> basis_dx;
            std::vector<double> basis_dy;
            double x = 0.0;
            double y = 0.0;
            double jacobian = 0.0;
            naive_element.GenerateElementSingleQP(points_x[qx], points_y[qy],
                element_control_points, basis, basis_dx, basis_dy,
                x, y, jacobian);
            for (int i = 0; i < n_basis; ++i)
            {
                for (int j = 0; j < n_basis; ++j)
                {
                    const int pair = i*n_basis+j;
                    basis_pair[pair*n_quadrature+qp] =
                        -(basis_dx[i]*basis_dx[j] + basis_dy[i]*basis_dy[j]);
                }
            }
        }
    }

    // Precompute all tensor bases and geometry factors used by MFSF.
    ElementMFSF mfsf_element(p, q);
    mfsf_element.SetElement(element_extraction_x, element_extraction_y,
        elem_size_x[representative_x], elem_size_y[representative_y]);
    std::vector<double> basis_x;
    std::vector<double> basis_y;
    std::vector<double> basis_derivative_x;
    std::vector<double> basis_derivative_y;
    std::vector<double> weight_sum;
    std::vector<double> jacobian;
    std::vector<double> weight_derivative_x;
    std::vector<double> weight_derivative_y;
    std::vector<double> dxi_dx;
    std::vector<double> dxi_dy;
    std::vector<double> deta_dx;
    std::vector<double> deta_dy;
    mfsf_element.GenerateElement(&quad_x, &quad_y, element_control_points,
        basis_x, basis_y, basis_derivative_x, basis_derivative_y,
        weight_sum, jacobian, weight_derivative_x, weight_derivative_y,
        dxi_dx, dxi_dy, deta_dx, deta_dy);

    std::vector<double> quadrature_scale(n_elements*n_quadrature, 0.0);
    std::vector<double> d00(n_elements*n_quadrature, 0.0);
    std::vector<double> d01(n_elements*n_quadrature, 0.0);
    std::vector<double> d02(n_elements*n_quadrature, 0.0);
    std::vector<double> d11(n_elements*n_quadrature, 0.0);
    std::vector<double> d12(n_elements*n_quadrature, 0.0);
    std::vector<double> d22(n_elements*n_quadrature, 0.0);
    for (int element = 0; element < n_elements; ++element)
    {
        const double geometry_factor = 1.0 + 1.0e-12*(element%17);
        for (int qp = 0; qp < n_quadrature; ++qp)
        {
            const int qx = qp%n;
            const int qy = qp/n;
            quadrature_scale[element*n_quadrature+qp] = geometry_factor
                * jacobian[qp]*weights_x[qx]*weights_y[qy];
            const double scale = quadrature_scale[element*n_quadrature+qp];
            const double metric00 = scale
                * (dxi_dx[qp]*dxi_dx[qp] + dxi_dy[qp]*dxi_dy[qp]);
            const double metric01 = scale
                * (dxi_dx[qp]*deta_dx[qp] + dxi_dy[qp]*deta_dy[qp]);
            const double metric11 = scale
                * (deta_dx[qp]*deta_dx[qp] + deta_dy[qp]*deta_dy[qp]);
            const int index = element*n_quadrature+qp;
            const double inverse_weight = 1.0/weight_sum[qp];
            const double inverse_weight_squared = inverse_weight*inverse_weight;
            const double inverse_weight_cubed =
                inverse_weight_squared*inverse_weight;
            const double inverse_weight_fourth =
                inverse_weight_squared*inverse_weight_squared;
            const double dw_xi = weight_derivative_x[qp];
            const double dw_eta = weight_derivative_y[qp];
            d00[index] =
                (dw_xi*dw_xi*metric00 + 2.0*dw_xi*dw_eta*metric01
                 + dw_eta*dw_eta*metric11) * inverse_weight_fourth;
            d01[index] = -(dw_xi*metric00 + dw_eta*metric01)
                * inverse_weight_cubed;
            d02[index] = -(dw_xi*metric01 + dw_eta*metric11)
                * inverse_weight_cubed;
            d11[index] = metric00*inverse_weight_squared;
            d12[index] = metric01*inverse_weight_squared;
            d22[index] = metric11*inverse_weight_squared;
        }
    }

    constexpr int input_variants = 16;
    std::vector<double> inputs(input_variants*n_basis, 0.0);
    for (int variant = 0; variant < input_variants; ++variant)
        for (int i = 0; i < n_basis; ++i)
            inputs[variant*n_basis+i] = std::sin(0.17*(i+1+variant))
                + 0.25*std::cos(0.11*(i+1+2*variant));

    std::vector<double> stiffness(n_basis*n_basis, 0.0);
    std::vector<double> value_x(n_basis, 0.0);
    std::vector<double> derivative_x_stage(n_basis, 0.0);
    std::vector<double> scale_xi(n_quadrature, 0.0);
    std::vector<double> scale_eta(n_quadrature, 0.0);
    std::vector<double> scale_value(n_quadrature, 0.0);
    std::vector<double> backward_value(n_basis, 0.0);
    std::vector<double> backward_eta(n_basis, 0.0);
    std::vector<double> output(n_basis, 0.0);
    volatile double sink = 0.0;

    // Timed region: multiplication and accumulation only; all basis products are ready.
    PetscLogDouble start = 0.0;
    PetscLogDouble end = 0.0;
    PetscTime(&start);
    for (int element = 0; element < n_elements; ++element)
    {
        const double *element_scale = quadrature_scale.data()+element*n_quadrature;
        for (int i = 0; i < n_basis; ++i)
        {
            for (int j = 0; j < n_basis; ++j)
            {
                const int pair = i*n_basis+j;
                const double *pair_values =
                    basis_pair.data()+pair*n_quadrature;
                double value = 0.0;
                for (int qp = 0; qp < n_quadrature; ++qp)
                {
                    value += element_scale[qp]*pair_values[qp];
                }
                stiffness[pair] = value;
            }
        }
        sink += stiffness[element%(n_basis*n_basis)];
    }
    PetscTime(&end);
    const double naive_seconds = end-start;

    // Timed region: forward contractions, pointwise products, backward contractions.
    PetscTime(&start);
    for (int element = 0; element < n_elements; ++element)
    {
        const double *input = inputs.data()+(element%input_variants)*n_basis;
        for (int qx = 0; qx < n; ++qx)
        {
            for (int j = 0; j < n; ++j)
            {
                double value = 0.0;
                double derivative = 0.0;
                for (int i = 0; i < n; ++i)
                {
                    value += basis_x[qx*n+i]*input[j*n+i];
                    derivative += basis_derivative_x[qx*n+i]*input[j*n+i];
                }
                value_x[qx*n+j] = value;
                derivative_x_stage[qx*n+j] = derivative;
            }
        }

        for (int qy = 0; qy < n; ++qy)
        {
            for (int qx = 0; qx < n; ++qx)
            {
                const int qp = qy*n+qx;
                double value = 0.0;
                double derivative_xi = 0.0;
                double derivative_eta = 0.0;
                for (int j = 0; j < n; ++j)
                {
                    value += basis_y[qy*n+j]*value_x[qx*n+j];
                    derivative_xi += basis_y[qy*n+j]*derivative_x_stage[qx*n+j];
                    derivative_eta += basis_derivative_y[qy*n+j]*value_x[qx*n+j];
                }
                const int d_index = element*n_quadrature+qp;
                scale_value[qp] = d00[d_index]*value
                    + d01[d_index]*derivative_xi + d02[d_index]*derivative_eta;
                scale_xi[qp] = d01[d_index]*value
                    + d11[d_index]*derivative_xi + d12[d_index]*derivative_eta;
                scale_eta[qp] = d02[d_index]*value
                    + d12[d_index]*derivative_xi + d22[d_index]*derivative_eta;
            }
        }

        for (int qy = 0; qy < n; ++qy)
        {
            for (int i = 0; i < n; ++i)
            {
                double value = 0.0;
                double eta_value = 0.0;
                for (int qx = 0; qx < n; ++qx)
                {
                    const int qp = qy*n+qx;
                    value += basis_x[qx*n+i]*scale_value[qp]
                        + basis_derivative_x[qx*n+i]*scale_xi[qp];
                    eta_value += basis_x[qx*n+i]*scale_eta[qp];
                }
                backward_value[qy*n+i] = value;
                backward_eta[qy*n+i] = eta_value;
            }
        }

        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                double value = 0.0;
                for (int qy = 0; qy < n; ++qy)
                {
                    value += basis_y[qy*n+j]*backward_value[qy*n+i]
                        + basis_derivative_y[qy*n+j]*backward_eta[qy*n+i];
                }
                output[j*n+i] = -value;
            }
        }
        sink += output[element%n_basis];
    }
    PetscTime(&end);
    const double mfsf_seconds = end-start;

    const int last_element = n_elements-1;
    const double *validation_input =
        inputs.data()+(last_element%input_variants)*n_basis;
    std::vector<double> reference(n_basis, 0.0);
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j)
            reference[i] += stiffness[i*n_basis+j]*validation_input[j];
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    for (int i = 0; i < n_basis; ++i)
    {
        const double difference = reference[i]-output[i];
        difference_squared += difference*difference;
        reference_squared += reference[i]*reference[i];
    }
    const double relative_error = std::sqrt(difference_squared/reference_squared);

    std::cout << std::setprecision(15)
              << "p: " << p << '\n'
              << "n_basis_1d: " << n << '\n'
              << "arithmetic kernel calls: " << n_elements << '\n'
              << "naive arithmetic seconds: " << naive_seconds << '\n'
              << "MFSF arithmetic seconds: " << mfsf_seconds << '\n'
              << "MFSF/naive arithmetic ratio: " << mfsf_seconds/naive_seconds << '\n'
              << "relative operator error: " << relative_error << '\n'
              << "checksum: " << sink << '\n';
    std::cout << "ARITHMETIC_CSV,"
              << p << ',' << n << ',' << n_elem_x << ',' << n_elem_y << ','
              << n_elements << ',' << naive_seconds << ','
              << mfsf_seconds << ',' << mfsf_seconds/naive_seconds << ','
              << relative_error << '\n';

    PetscFinalize();
    return 0;
}
