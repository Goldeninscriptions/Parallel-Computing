#include "ElementMFSF.hpp"

void ElementMFSF::GenerateBSplineBasis1D(const double &xi,
    std::vector<double> &N1, std::vector<double> &N2) const
{
    BernsteinBasis * bern1 = new BernsteinBasis(p);
    BernsteinBasis * bern2 = new BernsteinBasis(q);

    RefElement * ref = new RefElement();

    N1 = ref->GenerateBasis1DSingleQP(bern1, extraction1, xi);
    N2 = ref->GenerateBasis1DSingleQP(bern2, extraction2, xi);

    delete bern1;
    delete bern2;
    delete ref;
}

void ElementMFSF::GenerateBSplineBasis1D(const double &xi,
    std::vector<double> &N1, std::vector<double> &N2,
    std::vector<double> &dN1, std::vector<double> &dN2) const
{
    BernsteinBasis * bern1 = new BernsteinBasis(p);
    BernsteinBasis * bern2 = new BernsteinBasis(q);

    RefElement * ref = new RefElement();

    N1 = ref->GenerateBasis1DSingleQP(bern1, extraction1, xi);
    dN1 = ref->GenerateBasisDerivative1DSingleQP(bern1, extraction1, xi, hx);
    N2 = ref->GenerateBasis1DSingleQP(bern2, extraction2, xi);
    dN2 = ref->GenerateBasisDerivative1DSingleQP(bern2, extraction2, xi, hy);

    delete bern1;
    delete bern2;
    delete ref;
}

void ElementMFSF::GenerateElementSingleQP(const std::vector<double> &eCP,
    const std::vector<double> &N1, const std::vector<double> &N2,
    const std::vector<double> &dN1, const std::vector<double> &dN2,
    std::vector<double> &R,
    double &jacobian) const
{
    std::vector<double> N{};
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};

    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;

    for (int j = 0; j < q+1; ++j)
    {
        for (int i = 0; i < p+1; ++i)
        {
            N.push_back(N1[i] * N2[j]);
            w += N.back();
            dN_dxi.push_back(dN1[i] * N2[j]);
            dw_dxi += dN_dxi.back();
            dN_deta.push_back(N1[i] * dN2[j]);
            dw_deta += dN_deta.back();
        }
    }

    R.clear();
    std::vector<double> dR_dxi{};
    std::vector<double> dR_deta{};

    for (int j = 0; j < q+1; ++j)
    {
        for (int i = 0; i < p+1; ++i)
        {
            R.push_back(N[j*(p+1)+i]/w);
            dR_dxi.push_back((dN_dxi[j*(p+1)+i]-dw_dxi*R.back())/w);
            dR_deta.push_back((dN_deta[j*(p+1)+i]-dw_deta*R.back())/w);
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    jacobian = 0.0;

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            dx_dxi += eCP[2*(j*(p+1)+i)] * dR_dxi[j*(p+1)+i];
            dx_deta += eCP[2*(j*(p+1)+i)] * dR_deta[j*(p+1)+i];
            dy_dxi += eCP[2*(j*(p+1)+i)+1] * dR_dxi[j*(p+1)+i];
            dy_deta += eCP[2*(j*(p+1)+i)+1] * dR_deta[j*(p+1)+i];
        }
    }

    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;
    jacobian *= hx*hy;
}

void ElementMFSF::GenerateElementSingleQP(const std::vector<double> &eCP,
    const std::vector<double> &N1, const std::vector<double> &N2,
    const std::vector<double> &dN1, const std::vector<double> &dN2,
    double &w, double &jacobian,
    double &dw_dx, double &dw_dy) const
{
    std::vector<double> N{};
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            N.push_back(N1[i] * N2[j]);
            w += N.back();
            dN_dxi.push_back(dN1[i] * N2[j]);
            dw_dxi += dN_dxi.back();
            dN_deta.push_back(N1[i] * dN2[j]);
            dw_deta += dN_deta.back();
        }
    }

    std::vector<double> R{};
    std::vector<double> dR_dxi{};
    std::vector<double> dR_deta{};

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            R.push_back(N[j*(p+1)+i]/w);
            dR_dxi.push_back((dN_dxi[j*(p+1)+i]-dw_dxi*R.back())/w);
            dR_deta.push_back((dN_deta[j*(p+1)+i]-dw_deta*R.back())/w);
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    jacobian = 0.0;

    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    double x = 0.0;
    double y = 0.0;

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            x += eCP[2*(j*(p+1)+i)] * R[j*(p+1)+i];
            y += eCP[2*(j*(p+1)+i)+1] * R[j*(p+1)+i];
            dx_dxi += eCP[2*(j*(p+1)+i)] * dR_dxi[j*(p+1)+i];
            dx_deta += eCP[2*(j*(p+1)+i)] * dR_deta[j*(p+1)+i];
            dy_dxi += eCP[2*(j*(p+1)+i)+1] * dR_dxi[j*(p+1)+i];
            dy_deta += eCP[2*(j*(p+1)+i)+1] * dR_deta[j*(p+1)+i];
        }
    }

    jacobian = dx_dxi*dy_deta - dx_deta*dy_dxi;

    dxi_dx = dy_deta/jacobian;
    dxi_dy = -dx_deta/jacobian;
    deta_dx = -dy_dxi/jacobian;
    deta_dy = dx_dxi/jacobian;

    jacobian *= hx*hy;
    dw_dx = dw_dxi*dxi_dx + dw_deta*deta_dx;
    dw_dy = dw_dxi*dxi_dy + dw_deta*deta_dy;
}

void ElementMFSF::GenerateElement(const QuadraturePoint * const &quad1,
    const QuadraturePoint * const &quad2,
    const std::vector<double> &eCP,
    std::vector<double> &R,
    std::vector<double> &J)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};

    R.clear();
    J.clear();

    for (int ii = 0; ii < nqp1; ++ii)
    {
        std::vector<double> N1{};
        std::vector<double> dN1{};
        std::vector<double> N2{};
        std::vector<double> dN2{};
        GenerateBSplineBasis1D(qp1[ii], N1, N2, dN1, dN2);

        B1.insert(B1.end(), N1.begin(), N1.end());
        B2.insert(B2.end(), N2.begin(), N2.end());
        dB1.insert(dB1.end(), dN1.begin(), dN1.end());
        dB2.insert(dB2.end(), dN2.begin(), dN2.end());
    }

    for (int jj = 0; jj < nqp2; ++jj)
    {
        std::vector<double> N2{};
        std::vector<double> dN2{};
        N2.assign(B2.begin() + jj*(q+1), B2.begin() + (jj+1)*(q+1));
        dN2.assign(dB2.begin() + jj*(q+1), dB2.begin() + (jj+1)*(q+1));
        for (int ii = 0; ii < nqp1; ++ii)
        {
            std::vector<double> N1{};
            std::vector<double> dN1{};
            N1.assign(B1.begin() + ii*(p+1), B1.begin() + (ii+1)*(p+1));
            dN1.assign(dB1.begin() + ii*(p+1), dB1.begin() + (ii+1)*(p+1));

            std::vector<double> Rtemp{};
            double Jtemp = 0.0;
            GenerateElementSingleQP(eCP, N1, N2, dN1, dN2, Rtemp, Jtemp);
            
            R.insert(R.end(), Rtemp.begin(), Rtemp.end());
            J.push_back(Jtemp);
        }
    }
}

void ElementMFSF::GenerateElement(const QuadraturePoint * const &quad1,
    const QuadraturePoint * const &quad2,
    const std::vector<double> &eCP,
    std::vector<double> &B1, std::vector<double> &B2,
    std::vector<double> &dB1, std::vector<double> &dB2,
    std::vector<double> &W, std::vector<double> &J,
    std::vector<double> &dW_dx, std::vector<double> &dW_dy) const
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();

    B1.clear();
    B2.clear();
    dB1.clear();
    dB2.clear();
    W.clear();
    J.clear();
    dW_dx.clear();
    dW_dy.clear();

    for (int ii = 0; ii < nqp1; ++ii)
    {
        std::vector<double> N1{};
        std::vector<double> dN1{};
        std::vector<double> N2{};
        std::vector<double> dN2{};
        GenerateBSplineBasis1D(qp1[ii], N1, N2, dN1, dN2);

        B1.insert(B1.end(), N1.begin(), N1.end());
        B2.insert(B2.end(), N2.begin(), N2.end());
        dB1.insert(dB1.end(), dN1.begin(), dN1.end());
        dB2.insert(dB2.end(), dN2.begin(), dN2.end());
    }

    for (int jj = 0; jj < nqp2; ++jj)
    {
        std::vector<double> N2{};
        std::vector<double> dN2{};
        N2.assign(B2.begin() + jj*(q+1), B2.begin() + (jj+1)*(q+1));
        dN2.assign(dB2.begin() + jj*(q+1), dB2.begin() + (jj+1)*(q+1));
        for (int ii = 0; ii < nqp1; ++ii)
        {
            std::vector<double> N1{};
            std::vector<double> dN1{};
            N1.assign(B1.begin() + ii*(p+1), B1.begin() + (ii+1)*(p+1));
            dN1.assign(dB1.begin() + ii*(p+1), dB1.begin() + (ii+1)*(p+1));

            double w, jacobian, dw_dx, dw_dy;
            GenerateElementSingleQP(eCP,
                N1, N2, dN1, dN2,
                w, jacobian, dw_dx, dw_dy);

            W.push_back(w);
            J.push_back(jacobian);
            dW_dx.push_back(dw_dx);
            dW_dy.push_back(dw_dy);
        }
    }
}