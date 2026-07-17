#include "ElementMF.hpp"

void ElementMF::GenerateElementSingleQP(const double &xi, const double &eta,
    const std::vector<double> &eCP,
    std::vector<double> &R, std::vector<double> &dR_dx, std::vector<double> &dR_dy,
    double &jacobian) const
{
    BernsteinBasis * bern1 = new BernsteinBasis(p);
    BernsteinBasis * bern2 = new BernsteinBasis(q);

    RefElement * ref = new RefElement();

    std::vector<double> N1 = ref->GenerateBasis1DSingleQP(bern1, extraction1, xi);
    std::vector<double> dN1 = ref->GenerateBasisDerivative1DSingleQP(bern1, extraction1, xi, hx);
    std::vector<double> N2 = ref->GenerateBasis1DSingleQP(bern2, extraction2, eta);
    std::vector<double> dN2 = ref->GenerateBasisDerivative1DSingleQP(bern2, extraction2, eta, hy);

    delete bern1;
    delete bern2;
    delete ref;

    std::vector<double> N{};
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    double w = 0.0;
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

    R.clear();
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

    dR_dx.clear();
    dR_dy.clear();

    for (int j = 0; j<q+1; ++j)
    {
        for (int i = 0; i<p+1; ++i)
        {
            dR_dx.push_back(dxi_dx*dR_dxi[j*(p+1)+i] + deta_dx*dR_deta[j*(p+1)+i]);
            dR_dy.push_back(dxi_dy*dR_dxi[j*(p+1)+i] + deta_dy*dR_deta[j*(p+1)+i]);
        }
    }

    jacobian *= hx*hy;
}

void ElementMF::GenerateElement(const QuadraturePoint * const &quad1,
    const QuadraturePoint * const &quad2,
    const std::vector<double> &eCP)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();

    RR.clear();
    dRR_dx.clear();
    dRR_dy.clear();

    JxW.clear();

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            double jacobian;
            std::vector<double> R{};
            std::vector<double> dR_dx{};
            std::vector<double> dR_dy{};

            GenerateElementSingleQP(qp1[ii], qp2[jj], eCP, R, dR_dx, dR_dy, jacobian);

            double J_W = w1[ii]*w2[jj]*jacobian;

            for (int kk = 0; kk < (p+1)*(q+1); ++kk)
            {
                RR.push_back(R[kk]);
                dRR_dx.push_back(dR_dx[kk]);
                dRR_dy.push_back(dR_dy[kk]);
            }
            JxW.push_back(J_W);
        }
    }
}