#include "d2V_linear_tetrahedron_dq2.h"
#include "dphi_linear_tetrahedron_dX.h"
#include "d2psi_neo_hooken_dF2.h"
#include "quadrature_single_point.h"

void d2V_linear_tetrahedron_dq2(
    Eigen::Matrix1212d &H,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    double C,
    double D)
{
    // todo

    auto neohookean_linear_tet = [&](
                                     Eigen::Matrix1212d &dV,
                                     Eigen::Ref<const Eigen::VectorXd> q,
                                     Eigen::Ref<const Eigen::RowVectorXi> element)
    {
        Eigen::Vector3d X = V.row(element(0)).transpose();

        Eigen::Matrix34d x;
        x << q.segment<3>(element(0) * 3),
            q.segment<3>(element(1) * 3),
            q.segment<3>(element(2) * 3),
            q.segment<3>(element(3) * 3);

        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);

        Eigen::Matrix99d d2psi;
        d2psi_nei_hookean_dF2(d2psi, x * dphi, C, D);

        Eigen::MatrixXd B(9, 12);
        B.setZero();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                B.block<3, 3>(i * 3, j * 3) = Eigen::Matrix3d::Identity() * dphi(j, i);
            }
        }

        dV = B.transpose() * d2psi * B;
    };

    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);

    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();

    for (int i = 0; i < 12; i++)
    {
        if (es.eigenvalues()[i] < 1e-6)
        {
            DiagEval(i, i) = 1e-3;
        }
    }
    H = Evec * DiagEval * Evec.transpose();
}