#include "dV_linear_tetrahedron_dq.h"
#include "dphi_linear_tetrahedron_dX.h"
#include "dpsi_neo_hookean_dF.h"
#include "quadrature_single_point.h"

void dV_linear_tetrahedron_dq(
    Eigen::Vector12d &dV,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    double C,
    double D)
{
    auto neohookean_linear_tet = [&](
                                     Eigen::Vector12d &dV,
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

        Eigen::Vector9d dpsi;
        dpsi_neo_hookean_dF(dpsi, x * dphi, C, D);

        Eigen::MatrixXd B(9, 12);
        B.setZero();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                B.block<3, 3>(i * 3, j * 3) = Eigen::Matrix3d::Identity() * dphi(j, i);
            }
        }

        dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
}