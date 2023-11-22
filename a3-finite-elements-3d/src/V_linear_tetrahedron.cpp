#include "V_linear_tetrahedron.h"
#include "dphi_linear_tetrahedron_dX.h"
#include "psi_neo_hookean.h"
#include "quadrature_single_point.h"

void V_linear_tetrahedron(
    double &energy,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    double C,
    double D)
{
    auto neohookean_linear_tet = [&](
                                     double &e,
                                     Eigen::Ref<const Eigen::VectorXd> q,
                                     Eigen::Ref<const Eigen::RowVectorXi> element)
    {
        Eigen::Vector3d X = V.row(element(0)).transpose();

        Eigen::Matrix34d x;
        Eigen::Matrix43d dphi;

        x << q.segment<3>(element(0) * 3),
            q.segment<3>(element(1) * 3),
            q.segment<3>(element(2) * 3),
            q.segment<3>(element(3) * 3);
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        psi_neo_hookean(e, x * dphi, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);
}