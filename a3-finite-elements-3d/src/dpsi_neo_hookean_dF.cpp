#include "dpsi_neo_hookean_dF.h"

void dpsi_neo_hookean_dF(
    Eigen::Vector9d &psi,
    Eigen::Ref<const Eigen::Matrix3d> F,
    double C,
    double D)
{
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    // Compute the Green deformation tensor
    Eigen::Matrix3d FTF = F.transpose() * F;

    // Compute the determinant of the deformation gradient
    double J = F.determinant();

    // Compute the trace of the Green deformation tensor
    double trFTF = FTF.trace();

    // Compute the partial derivative of the trace term
    Eigen::Matrix3d d_trace_term = -2.0 / 3.0 * C * std::pow(J, -5.0 / 3.0) * trFTF * F.inverse().transpose() + 2.0 * C * std::pow(J, -2.0 / 3.0) * F;

    // Compute the partial derivative of the determinant term
    Eigen::Matrix3d d_det_term = 2.0 * D * (J - 1) * F.inverse().transpose();

    // Combine the partial derivatives and multiply by C and D
    Eigen::Matrix3d term1 = C * d_trace_term;
    Eigen::Matrix3d term2 = 2.0 * D * (J - 1) * F.inverse().transpose();

    Eigen::Matrix3d res = term1 + term2;

    // Flatten the resulting matrix into a vector
    psi << res(0, 0), res(1, 0), res(2, 0),
        res(0, 1), res(1, 1), res(2, 1),
        res(0, 2), res(1, 2), res(2, 2);
}