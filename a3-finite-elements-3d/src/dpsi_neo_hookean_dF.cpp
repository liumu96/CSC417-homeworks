#include "dpsi_neo_hookean_dF.h"

void dpsi_neo_hookean_dF(
    Eigen::Vector9d &dw,
    Eigen::Ref<const Eigen::Matrix3d> F,
    double C,
    double D)
{
    Eigen::Vector9d dpsi;
    dpsi.setZero();
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    // Compute the Green deformation tensor
    Eigen::Matrix3d FTF = F.transpose() * F;

    // Compute the determinant of the deformation gradient
    double J = F.determinant();

    // Compute the trace of the Green deformation tensor
    double trFTF = FTF.trace();

    Eigen::Matrix3d adjF = F.inverse().transpose();

    // Compute the partial derivative of the trace term
    Eigen::Matrix3d d_trace_term =
        -2.0 / 3.0 * std::pow(J, -5.0 / 3.0) * trFTF * adjF + 2.0 * std::pow(J, -2.0 / 3.0) * F;

    // Compute the partial derivative of the determinant term
    Eigen::Matrix3d d_det_term = 2.0 * D * (J - 1) * adjF;

    // Combine the partial derivatives and multiply by C and D
    Eigen::Matrix3d term1 = C * d_trace_term;
    Eigen::Matrix3d term2 = 2.0 * D * (J - 1) * adjF;

    Eigen::Matrix3d res = term1 + term2;

    // Flatten the resulting matrix into a vector
    // dpsi << res(0, 0), res(1, 0), res(2, 0),
    //     res(0, 1), res(1, 1), res(2, 1),
    //     res(0, 2), res(1, 2), res(2, 2);
    dpsi << res(0, 0), res(0, 1), res(0, 2),
        res(1, 0), res(1, 1), res(1, 2),
        res(2, 0), res(2, 1), res(2, 2);

    dw.setZero();

    dw(0) = C * (F(0, 0) / pow(J, 2.0 / 3.0) * 2.0 - (F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) - D * (F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1)) * (-J + 1.0) * 2.0;
    dw(1) = C * (F(1, 0) / pow(J, 2.0 / 3.0) * 2.0 + (F(0, 1) * F(2, 2) - F(0, 2) * F(2, 1)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) + D * (F(0, 1) * F(2, 2) - F(0, 2) * F(2, 1)) * (-J + 1.0) * 2.0;
    dw(2) = C * (F(2, 0) / pow(J, 2.0 / 3.0) * 2.0 - (F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) - D * (F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1)) * (-J + 1.0) * 2.0;
    dw(3) = C * (F(0, 1) / pow(J, 2.0 / 3.0) * 2.0 + (F(1, 0) * F(2, 2) - F(1, 2) * F(2, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) + D * (F(1, 0) * F(2, 2) - F(1, 2) * F(2, 0)) * (-J + 1.0) * 2.0;
    dw(4) = C * (F(1, 1) / pow(J, 2.0 / 3.0) * 2.0 - (F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) - D * (F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0)) * (-J + 1.0) * 2.0;
    dw(5) = C * (F(2, 1) / pow(J, 2.0 / 3.0) * 2.0 + (F(0, 0) * F(1, 2) - F(0, 2) * F(1, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) + D * (F(0, 0) * F(1, 2) - F(0, 2) * F(1, 0)) * (-J + 1.0) * 2.0;
    dw(6) = C * (F(0, 2) / pow(J, 2.0 / 3.0) * 2.0 - (F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) - D * (F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0)) * (-J + 1.0) * 2.0;
    dw(7) = C * (F(1, 2) / pow(J, 2.0 / 3.0) * 2.0 + (F(0, 0) * F(2, 1) - F(0, 1) * F(2, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) + D * (F(0, 0) * F(2, 1) - F(0, 1) * F(2, 0)) * (-J + 1.0) * 2.0;
    dw(8) = C * (F(2, 2) / pow(J, 2.0 / 3.0) * 2.0 - (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0)) * 1.0 / pow(J, 5.0 / 3.0) * (trFTF) * (2.0 / 3.0)) - D * (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0)) * (-J + 1.0) * 2.0;
}