#include "psi_neo_hookean.h"

void psi_neo_hookean(
    double &psi,
    Eigen::Ref<const Eigen::Matrix3d> F,
    double C,
    double D)
{
    double J = F.determinant();
    psi = C * (pow(J, -2.0 / 3.0) * (F.transpose() * F).trace() - 3) + D * pow(J - 1, 2.0);
}