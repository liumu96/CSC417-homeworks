#include "T_cloth.h"

// The kinetic energy of the whole cost mesh.

void T_cloth(
    double &T,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    Eigen::SparseMatrixd &M)
{
    T = 0.5 * qdot.transpose() * M * qdot;
}