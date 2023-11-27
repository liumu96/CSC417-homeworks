#include "dV_cloth_gravity_dq.h"

// Gradient of potential energy due to gravity
void dV_cloth_gravity_dq(
    Eigen::VectorXd &fg,
    Eigen::SparseMatrixd &M,
    Eigen::Ref<const Eigen::Vector3d> g)
{
    fg.setZero();

    Eigen::VectorXd replicatedGravity = -g.replicate(M.rows() / g.size(), 1);

    fg = M * replicatedGravity;
}