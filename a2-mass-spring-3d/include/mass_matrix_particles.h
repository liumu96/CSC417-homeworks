#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void mass_matrix_particles(
    Eigen::SparseMatrixd &M,
    Eigen::Ref<const Eigen::VectorXd> q,
    double mass);