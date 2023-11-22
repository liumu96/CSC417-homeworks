#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <dV_linear_tetrahedron_dq.h>

/**
 * @param f - the vector 3xn vector of forces acting on each node of the mass-spring system
 * @param q - generalized coordinates for the FEM system
 * @param qdot - generalized velocity for the FEM system
 * @param V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position
 * @param T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
 * @param v0 - the mx1 vector of undeformed tetrahedron volumes
 * @param C, D - material parameters for the Neo-Hookean model
 */
void assemble_forces(
    Eigen::VectorXd &f,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> T,
    Eigen::Ref<const Eigen::VectorXd> v0,
    double C,
    double D);