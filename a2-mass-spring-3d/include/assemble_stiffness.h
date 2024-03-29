#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <d2V_spring_particle_particle_dq2.h>

/**
 * @param q - generalized coordinates for the mass-spring system
 * @param qdot - generalized velocity for the mass-spring system
 * @param V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
 * @param E - the mx2 spring connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
 * @param l0 - the mx1 vector of undeformed spring lengths
 * @param m - the mass of each particle in the mass-spring system
 * @param k - the stiffness of each spring in the mass-spring system
 * @param K - the 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function
 */
void assemble_stiffness(
    Eigen::SparseMatrixd &K,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> E,
    Eigen::Ref<const Eigen::VectorXd> l0,
    double k);