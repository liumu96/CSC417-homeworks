#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param M - dense 12x12 per-tetrahedron mass matrix
 * @param qdot - generalized velocity for the FEM system
 * @param element - the 1x4 vertex indices for this tetrahedron
 * @param density - density of material
 * @param volume - the undeformed tetrahedron volume
 */
void mass_matrix_linear_tetrahedron(
    Eigen::Matrix1212d &M,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double density,
    double volume);