#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * @param q - generalized coordinates of FEM system
 * @param V - vertex matrix for the mesh
 * @param element - vertex indices of the element
 * @param volume - volume of tetrahedron
 * @param C,D - material parameters
 * Output:
 * @param  energy - potential energy of tetrahedron
 */
void V_linear_tetrahedron(
    double &energy,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    double C,
    double D);