#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * @param qdot - generalied velocity of FEM system
 * @param element - vertex indices of the element
 * @param density - material density
 * @param volume - volume of tetrahedron
 * Output:
 * @param  T - kinetic energy of tetrahedron
 */
void T_linear_tetrahedron(
    double &T,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double density,
    double volume);