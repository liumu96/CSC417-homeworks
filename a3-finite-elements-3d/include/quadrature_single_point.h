#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * Input:
 * q - generalized coordinates of the FEM system
 * element - vertex indices for the tetrahedron
 * volume - volume of the tetrahedron
 * integrand(out, q, X) - function to be integrated, returns value in out.
 * Output:
 * integrated - the value of the integrated function
 */
template <typename Ret, typename Integrand_Func>
inline void quadrature_single_point(
    Ret &&integrated,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::RowVectorXi> element,
    double volume,
    Integrand_Func integrand)
{
    // X provided in lambda, hence no need to provide it here

    // Call the integrand function to compute the value at the given coordinates
    integrand(integrated, q, element);

    // Multiply the result by the volume of the tetrahedron
    integrated *= volume;
}