#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"
#include "rodrigues.h"

/**
 * Input:
 * @prams: q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles.
 * @prams:     The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
 * @prams: qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and
 * @prams:        the second 3 are the world space linear velocity.
 * @prams: dt - the integration time step
 * @prams: masses - a vector to mass matrices for each rigid body
 * @prams: forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
 * @prams:          while the second 3 doubles are the linear forces.
 * Output:
 * @prams: q - updated generalized coordinates
 * @prams: qdot - updated generalized velocities
 */
inline void exponential_euler(
    Eigen::VectorXd &q,
    Eigen::VectorXd &qdot,
    double dt,
    std::vector<Eigen::Matrix66d> &masses,
    Eigen::Ref<const Eigen::VectorXd> forces)
{
    // todo
}