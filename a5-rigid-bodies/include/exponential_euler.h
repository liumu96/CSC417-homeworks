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
    for (int i = 0; i < qdot.size() / 6; i++)
    {
        // Get all the current parameters
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * i).data());
        Eigen::Vector3d p = q.segment<3>(12 * i + 9);
        Eigen::Vector3d omega = qdot.segment<3>(6 * i);
        Eigen::Vector3d pdot = qdot.segment<3>(6 * i + 3);
        Eigen::Matrix3d I = masses[i].block<3, 3>(0, 0);
        double mass = masses[i](3, 3);
        Eigen::Vector3d t_torq = forces.segment<3>(6 * i);
        Eigen::Vector3d f_ext = forces.segment<3>(6 * i + 3);

        Eigen::Matrix3d update_rotation, R_next_t;
        rodrigues(update_rotation, omega * dt);
        R_next_t = update_rotation * R;
        q.segment<9>(12 * i) = Eigen::Map<const Eigen::Vector9d>(R_next_t.data());
        q.segment<3>(12 * i + 9) = p + dt * pdot;

        qdot.segment<3>(6 * i) = (R * I * R.transpose()).inverse() * ((R * I * R.transpose()) * omega - dt * omega.cross((R * I * R.transpose()) * omega) + dt * t_torq);
        qdot.segment<3>(6 * i + 3) = (mass * pdot + dt * f_ext) / mass;
    }
}