#include "d2V_spring_particle_particle_dq2.h"

#include <iostream>
void d2V_spring_particle_particle_dq2(
    Eigen::Ref<Eigen::Matrix66d> H,
    Eigen::Ref<const Eigen::Vector3d> q0,
    Eigen::Ref<const Eigen::Vector3d> q1,
    double l0,
    double stiffness)
{
    H.setZero();

    Eigen::Vector3d delta_q = q1 - q0;
    double norm_delta_q = delta_q.norm();
    double norm_delta_q_cubic = norm_delta_q * norm_delta_q * norm_delta_q;

    for (int i = 0; i < 6; i++)
    {
        for (int j = i; j < 6; j++)
        {

            H(i, j) = delta_q[i % 3] * delta_q[j % 3] * l0 / norm_delta_q_cubic;
            if (i == j)
            {
                H(i, j) += (1 - l0 / norm_delta_q);
            }
            else if (i / 3 < 1 && j / 3 > 0)
            {
                H(i, j) = -H(i, j);
            }
            if (i != j)
            {
                H(j, i) = H(i, j);
            }
        }
    }

    H *= stiffness;
}
