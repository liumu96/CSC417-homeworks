#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0,
                     double mass, double k)
{
    f.resize(q.size());
    f.setZero();
    // Since mass is uniform for all particles, gravitational force will be uniform too
    Eigen::Vector3d g(0, 0, 0);
    Eigen::Vector3d gf;
    dV_gravity_particle_dq(gf, mass, g);
    for (int spring_i = 0; spring_i < E.rows(); spring_i++)
    {
        Eigen::Vector3d q0 = q.segment<3>(3 * E(spring_i, 0));
        Eigen::Vector3d q1 = q.segment<3>(3 * E(spring_i, 1));
        Eigen::Vector6d f_spring_i;
        dV_spring_particle_particle_dq(f_spring_i, q0, q1, l0(spring_i), k);
        f.segment<3>(3 * E(spring_i, 0)) -= f_spring_i.head<3>() + gf;
        f.segment<3>(3 * E(spring_i, 1)) -= f_spring_i.tail<3>() + gf;
    }
};