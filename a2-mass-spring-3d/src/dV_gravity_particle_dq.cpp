#include <dV_gravity_particle_dq.h>

// Compute the gradient of the gravitational potential energy for a single particle.
void dV_gravity_particle_dq(Eigen::Ref<Eigen::Vector3d> f, double mass, Eigen::Ref<const Eigen::Vector3d> g)
{
    f = mass * g;
}