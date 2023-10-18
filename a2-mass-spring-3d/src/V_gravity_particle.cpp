#include <V_gravity_particle.h>

// Compute the gravitational potental energy of a single mass particle.
void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q, double mass, Eigen::Ref<const Eigen::Vector3d> g)
{

    V = mass * std::abs(g.dot(q));
}