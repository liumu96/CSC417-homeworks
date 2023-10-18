#include <mass_matrix_particles.h>

// Compute the sparse, diagonal mass matrix that stores the mass of each particle in the mass-spring on its diagonal.
void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass)
{
    M.resize(q.size(), q.size());
    M.setIdentity();
    M = mass * M;
}
