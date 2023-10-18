#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0,
                        double k)
{
    typedef Eigen::Triplet<double> T;
    K.resize(q.size(), q.size());
    K.setZero();
    Eigen::Matrix66d spring_i_H;
    std::vector<T> K_entries;
    int q0_index, q1_index;
    for (int spring_i = 0; spring_i < E.rows(); spring_i++)
    {
        q0_index = 3 * E(spring_i, 0);
        q1_index = 3 * E(spring_i, 1);
        Eigen::Vector3d q0 = q.segment<3>(q0_index);
        Eigen::Vector3d q1 = q.segment<3>(q1_index);
        d2V_spring_particle_particle_dq2(spring_i_H, q0, q1, l0(spring_i), k);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                K_entries.push_back(T(q0_index + i, q0_index + j, -spring_i_H(i, j)));
                K_entries.push_back(T(q0_index + i, q1_index + j, -spring_i_H(i, j + 3)));
                K_entries.push_back(T(q1_index + i, q0_index + j, -spring_i_H(i + 3, j)));
                K_entries.push_back(T(q1_index + i, q1_index + j, -spring_i_H(i + 3, j + 3)));
            }
        }
    }
    K.setFromTriplets(K_entries.begin(), K_entries.end());
};