#include "assemble_stiffness.h"
#include <vector>

void assemble_stiffness(
    Eigen::SparseMatrixd &K,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::VectorXd> qdot,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> T,
    Eigen::Ref<const Eigen::VectorXd> v0,
    double C,
    double D)
{
    // todo
    K.resize(q.size(), q.size());
    K.setZero();
    std::vector<Eigen::Triplet<double>> K_entries;

    for (int i = 0; i < T.rows(); i++)
    {
        Eigen::RowVector4i current_tetrahedron = T.row(i);

        Eigen::Matrix1212d d2V;
        d2V_linear_tetrahedron_dq2(d2V, q, V, current_tetrahedron, v0(i), C, D);

        for (int vertex_i = 0; vertex_i < 4; vertex_i++)
        {
            for (int vertex_j = 0; vertex_j < 4; vertex_j++)
            {
                for (int xyz_i = 0; xyz_i < 3; xyz_i++)
                {
                    for (int xyz_j = 0; xyz_j < 3; xyz_j++)
                    {
                        K_entries.push_back(Eigen::Triplet<double>(
                            current_tetrahedron(vertex_i) * 3 + xyz_i,
                            current_tetrahedron(vertex_j) * 3 + xyz_j,
                            -d2V(vertex_i * 3 + xyz_i, vertex_j * 3 + xyz_j)));
                    }
                }
            }
        }
    }

    K.setFromTriplets(K_entries.begin(), K_entries.end());
}