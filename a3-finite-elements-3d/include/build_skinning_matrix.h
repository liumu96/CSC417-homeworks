#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>

/**
 * @param N - the lxn sparse skinning matrix
 * @param V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position
 * @param T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
 * @param V_skin - lx3 matrix of vertices of the display mesh
 */
void build_skinning_matrix(
    Eigen::SparseMatrixd &N,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> T,
    Eigen::Ref<const Eigen::MatrixXd> V_skin);