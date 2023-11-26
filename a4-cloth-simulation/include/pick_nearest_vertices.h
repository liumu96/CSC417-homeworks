#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"

/**
 * Input:
 * @param: win - window coordinate of mouse click (x_window, y_window, 0)
 * @param: view - view transformation matrix
 * @param: proj - projection matrix
 * @param: viewport - viewport coordinates .
 * @param: V - 3xn dense matrix of mesh vertices, each row of the matrix is a single vertex.
 * @param: radius - selection radius for vertex picking
 * Output:
 * @param: verts - vertex ids (rows in V) of selected vertices
 */

bool pick_nearest_vertices(
    std::vector<unsigned int> &verts,
    Eigen::Ref<const Eigen::Vector3d> win,
    Eigen::Ref<const Eigen::Matrix44f> view,
    Eigen::Ref<const Eigen::Matrix44f> proj,
    Eigen::Vector4f viewport,
    Eigen::Ref<const Eigen::MatrixXd> V,
    Eigen::Ref<const Eigen::MatrixXi> F,
    double radius);