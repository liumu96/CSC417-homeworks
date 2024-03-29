#pragma once

#define IMGUI_DEFINE_MATH_OPERATORS

#include <igl/unproject.h>
#include <pick_nearest_vertices.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <deque>
#include <Eigen/Dense>

namespace Visualize
{
    // custom phase space plot
    bool plot_phase_space(
        const char *label,
        ImVec2 q_bounds,
        ImVec2 q_dot_bounds,
        const Eigen::VectorXd &q,
        const Eigen::VectorXd &q_dot);

    void add_energy(float t, float T, float V);
    bool plot_energy(
        const char *label,
        unsigned int type,
        ImVec2 T_bounds,
        ImVec2 E_bounds,
        ImU32 plot_col);

    void setup(
        const Eigen::VectorXd &q,
        const Eigen::VectorXd &qdot,
        bool ps_plot = false);

    void add_object_to_scene(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::RowVector3d color);

    igl::opengl::glfw::Viewer &viewer();

    // animate geometry using physics simulation
    void rigid_transform_1d(unsigned int id, double x);

    void scale_x(unsigned int id, double x);

    void update_vertex_positions(
        unsigned int id,
        Eigen::Ref<const Eigen::VectorXd> pos);

    // UI methods
    bool mouse_down(
        igl::opengl::glfw::Viewer &viewer, int x, int y);

    bool mouse_up(
        igl::opengl::glfw::Viewer &viewer, int x, int y);

    bool mouse_move(
        igl::opengl::glfw::Viewer &viewer, int x, int y);

    const Eigen::Vector3d &mouse_world();

    const Eigen::Vector3d &mouse_drag_world();

    const std::vector<unsigned int> &picked_vertices();

    bool is_mouse_dragging();
}