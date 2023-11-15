#include "visualization.h"

namespace Visualize
{
    igl::opengl::glfw::Viewer g_viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;

    // meshes in the scene
    std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> g_geometry;
    std::vector<unsigned int> g_id; // id into libigl for these meshes

    // pointers to q and qdot
    Eigen::VectorXd const *g_q;
    Eigen::VectorXd const *g_qdot;

    // cache for phase space data
    std::deque<std::pair<float, float>> g_state;
}

igl::opengl::glfw::Viewer &Visualize::viewer()
{
    return g_viewer;
}

bool Visualize::plot_phase_space(
    const char *label,
    ImVec2 q_bounds,
    ImVec2 q_dot_bounds,
    const Eigen::VectorXd &q,
    const Eigen::VectorXd &q_dot)
{
    using namespace ImGui;

    unsigned int cache_size = 10000;
    unsigned int num_lines = 5;

    ImGuiContext &g = *GImGui;
    const ImGuiStyle &style = g.Style;

    if (g_state.size() > cache_size)
    {
        g_state.pop_front();
    }

    // update plotting cache
    g_state.push_back(std::make_pair(q(0), q_dot(0)));

    const ImGuiStyle &Style = GetStyle();
    const ImGuiIO &IO = GetIO();
    ImDrawList *DrawList = GetWindowDrawList();
    ImGuiWindow *Window = GetCurrentWindow();

    if (Window->SkipItems)
        return false;

    // header ans spacing
    int hovered = IsItemActive() || IsItemHovered();
    Dummy(ImVec2(0, 3));

    // prepare canvas
    ImVec2 avail = GetContentRegionAvail();
    ImVec2 Canvas(ImMin(avail.x, avail.y), ImMin(avail.x, avail.y));

    Canvas = CalcItemSize(Canvas, style.FramePadding.x * 2.0f, style.FramePadding.y * 2.0f);
    ImRect bb(Window->DC.CursorPos, Window->DC.CursorPos + Canvas);

    const ImGuiID id = Window->GetID(label);

    RenderFrame(bb.Min, bb.Max, GetColorU32(ImGuiCol_FrameBg, 1), true, Style.FrameRounding);

    // local grid coordinates are -1, 1 in both directions
    auto pix_to_normalized = [&bb, &Canvas](ImVec2 pix)
    { return ImVec2((pix.x - bb.Min.x) / Canvas.x, (pix.y - bb.Min.y) / Canvas.y); };
    auto normalized_to_pix = [&bb, &Canvas](ImVec2 norm)
    { return ImVec2(norm.x * Canvas.x + bb.Min.x, norm.y * Canvas.y + bb.Min.y); };
    auto data_to_normalized = [&q_bounds, &q_dot_bounds](ImVec2 state)
    { return ImVec2((state.x - q_bounds.x) / (q_bounds.y - q_bounds.x), (state.y - q_dot_bounds.x) / (q_dot_bounds.y - q_dot_bounds.x)); };

    // background grid centered on origin
    for (float i = 0.f; i <= 1.f; i += 1.f / static_cast<float>(num_lines - 1))
    {
        DrawList->AddLine(
            normalized_to_pix(ImVec2(i, 0.f)),
            normalized_to_pix(ImVec2(i, 1.f)),
            GetColorU32(ImGuiCol_TextDisabled), 1.2);
    }

    for (float i = 0.f; i <= 1.f; i += 1.f / static_cast<float>(num_lines - 1))
    {
        DrawList->AddLine(
            normalized_to_pix(ImVec2(0.f, i)),
            normalized_to_pix(ImVec2(1.f, i)),
            GetColorU32(ImGuiCol_TextDisabled), 1.2);
    }

    // plot phase space trajectory
    bool clip_p1;
    bool clip_p2;
    for (unsigned int i = 0; i < g_state.size() - 1; ++i)
    {

        clip_p1 = false;
        clip_p2 = false;

        ImVec2 p1 = data_to_normalized(ImVec2(g_state[i].first, g_state[i].second));
        ImVec2 p2 = data_to_normalized(ImVec2(g_state[i + 1].first, g_state[i + 1].second));

        if (p1.x < 0.f || p1.x > 1.f || p1.y < 0.f || p1.y > 1.f)
        {
            clip_p1 = true;
        }

        if (p2.x < 0.f || p2.x > 1.f || p2.y < 0.f || p2.y > 1.f)
        {
            clip_p2 = true;
        }

        p1.x = ImMin(ImMax(p1.x, 0.f), 1.f);
        p1.y = ImMin(ImMax(p1.y, 0.f), 1.f);
        p2.x = ImMin(ImMax(p2.x, 0.f), 1.f);
        p2.y = ImMin(ImMax(p2.y, 0.f), 1.f);

        if (!clip_p1 || !clip_p2)
        {
            DrawList->AddLine(
                normalized_to_pix(p1),
                normalized_to_pix(p2),
                4290733594, 2);
        }
    }

    // label axes

    return true;
}

void Visualize::setup(const Eigen::VectorXd &q,
                      const Eigen::VectorXd &qdot)
{
    g_q = &q;
    g_qdot = &qdot;

    // add new menu for phase space plotting
    Visualize::g_viewer.plugins.push_back(&plugin);
    plugin.widgets.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]()
    {
        // menu.draw_viewer_menu();
    };

    // Add content to the default menu window
    menu.callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        // ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin(
            "Pahse Space Plot", nullptr,
            ImGuiWindowFlags_NoSavedSettings);

        ImVec2 min = ImGui::GetWindowContentRegionMin();
        ImVec2 max = ImGui::GetWindowContentRegionMax();

        max.x = (max.x - min.x) / 2;
        max.y -= min.y + ImGui::GetTextLineHeightWithSpacing() * 3;

        Visualize::plot_phase_space("example", ImVec2(-15, 15), ImVec2(-15, 15), *g_q, *g_qdot);

        ImGui::End();
    };

    Visualize::g_viewer.core().background_color.setConstant(1.0);
    Visualize::g_viewer.core().is_animating = false;
}

void Visualize::add_object_to_scene(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::RowVector3d color)
{
    // add mesh to libigl and store id for access later
    if (g_geometry.size() == 0)
    {
        g_id.push_back(0);
    }
    else
    {
        g_id.push_back(g_viewer.append_mesh());
    }

    g_viewer.data().set_mesh(V, F);
    g_viewer.data().set_colors(color);

    // add mesh to geometry vector
    g_geometry.push_back(std::make_pair(V, F));
}

void Visualize::rigid_transform_1d(unsigned int id, double x)
{
    // reset vertex positions
    for (unsigned int ii = 0; ii < g_geometry[id].first.rows(); ii++)
    {
        g_viewer.data_list[g_id[id]].V(ii, 0) = g_geometry[id].first(ii, 0) + x;
    }

    // tell viewer to update
    g_viewer.data_list[g_id[id]].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
}

void Visualize::scale_x(unsigned int id, double x)
{
    // reset vertex positions
    for (unsigned int ii = 0; ii < g_geometry[id].first.rows(); ++ii)
    {
        g_viewer.data_list[g_id[id]].V(ii, 0) = x * g_geometry[id].first(ii, 0);
    }

    // tell viewer to update
    g_viewer.data_list[g_id[id]].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
}