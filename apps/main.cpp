#include <stdio.h>
#include <thread>
#include <vector>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
using namespace std;
#include "app.h"
#include <utils/drawing_circle.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

//
#include "editing.h"
#include "playback.h"
inline std::tuple<vector<float>, int> normalized_weights(const App &app) {
  auto result = vector<float>(app.w.size());
  auto seed = -1;
  auto max_w = flt_min;
  auto sum = 0.f;
  for (auto i = 0; i < app.w.size(); ++i) {
    sum += app.w[i];
  }
  for (auto i = 0; i < app.w.size(); ++i) {
    if (app.w[i] > max_w) {
      seed = i;
      max_w = app.w[i];
    }
    result[i] = app.w[i] / sum;
  }
  return {result, seed};
}
vector<vec3f> concatenate_curve(const vector<vector<vec3f>> &curve) {
  auto result = vector<vec3f>{};
  for (auto i = 0; i < curve.size(); ++i)
    result.insert(result.end(), curve[i].begin(), curve[i].end());

  return result;
}
void set_common_uniforms(const App &app, const ogl_program *program) {
  auto &view = app.matrices.view;
  auto &projection = app.matrices.projection;
  set_uniform(program, "frame", identity4x4f);
  set_uniform(program, "view", view);
  set_uniform(program, "projection", projection);
  set_uniform(program, "eye", app.camera->frame.o);
  set_uniform(program, "envlight", (int)app.envlight);
  set_uniform(program, "gamma", app.shade_params.gamma);
  set_uniform(program, "exposure", app.shade_params.exposure);
  // set_uniform(program, "size", app.line_size);
  if (app.scene->environments.size()) {
    auto &env = app.scene->environments.front();
    if (env->envlight_diffuse)
      set_uniform(program, "envlight_irradiance", env->envlight_diffuse, 6);
    if (env->envlight_specular)
      set_uniform(program, "envlight_reflection", env->envlight_specular, 7);
    if (env->envlight_brdflut)
      set_uniform(program, "envlight_brdflut", env->envlight_brdflut, 8);
  }
}
vector<string> weigths_names(const App &app) {
  if (app.control_points.size() == 0)
    return {"Not yet available"};

  auto n = app.control_points.size();
  auto result = vector<string>(n);
  for (auto i = 0; i < n; ++i) {
    result[i] = "w" + std::to_string(i);
  }
  return result;
}
void draw_scene(const App &app, const vec4i &viewport) {
  clear_ogl_framebuffer(vec4f{0, 0, 0, 1});

  // Draw mesh and environment.
  draw_scene(app.scene, app.camera, viewport, app.shade_params);

  if (app.show_points) {
    auto program = &app.shaders.at("points");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "size", 3.0f * 0.0015f * app.line_size);

    set_uniform(program, "color", vec3f{0, 1, 0});
  }

  if (app.temp_levels > 0)
    draw_shape(app.temp_points[app.temp_levels]);
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(app.camera->film /
                             (camera_aspect * 2 * app.camera->lens)))
          : (2 * yocto::atan(app.camera->film / (2 * app.camera->lens)));
  auto view = frame_to_mat(inverse(app.camera->frame));
  auto projection = perspective_mat(
      camera_yfov, camera_aspect, app.shade_params.near, app.shade_params.far);

  if (app.gpu_shapes.find("edges") != app.gpu_shapes.end())
    gpu::draw_shape(app.gpu_shapes.at("edges"), app.gpu_shaders.at("points"),
                    gpu::Uniform("color", vec3f{0, 0, 0}));
  gpu::set_point_size(10);
  if (app.gpu_shapes.find("selected_points") != app.gpu_shapes.end()) {

    gpu::draw_shape(
        app.gpu_shapes.at("selected_points"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));
  }

  if (app.gpu_shapes.find("vector_field") != app.gpu_shapes.end())
    gpu::draw_shape(
        app.gpu_shapes.at("vector_field"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));

  if (app.gpu_shapes.find("vector_field_2") != app.gpu_shapes.end())
    gpu::draw_shape(
        app.gpu_shapes.at("vector_field_2"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{1, 0, 0}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));

  if (app.show_edges) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{0, 0, 0});
    draw_shape(&app.edges_shape);
  }
}

inline void sleep(int ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

inline bool is_pressing(gui_button button) {
  return button.state == gui_button::state::pressing;
}
inline bool is_releasing(gui_button button) {
  return button.state == gui_button::state::releasing;
}
inline bool is_down(gui_button button) {
  return button.state == gui_button::state::down ||
         button.state == gui_button::state::pressing;
}
inline bool is_pressing(const gui_input &input, gui_key key) {
  return is_pressing(input.key_buttons[(int)key]);
}

bool process_camera_move(App &app, const gui_input &input) {

  auto rotating = input.modifier_shift;
  auto panning = input.modifier_alt;
  auto &camera = *app.camera;

  auto update_camera_frame = [&](frame3f &frame, float &focus, bool rotating,
                                 bool panning, bool zooming) {
    auto last_pos = input.mouse_last;
    auto mouse_pos = input.mouse_pos;
    auto mouse_left = is_down(input.mouse_left);
    auto mouse_right = is_down(input.mouse_right);
    // handle mouse and keyboard for navigation
    if (mouse_left) {
      auto dolly = 0.0f;
      auto pan = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left)
          rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right)
          dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_left)
          pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    }
  };

  if (is_down(input.mouse_left) && (rotating || panning)) {
    update_camera_frame(camera.frame, app.camera_focus, rotating, panning,
                        false);
    return true;
  }

  // Zoom-in/out by scrolling;
  float zoom = input.scroll.y * 0.1;
  if (zoom != 0) {
    update_turntable(camera.frame, app.camera_focus, zero2f, zoom, zero2f);
    return true;
  }

  return false;
}
void process_key_input(App &app, const gui_input &input) {

  for (int i = 0; i < input.key_buttons.size(); i++) {
    auto key = gui_key(i);
    if (!is_pressing(input, key))
      continue;

    printf("%c pressed!\n", (char)key);

    if (key == gui_key::enter) {
      app.bezier_curve_shape = nullptr;
      app.shortest_geodesic_shape = nullptr;
      app.ground_truth_shape = nullptr;
      app.control_points.clear();
    }
  }
}
bool process_user_input(App &app, const gui_input &input) {
  //  static bool yyy = false;
  if (is_active(app.widget))
    return false;
  if (process_camera_move(app, input)) {
    update_camera_info(app, input);
    return false;
  }
  auto mouse = input.mouse_pos;
  auto size = vec2f{(float)input.window_size.x, (float)input.window_size.y};
  mouse = vec2f{2 * (mouse.x / size.x) - 1, 1 - 2 * (mouse.y / size.y)};
  auto editing = input.modifier_alt || input.modifier_shift;

  auto point = intersect_mesh(app, mouse);
  if (is_pressing(input.mouse_right))
    app.selected_point = intersect_control_points(app, mouse);

  if (is_releasing(input.mouse_left) && app.selected_point == -1 && !editing &&
      point.face != -1) {

    app.point_to_shape.insert(
        app.point_to_shape.end(),
        {(int)app.control_points.size(), (int)app.added_points.size()});
    app.control_points.push_back(point);
    app.w.push_back(1);
    app.point_moved.push_back(true);
    app.control_points_shape.push_back(
        add_points_shape(app, {point}, 0.0030, {0, 0, 0}));
  }
  auto drag = input.mouse_pos - input.mouse_last;
  if (is_down(input.mouse_right) && length(drag) > 0.05 &&
      app.selected_point != -1) {
    auto point = intersect_mesh(app, mouse);
    if (point.face != -1) {
      app.control_points[app.selected_point] = point;
      app.point_moved[app.selected_point] = true;
      update_points_shape(app.added_points, app.mesh, {point},
                          app.point_to_shape.at(app.selected_point));
      update_fields_from_lndmrks(app.mesh, app.topology, app.solver,
                                 app.control_points, app.type_of_solver,
                                 app.point_moved, app.fields);
      update_grads_from_lndmrks(app.operators, app.fields, app.point_moved,
                                app.grads);
      for (auto i = 0; i < app.point_moved.size(); ++i)
        app.point_moved[i] = false;

      if (app.bspline_curve_shape != nullptr) {
        update_rational_spline(
            app.mesh, app.topology, app.operators, app.dual_solver, app.fields,
            app.grads, app.control_points, app.w, app.subdiv,
            app.selected_point, app.polyline_pos, app.polyline);

        update_path_shape(app.bspline_curve_shape->instance->shape, app.mesh,
                          concatenate_curve(app.polyline_pos), 0.001);
      }
      if (app.bezier_curve_shape != nullptr) {
        auto pos = rational_bézier_curve(
            app.mesh, app.topology, app.operators, app.dual_solver, app.solver,
            app.fields, app.grads, app.control_points, app.w, app.subdiv);
        update_path_shape(app.bezier_curve_shape->instance->shape, app.mesh,
                          pos, 0.001);
      }
    }
  }
  if (is_releasing(input.mouse_right))
    app.selected_point = -1;
  return false;
}

void update_app(App &app, const gui_input &input) {
  // process_gui_input(app, input); TODO(giacomo)

  if (is_active(app.widget))
    return;

  app.window_size = input.window_size;
  process_key_input(app, input);
  process_user_input(app, input);

  auto tasks = vector<vec2i>{};
}

void draw(const gui_input &input, void *data) {
  auto &app = *(App *)data;
  app.started = true;
  auto mesh = &app.mesh;
  auto geometry = &app.topology;
  auto op = &app.operators;
  update_camera_info(app, input);

  // Do everything
  update_app(app, input);

  draw_scene(app, input.framebuffer_viewport);

  auto widget = app.widget;
  begin_widget(widget, "RCM on meshes");

  draw_bullet_text(widget, "Distance Field");
  vector<string> solver_names = {"Exact", "Graph Based"};
  if (draw_combobox(widget, "Solver", app.type_of_solver, solver_names)) {
    for (auto i = 0; i < app.point_moved.size(); ++i) {
      app.point_moved[i] = true;
    }
  }
  draw_checkbox(widget, "Re-Compute fields", app.recompute_fields);
  draw_separator(widget);
  draw_bullet_text(widget, "Parameters");
  draw_slider(widget, "Subdivision Steps", app.subdiv, 1, 12);

  if (draw_combobox(widget, "Weights", app.selected_weights,
                    weigths_names(app))) {
    if (app.w.size() > app.selected_weights) {
      app.curr_w = app.w[app.selected_weights];
    }
  }
  if (draw_slider(widget, "W", app.curr_w, 0.05, 5)) {
    if (app.w.size() > app.selected_weights) {
      app.w[app.selected_weights] = app.curr_w;
    }
    if (app.fields.size() > 0 && app.bezier_curve_shape != nullptr) {
      auto pos = rational_bézier_curve(*mesh, *geometry, *op, app.dual_solver,
                                       app.solver, app.fields, app.grads,
                                       app.control_points, app.w, app.subdiv);
      update_path_shape(app.bezier_curve_shape->instance->shape, app.mesh, pos,
                        0.001);
    }
    if (app.fields.size() > 0 && app.béz_interp_curve_shape != nullptr) {
      auto pos = béz_interp_rational_curve(
          *mesh, *geometry, *op, app.solver, app.dual_solver, app.fields,
          app.grads, app.control_points, app.w, app.subdiv);
      update_path_shape(app.béz_interp_curve_shape->instance->shape, app.mesh,
                        pos, 0.001);
    }
    if (app.fields.size() > 0 && app.bspline_curve_shape != nullptr) {

      std::tie(app.polyline_pos, app.polyline) = trace_rational_bspline(
          *mesh, *geometry, *op, app.dual_solver, app.fields, app.grads,
          app.control_points, app.w, app.subdiv);
      update_path_shape(app.bspline_curve_shape->instance->shape, app.mesh,
                        concatenate_curve(app.polyline_pos), 0.001);
    }
  }
  draw_separator(widget);
  draw_bullet_text(widget, "Curve Tracing");
  if (draw_button(widget, "Trace Circle")) {
    if (app.control_points.size() > 0 && app.control_points.size() > 0) {
      auto points = circle_control_points(app.mesh, app.topology,
                                          app.control_points[0], 0.125, 7);
      add_points_shape(app, points, 0.0030, zero3f);
      auto pos = vector<vec3f>{};
      auto curve = trace_circle(*mesh, *geometry, *op, app.solver,
                                app.dual_solver, app.control_points[0], 0.125,
                                7, app.subdiv, app.type_of_solver);
      for (auto i = 0; i < curve.size(); ++i) {
        pos.insert(pos.end(), curve[i].begin(), curve[i].end());
      }
      if (app.bezier_curve_shape == nullptr) {
        app.bezier_curve_shape = add_path_shape(app, pos, 0.001, {1, 0, 0});
      } else
        update_path_shape(app.bezier_curve_shape->instance->shape, app.mesh,
                          pos, 0.001);
    }
  }
  if (draw_button(widget, "Interpolant Bézier Curve")) {
    if (app.control_points.size() > 0 && app.control_points.size() > 2) {

      if (app.recompute_fields) {
        app.fields =
            fields_from_lndmrks(*mesh, *geometry, app.solver,
                                app.control_points, app.type_of_solver);
        app.grads = grads_from_lndmrks(app.operators, app.fields);
      } else {
        update_fields_from_lndmrks(*mesh, *geometry, app.solver,
                                   app.control_points, app.type_of_solver,
                                   app.point_moved, app.fields);
        update_grads_from_lndmrks(*op, app.fields, app.point_moved, app.grads);
        for (auto i = 0; i < app.point_moved.size(); ++i)
          app.point_moved[i] = false;
      }
      auto pos = béz_interp_rational_curve(
          *mesh, *geometry, *op, app.solver, app.dual_solver, app.fields,
          app.grads, app.control_points, app.w, app.subdiv);

      if (app.béz_interp_curve_shape == nullptr)
        app.béz_interp_curve_shape = add_path_shape(app, pos, 0.001, {1, 0, 0});
      else
        update_path_shape(app.béz_interp_curve_shape->instance->shape, *mesh,
                          pos, 0.001);
    }
  }

  if (draw_button(widget, "Rational Bézier Curve")) {
    if (app.control_points.size() > 0 && app.control_points.size() == 4) {
      auto control_points = vector<int>(app.control_points.size());
      for (auto i = 0; i < app.control_points.size(); ++i) {
        control_points[i] =
            forced_vert_from_point(app.mesh.triangles, app.control_points[i]);
        app.control_points[i] = point_from_vert(
            mesh->triangles, control_points[i], app.control_points[i].face);
        update_points_shape(app.control_points_shape[i]->instance->shape,
                            {app.mesh.positions[control_points[i]]}, 0.002);
      }

      if (app.recompute_fields) {
        app.fields =
            fields_from_lndmrks(*mesh, *geometry, app.solver,
                                app.control_points, app.type_of_solver);
        app.grads =
            // grads_from_lndmrks_xu(*mesh, *geometry, *op, app.fields);
            grads_from_lndmrks(app.operators, app.fields);
      } else {
        update_fields_from_lndmrks(*mesh, *geometry, app.solver,
                                   app.control_points, app.type_of_solver,
                                   app.point_moved, app.fields);
        update_grads_from_lndmrks(*op, app.fields, app.point_moved, app.grads);
        for (auto i = 0; i < app.point_moved.size(); ++i)
          app.point_moved[i] = false;
      }

      auto pos = vector<vec3f>{};

      pos = rational_bézier_curve(*mesh, *geometry, *op, app.dual_solver,
                                  app.solver, app.fields, app.grads,
                                  app.control_points, app.w, app.subdiv);
      if (app.bezier_curve_shape == nullptr)
        app.bezier_curve_shape = add_path_shape(app, pos, 0.001, {1, 0, 0});
      else
        update_path_shape(app.bezier_curve_shape->instance->shape, app.mesh,
                          pos, 0.001);
      if (app.bezier_curve_shape != nullptr) {
        clear_shape(app.bezier_curve_shape->instance->shape);
        app.bezier_curve_shape = nullptr;
      }
      app.bezier_curve_shape =
          add_path_shape(app, pos, 0.001, {1, 0, 0},
                         app.threshold_for_jumps * mesh->avg_edge_len);
    }
  }

  if (draw_button(widget, "Trace Rational B-spline")) {
    if (app.control_points.size() > 0 && app.control_points.size() >= 4) {
      auto control_points = vector<int>(app.control_points.size());
      for (auto i = 0; i < app.control_points.size(); ++i) {
        control_points[i] =
            forced_vert_from_point(app.mesh.triangles, app.control_points[i]);
        app.control_points[i] = point_from_vert(
            mesh->triangles, control_points[i], app.control_points[i].face);
        update_points_shape(app.control_points_shape[i]->instance->shape,
                            {app.mesh.positions[control_points[i]]}, 0.003);
      }

      if (app.recompute_fields) {
        app.fields =
            fields_from_lndmrks(*mesh, *geometry, app.solver,
                                app.control_points, app.type_of_solver);
        app.grads = grads_from_lndmrks(app.operators, app.fields);
      } else {
        update_fields_from_lndmrks(*mesh, *geometry, app.solver,
                                   app.control_points, app.type_of_solver,
                                   app.point_moved, app.fields);
        update_grads_from_lndmrks(*op, app.fields, app.point_moved, app.grads);
        for (auto i = 0; i < app.point_moved.size(); ++i)
          app.point_moved[i] = false;
      }

      std::tie(app.polyline_pos, app.polyline) = trace_rational_bspline(
          *mesh, *geometry, *op, app.dual_solver, app.fields, app.grads,
          app.control_points, app.w, app.subdiv);
      if (app.bspline_curve_shape != nullptr) {
        clear_shape(app.bspline_curve_shape->instance->shape);
        app.bspline_curve_shape = nullptr;
      }
      app.bspline_curve_shape = add_path_shape(
          app, concatenate_curve(app.polyline_pos), 0.001, {1, 0, 0});
    }
  }
  draw_separator(widget);
  draw_bullet_text(widget, "I/O");
  draw_textinput(widget, "Scene Name", app.models_name);
  if (draw_button(widget, "Save Landmarks")) {

    if (app.control_points.size() > 0) {
      save_control_points(app);
    }
  }
  if (draw_button(widget, "Load Landmarks")) {

    load_control_points(app);
  }
  if (draw_button(widget, "Save Curve")) {
    if (app.polyline.size() == 0)
      std::cout << "no curve to export" << std::endl;
    else
      save_curve(app);
  }
  if (draw_button(widget, "Load Curve")) {
    auto pos = load_curve(app);
    add_path_shape(app, pos, 0.001, {1, 0, 0});
  }
  draw_separator(widget);
  draw_checkbox(widget, "show edges", app.show_edges);
  draw_coloredit(widget, "Mesh Color", app.mesh_material->color);

  app.shade_params.faceted = app.show_edges;

  if (draw_button(widget, " Reset")) {
    app.control_points.clear();
    app.control_points.clear();
    app.control_points_shape.clear();
    for (auto path : app.added_paths)
      clear_shape(path->instance->shape);

    for (auto point : app.added_points)
      clear_shape(point->instance->shape);

    if (app.target_shape != nullptr)
      clear_shape(app.target_shape->instance->shape);
    if (app.source_shape != nullptr)
      clear_shape(app.source_shape->instance->shape);
    if (app.shortest_geodesic_shape != nullptr)
      clear_shape(app.shortest_geodesic_shape->instance->shape);
    if (app.ground_truth_shape != nullptr)
      clear_shape(app.ground_truth_shape->instance->shape);
    if (app.bezier_curve_shape != nullptr)
      clear_shape(app.bezier_curve_shape->instance->shape);
    if (app.bspline_curve_shape != nullptr)
      clear_shape(app.bspline_curve_shape->instance->shape);

    app.shortest_geodesic_shape = nullptr;
    app.ground_truth_shape = nullptr;
    app.bezier_curve_shape = nullptr;
    app.bspline_curve_shape = nullptr;
    app.vector_field.clear();
    app.added_paths.clear();
    app.field.clear();
    app.fields.clear();
    app.grads.clear();
    app.w.clear();
    app.point_moved.clear();
    app.point_to_shape.clear();
  }

  end_widget(widget);
}

int main(int num_args, const char *args[]) {
  auto app = App();

  int msaa = 1;
  string filename = args[1];

  if (!load_mesh(filename, app.mesh, app.topology, app.operators, app.solver,
                 app.dual_solver, app.error))
    print_fatal(app.error);
  init_bvh(app);

  app.filename = path_basename(filename);
  // Init window.
  auto win = new gui_window();
  win->msaa = msaa;
  init_window(win, {1080, 720}, "mesh viewer", true);
  win->user_data = &app;

  init_gpu(app, app.envlight);

  init_widget(app.widget, win);

  if (msaa > 1)
    set_ogl_msaa();

  run_ui(win, draw);

  // TODO(giacomo): delete app
  clear_window(win);
}
