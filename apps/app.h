#pragma once
#include <realtime/gpu.h>
#include <unordered_set>
#include <utils/mesh_io.h>
#include <utils/utilities.h>
#include <yocto/yocto_bvh.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

#define PROFILE 0

inline std::unordered_set<int> make_set(size_t n) {
  auto result = std::unordered_set<int>();
  result.reserve(n);
  for (int i = 0; i < n; i++) {
    result.insert(i);
  }
  return result;
}

struct Added_Path {
  geodesic_path path;
  vector<vec3f> positions;
  vec3f color;
  float radius;
  shade_instance *instance;
};

struct Added_Points {
  vector<mesh_point> points;
  vec3f color;
  float radius;
  shade_instance *instance;
};
enum struct editing_context : uint {
  is_doing_nothing,
  // is_rotating,
  is_editing_existing_curve,
  is_creating_new_curve,
  is_creating_first_curve,
};

static string context_names[4] = {
    "is_doing_nothing",
    "is_editing_existing_curve",
    "is_creating_new_curve",
    "is_creating_first_curve",
};

struct App {
  // TEMPORARY DATA
  string models_name = "test";
  vec3f prev_color = zero3f;
  vector<Added_Path *> added_paths = {};
  vector<Added_Points *> added_points = {};
  vector<Added_Path *> polyhedral_radial_geodesic_shape = {};
  vector<Added_Path *> smooth_radial_geodesic_shape = {};
  Added_Path *vector_field_shape = {};
  Added_Path *vector_field_shape2 = {};
  vector<Added_Path *> ground_truth_shapes = {};
  Added_Path *shortest_geodesic_shape = {};
  vector<Added_Path *> shortest_geodesic_shapes = {};
  Added_Path *straightest_geodesic_shape = {};
  Added_Path *polthier_straightest_shape = {};
  Added_Path *bezier_curve_shape = {};
  Added_Path *ground_truth_shape = {};
  Added_Path *mahen_curve_shape = {};
  Added_Path *bspline_curve_shape = {};
  Added_Points *source_shape = {};
  Added_Points *target_shape = {};
  Added_Points *samples_shape = {};
  vector<Added_Points *> control_points_shape = {};
  vec3f source_color = zero3f;
  vec3f cl_color = {1, 1, 0};
  int selected_point = -1;
  int selected_mesh = -1;
  int source_shape_entry = -1;
  int type_of_solver = graph;
  int subdiv = 6;
  float threshold_for_jumps = 10;
  vector<float> w = {};
  int selected_weights = 0;
  float curr_w = 1.0;
  bool show_gradient = false;
  vector<vector<mesh_point>> polyline = {};
  vector<vector<vec3f>> polyline_pos = {};
  vector<vector<int>> node_of_the_tree = {};
  vector<float> field = {};
  vector<vector<float>> fields = {};
  vector<vector<vec3f>> grads = {};
  float scaling_factor = 5;
  vector<vec3f> vector_field = {};
  vector<vec3f> glyph_normals{};
  vector<vec3f> vector_field2 = {};
  vector<ogl_shape *> temp_points = {};
  int temp_levels = -1;
  ogl_shape *eval_point_shape = new ogl_shape{};
  Svg svg = {};
  mesh_point xxx_point = {};
  vector<mesh_point> control_points = {};
  vector<bool> point_moved = {};
  unordered_map<int, int> point_to_shape;
  vector<int> verts = {};
  int selected_vert = 0;
  vec2i window_size = {};
  bool started = false;
  int playback_tick = 0;
  bool playback = false;

  shape_data mesh = {};
  shape_geometry topology = {};
  shape_op operators = {};
  geodesic_solver solver = {};
  dual_geodesic_solver dual_solver = {};

  struct {
    mat4f view = identity4x4f;
    mat4f projection = identity4x4f;
    mat4f projection_view = identity4x4f;
  } matrices;

  // bool recording_input = false;
  // vector<gui_input> input_record = {};
  // struct Editing_State {
  //   Gui_Input input = {};
  // };
  // editing_context context = editing_context::is_doing_nothing;
  // Editing_State state = {};

  // int editing_history_count = 0;
  // vector<Editing_State> editing_history = {};
  // int history_index = 0;

  // const Gui_Input &input() const { return state.input; }
  // Gui_Input &input() { return state.input; }

  // void commit_state() {
  //   editing_history.resize(history_index + 1);
  //   editing_history[history_index] = state;
  //   history_index += 1;
  //   printf("%s: %d\n", "commit_state()", history_index);
  // }

  string filename = "data/mesh.obj";
  string testname = "tests/test.json";
  string scene_name = "data/mesh.obj_gamma";
  string exported_scene_name = "scene.ply";
  float line_size = 1;
  float curve_size = 0.5;
  int scale_factor = 5;
  int vector_thickness_scale_factor = 3;
  float vector_thickness = 0.0001 * vector_thickness_scale_factor;
  int type_of_strip = 0;
  float vector_size = 0.001 * scale_factor;
  float lift_factor = 0;
  bool initial_strip = false;
  bool path_on_graph = false;
  bool recompute_fields = false;
  float time_of_last_click = -1;

  float angle = 0;

  string error = "";

  bool show_edges = false;
  bool show_points = true;
  bool envlight = false;
  bool show_control_polygon = false;

  ogl_shape edges_shape = {};
  ogl_shape branches_shape = {};
  ogl_shape co_branches_shape = {};

  gui_widget *widget = new gui_widget{};
  shade_scene *scene = nullptr;
  shade_material *spline_material = nullptr;
  shade_material *mesh_material = nullptr;
  shade_shape *mesh_shape = nullptr;
  vector<shade_material *> meshes_material = {};
  vector<shade_shape *> meshes_shape = {};
  shade_camera *camera = {};
  gpu::Camera gpu_camera;
  yocto::shade_params shade_params{};
  float camera_focus;
  shape_bvh bvh = {};

  // Data stored on the gpu for rendering.
  std::unordered_map<string, ogl_shape> shapes;
  std::unordered_map<string, ogl_program> shaders;

  std::unordered_map<string, gpu::Shape> gpu_shapes;
  std::unordered_map<string, gpu::Shader> gpu_shaders;
};

#include <thread>
template <typename F> inline void parallel_for(int size, F &&f) {
  auto num_threads = min(size, 16);
  auto threads = vector<std::thread>(num_threads);
  auto batch_size = (size + num_threads - 1) / num_threads;

  auto batch = [&](int k) {
    int from = k * batch_size;
    int to = min(from + batch_size, size);
    for (int i = from; i < to; i++)
      f(i);
  };

  for (int k = 0; k < num_threads; k++) {
    threads[k] = std::thread(batch, k);
  }
  for (int k = 0; k < num_threads; k++) {
    threads[k].join();
  }
}

template <typename F> inline void serial_for(int size, F &&f) {
  for (int i = 0; i < size; i++) {
    f(i);
  }
}

void init_bvh(App &app);
void init_bvhs(App &app);
shade_camera _make_framing_camera(const vector<vec3f> &positions);

void init_camera(App &app, const vec3f &from = vec3f{0, 0.5, 1.5},
                 const vec3f &to = {0, 0, 0});

vector<vec3f> make_normals(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions);

void init_bvh(App &app);

ray3f camera_ray(const App &app, vec2f mouse);

vec2f screenspace_from_worldspace(App &app, const vec3f &position);

mesh_point intersect_mesh(const App &app, vec2f mouse);

void init_gpu(App &app, bool envlight);

void delete_app(App &app);

bool load_program(ogl_program *program, const string &vertex_filename,
                  const string &fragment_filename);

void set_points_shape(ogl_shape *shape, const vector<vec3f> &positions);

void set_mesh_shape(ogl_shape *shape, const vector<vec3i> &triangles,
                    const vector<vec3f> &positions,
                    const vector<vec3f> &normals);
void set_polyline_shape(ogl_shape *shape, const vector<vec3f> &positions);

void save_curve(const App &app);
vector<vec3f> load_curve(App &app);
vector<vector<vec3f>> load_WA_curve(App &app);
void save_control_points(const App &app);
void load_control_points(App &app);

inline vec3f random_color() {
  auto result = zero3f;
  for (auto i = 0; i < 3; ++i) {
    result[i] = (float)rand() / RAND_MAX;
  }
  return result;
}

inline int intersect_control_points(App &app, const vec2f &mouse) {
  // Find index of clicked control point.
  float min_dist = flt_max;
  int selected_point = -1;
  float threshold = 0.1;

  for (int i = 0; i < app.control_points.size(); i++) {
    auto &mesh = app.mesh;
    // Skip handle points of non-selected anchors.
    auto point = app.control_points[i];
    auto pos = eval_position(mesh.triangles, mesh.positions, point);
    auto pos_ss = screenspace_from_worldspace(app, pos);
    float dist = length(pos_ss - mouse);
    if (dist < threshold && dist < min_dist) {
      selected_point = i;
      min_dist = dist;
    }
  }

  return selected_point;
}
inline void update_camera_info(App &app, const gui_input &input) {
  auto &camera = *app.camera;
  auto viewport = input.framebuffer_viewport;
  camera.aspect = (float)viewport.z / (float)viewport.w;

  auto camera_yfov =
      (camera.aspect >= 0)
          ? (2 * yocto::atan(camera.film / (camera.aspect * 2 * camera.lens)))
          : (2 * yocto::atan(camera.film / (2 * camera.lens)));

  app.matrices.view = frame_to_mat(inverse(camera.frame));
  app.matrices.projection = perspective_mat(
      camera_yfov, camera.aspect, app.shade_params.near, app.shade_params.far);
  app.matrices.projection_view = app.matrices.projection * app.matrices.view;
}

Added_Path *add_path_shape(App &app, const geodesic_path &path, float radius,
                           const vec3f &color);
Added_Path *add_path_shape(App &app, const vector<vec3f> &positions,
                           float radius, const vec3f &color,
                           const float &threshold = flt_max);
Added_Points *add_points_shape(App &app, const vector<mesh_point> &points,
                               float radius, const vec3f &color);

Added_Points *add_points_shape(App &app, const vector<vec3f> &points,
                               float radius, const vec3f &color);

Added_Path *add_vector_field_shape(App &app, const vector<vec3f> &vector_field,
                                   const float &scale, const float &radius,
                                   const vec3f &color);
Added_Path *add_glyph_shape(App &app, const vector<vec3f> &alpha,
                            const vector<vec3f> &instances,
                            const vector<vec3f> &normals, const float &scale,
                            const float &radius, const vec3f &color,
                            const float &offset);
void update_generic_vector_field_shape(shade_shape *shape,
                                       const vector<vec3f> &vector_field,
                                       const vector<vec3f> &instances,
                                       const vector<vec3f> &normals,
                                       const float &scale, const float &radius,
                                       const vec3f &color, const float &offset);
void update_glyph_shape(shade_shape *shape, const vector<vec3f> &alpha,
                        const vector<vec3f> &instances,
                        const vector<vec3f> &normals, const float &scale,
                        const float &radius, const vec3f &color,
                        const float &offset);

Added_Path *add_generic_vector_field_shape(
    App &app, const vector<vec3f> &vector_field, const vector<vec3f> &instances,
    const vector<vec3f> &normals, const float &scale, const float &radius,
    const vec3f &color, const float &offset);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3i> &adjacencies,
                       const geodesic_path &path, float radius);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius);
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius,
                       const float &treshold);
void update_points_shape(shade_shape *shape, const vector<vec3f> &positions,
                         float radius);
void update_points_shape(shade_shape *shape, const shape_data &mesh,
                         const vector<mesh_point> &points, float radius);
void update_points_shape(const vector<Added_Points *> &points,
                         const shape_data &data,
                         const vector<mesh_point> &new_points, const int entry);
void update_vector_field_shape(shade_shape *shape, shape_data &mesh,
                               const vector<vec3f> &vector_field,
                               const float &scale, const float &radius,
                               const vec3f &color);

void update_glpoints(App &app, const vector<vec3f> &positions,
                     const string &name);

void update_glpoints(App &app, const vector<mesh_point> &points,
                     const string &name = "selected_points");

void update_glvector_field(App &app, const vector<vec3f> &vector_field,
                           float &scale, const string &name = "vector_field");
