#include "app.h"

#include <utils/logging.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_modelio.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_shape.h>

void shift(vec3i &v, const int offset) {
  v.x += offset;
  v.y += offset;
  v.z += offset;
}
vector<vec3i> shift(const vector<vec3i> &v, const int offset) {
  auto result = v;
  for (auto &x : result) {
    shift(x, offset);
  }

  return result;
}

// vector<vec4f> bake_ambient_occlusion(App &app, int num_samples) {
//   auto result = vector<vec4f>(app.mesh.positions.size(), zero4f);
//   auto rng = rng_state{};

//   auto f = [&](int i) {
//     auto frame = basis_fromz(app.mesh.normals[i]);
//     for (int sample = 0; sample < num_samples; sample++) {
//       auto dir = sample_hemisphere_cos(rand2f(rng));
//       dir = transform_direction(frame, dir);
//       auto ray = ray3f{app.mesh.positions[i], dir};
//       auto isec = intersect_triangles_bvh(app.bvh, app.mesh.triangles,
//                                           app.mesh.positions, ray);
//       if (!isec.hit && dir.y > 0) {
//         result[i] += vec4f{1, 1, 1, 1};
//       } else {
//         result[i] += vec4f{0, 0, 0, 1};
//       }
//     }
//   };
//   parallel_for((int)app.mesh.positions.size(), f);

//   for (auto &r : result)
//     r /= r.w;
//   return result;
// }

// TODO(giacomo): rename after removing realtime lib
// TODO(giacomo): consider adding these functions to yocto_gui
shade_camera _make_lookat_camera(const vec3f &from, const vec3f &to,
                                 const vec3f &up = {0, 1, 0}) {
  auto camera = shade_camera{};
  camera.frame = lookat_frame(from, to, {0, 1, 0});
  camera.focus = length(from - to);
  return camera;
}

shade_camera _make_framing_camera(const vector<vec3f> &positions) {
  auto direction = vec3f{0, 1, 2};
  auto box = bbox3f{};
  for (auto &p : positions) {
    expand(box, p);
  }
  auto box_center = center(box);
  auto box_size = max(size(box));
  return _make_lookat_camera(direction * box_size + box_center, box_center);
}

void init_camera(App &app, const vec3f &from, const vec3f &to) {

  *app.camera = _make_framing_camera(app.mesh.positions);

  app.camera_focus = app.camera->focus;
}

void set_points_shape(ogl_shape *shape, const vector<vec3f> &positions) {
  auto sphere = make_uvsphere({16, 16}, 1, {1, 1});

  set_vertex_buffer(shape, sphere.positions, 0);
  set_index_buffer(shape, quads_to_triangles(sphere.quads));

  set_vertex_buffer(shape, positions, 1);
  set_instance_buffer(shape, 1, true);
}
// void set_points_shape(ogl_shape *shape, const bezier_mesh &mesh,
//                       const vector<mesh_point> &points) {
//   auto pos = vector<vec3f>(points.size());
//   for (int i = 0; i < points.size(); i++) {
//     pos[i] = eval_position(mesh, points[i]);
//   }
//   set_points_shape(shape, pos);
// }

void set_polyline_shape(ogl_shape *shape, const vector<vec3f> &positions) {
  if (positions.empty())
    return;
  auto cylinder = make_uvcylinder({16, 1, 1}, {1, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_vertex_buffer(shape, cylinder.positions, 0);
  set_vertex_buffer(shape, cylinder.normals, 1);
  set_index_buffer(shape, quads_to_triangles(cylinder.quads));

  auto froms = vector<vec3f>();
  auto tos = vector<vec3f>();
  auto colors = vector<vec3f>();
  froms.reserve(positions.size() - 1);
  tos.reserve(positions.size() - 1);
  // colors.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    if (positions[i] == positions[i + 1])
      continue;
    froms.push_back(positions[i]);
    tos.push_back(positions[i + 1]);
    // colors.push_back({sinf(i) * 0.5f + 0.5f, cosf(i) * 0.5f + 0.5f, 0});
  }

  // TODO(giacomo): solve rendering bug with degenerate polyline
  if (froms.empty()) {
    shape->num_instances = 0;
    set_vertex_buffer(shape, {}, 0);
    set_vertex_buffer(shape, {}, 1);
    set_vertex_buffer(shape, {}, 2);
    set_vertex_buffer(shape, {}, 3);
    set_index_buffer(shape, vector<vec3i>{});
  } else {
    set_vertex_buffer(shape, froms, 2);
    set_instance_buffer(shape, 2, true);
    set_vertex_buffer(shape, tos, 3);
    set_instance_buffer(shape, 3, true);
  }
}

void set_mesh_shape(ogl_shape *shape, const vector<vec3i> &triangles,
                    const vector<vec3f> &positions,
                    const vector<vec3f> &normals) {
  set_vertex_buffer(shape, positions, 0);
  set_vertex_buffer(shape, normals, 1);
  set_index_buffer(shape, triangles);
}

void _set_polyline_shape(ogl_shape *shape, const vector<vec3f> &positions,
                         const vector<vec3f> &normals) {
  set_vertex_buffer(shape, positions, 0);
  if (normals.size()) {
    set_vertex_buffer(shape, normals, 1);
  }
  shape->elements = ogl_element_type::lines;
}

vector<vec3f> make_normals(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto &normal : normals)
    normal = zero3f;
  for (auto &t : triangles) {
    auto normal =
        cross(positions[t.y] - positions[t.x], positions[t.z] - positions[t.x]);
    normals[t.x] += normal;
    normals[t.y] += normal;
    normals[t.z] += normal;
  }
  for (auto &normal : normals)
    normal = normalize(normal);
  return normals;
}

void init_bvh(App &app) {
  app.bvh = make_triangles_bvh(app.mesh.triangles, app.mesh.positions, {});
}

ray3f camera_ray(const App &app, vec2f mouse) {
  auto camera_ray = [](const frame3f &frame, float lens, const vec2f &film,
                       const vec2f &image_uv) {
    auto e = zero3f;
    auto q =
        vec3f{film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), lens};
    auto q1 = -q;
    auto d = normalize(q1 - e);
    auto ray = ray3f{transform_point(frame, e), transform_direction(frame, d)};
    return ray;
  };

  mouse += 1;
  mouse /= 2;
  mouse.y = 1 - mouse.y;
  auto &camera = *app.camera;
  return camera_ray(camera.frame, camera.lens,
                    {camera.film, camera.film / camera.aspect}, mouse);
}

vec2f screenspace_from_worldspace(App &app, const vec3f &position) {
  auto [x, y, z] = position;
  auto uv4f = app.matrices.projection_view * vec4f{x, y, z, 1};
  return vec2f{uv4f.x / uv4f.w, uv4f.y / uv4f.w};
};

mesh_point intersect_mesh(const App &app, vec2f mouse) {
  auto ray = camera_ray(app, mouse);
  auto isec = intersect_triangles_bvh(app.bvh, app.mesh.triangles,
                                      app.mesh.positions, ray);

  if (isec.hit) {
    return mesh_point{isec.element, isec.uv};
  } else {
    return {-1, {0, 0}};
  }
}

bool load_program(ogl_program *program, const string &vertex_filename,
                  const string &fragment_filename) {
  auto error = ""s;
  auto vertex_source = ""s;
  auto fragment_source = ""s;

  if (!load_text(vertex_filename, vertex_source, error)) {
    printf("error loading vertex shader (%s): \n%s\n", vertex_filename.c_str(),
           error.c_str());
    return false;
  }
  if (!load_text(fragment_filename, fragment_source, error)) {
    printf("error loading fragment shader (%s): \n%s\n",
           fragment_filename.c_str(), error.c_str());
    return false;
  }

  auto error_buf = ""s;
  if (!set_program(program, vertex_source, fragment_source, error, error_buf)) {
    printf("\nerror: %s\n", error.c_str());
    printf("    %s\n", error_buf.c_str());
    return false;
  }
  return true;
}

void init_gpu(App &app, bool envlight) {
  string error;
  init_ogl(error);

  app.scene = new shade_scene{};

  app.camera = add_camera(app.scene);
  // app.scene->cameras.push_back(&app.camera);
  app.spline_material =
      add_material(app.scene, {0, 0, 0}, {1, 0, 0}, 1, 0, 0.4);

  init_camera(app);

  // Init opengl mesh

  auto &mesh = app.mesh;
  // auto  colors     = bake_ambient_occlusion(app, 1024);
  app.mesh_shape = add_shape(app.scene, {}, {}, mesh.triangles, {},
                             mesh.positions, mesh.normals, {}, {});
  app.mesh_material =
      add_material(app.scene, {0, 0, 0}, {0.9, 0.9, 0.9}, 0.04, 0, 0.5);
  add_instance(app.scene, identity3x4f, app.mesh_shape, app.mesh_material);

  app.shade_params.hide_environment = true;
  app.shade_params.exposure = -0.5;
  app.shade_params.background = {1, 1, 1, 1};
  app.shade_params.lighting =
      envlight ? shade_lighting_type::envlight : shade_lighting_type::eyelight;

  init_scene(app.scene, true);

  // Init shaders.
  auto base = string(SHADERS_PATH); // defined in parent CMakeLists.txt

  load_program(&app.shaders["mesh"], base + "mesh.vert", base + "mesh.frag");
  load_program(&app.shaders["lines"], base + "points.vert",
               base + "points.frag");

  if (envlight) {
    load_program(&app.shaders["points"], base + "sphere.vert",
                 base + "sphere-envlight.frag");
    load_program(&app.shaders["polyline"], base + "polyline.vert",
                 base + "polyline-envlight.frag");
  } else {
    load_program(&app.shaders["points"], base + "sphere.vert",
                 base + "sphere-eyelight.frag");
    load_program(&app.shaders["polyline"], base + "polyline.vert",
                 base + "polyline-eyelight.frag");
  }
  load_program(&app.shaders["flat"], base + "flat.vert", base + "flat.frag");

  app.gpu_shaders["points"] =
      gpu::make_shader_from_file(base + "points.vert", base + "points.frag");
  app.gpu_shaders["mesh"] =
      gpu::make_shader_from_file(base + "mesh.vert", base + "mesh.frag");
  // init edge shape
  auto surface_offset = 0.00002f;
  auto positions = app.mesh.positions;
  auto normals = app.mesh.normals;
  auto triangles = app.mesh.triangles;
  for (int i = 0; i < positions.size(); i++) {
    positions[i] += normals[i] * surface_offset;
  }
  auto edges = vector<vec2i>();
  edges.reserve(positions.size() * 3);
  for (auto &t : triangles) {
    if (t.x < t.y) {
      edges.push_back({t.x, t.y});
      app.mesh.avg_edge_len += length(positions[t.x] - positions[t.y]);
    }
    if (t.y < t.z) {
      edges.push_back({t.y, t.z});
      app.mesh.avg_edge_len += length(positions[t.y] - positions[t.z]);
    }
    if (t.z < t.x) {
      edges.push_back({t.z, t.x});
      app.mesh.avg_edge_len += length(positions[t.z] - positions[t.x]);
    }
  }

  app.mesh.avg_edge_len /= edges.size();

  set_vertex_buffer(&app.edges_shape, positions, 0);
  set_index_buffer(&app.edges_shape, edges);
}

void delete_app(App &app) {
  for (auto &[name, shape] : app.gpu_shapes) {
    delete_shape(shape);
  }
  for (auto &[name, shader] : app.gpu_shaders) {
    gpu::delete_shader(shader);
  }
}

#include <yocto/yocto_sceneio.h>

//#include "scene_exporter.h"
#ifdef _WIN32
#undef near
#undef far
#endif

Added_Path *add_vector_field_shape(App &app, const vector<vec3f> &vector_field,
                                   const float &scale, const float &radius,
                                   const vec3f &color) {
  auto added_path = app.added_paths.emplace_back(new Added_Path{});
  added_path->color = color;
  added_path->radius = radius;
  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  auto frames = vector<frame3f>();
  frames.reserve(vector_field.size());
  if (vector_field.size() == app.mesh.positions.size()) {
    for (int i = 0; i < app.mesh.positions.size(); i++) {
      auto from = app.mesh.positions[i];
      auto to = from + scale * vector_field[i];
      if (from == to)
        continue;
      auto frame = frame_fromz(from, normalize(to - from));
      frame.z *= length(to - from);
      frames.push_back(frame);
    }
  } else {
    for (int i = 0; i < app.mesh.triangles.size(); i++) {
      auto from = tid_centroid(app.mesh.triangles, app.mesh.positions, i);
      auto to = from + scale * vector_field[i];
      if (from == to)
        continue;
      auto frame = frame_fromz(from, normalize(to - from));
      frame.z *= length(to - from);
      frames.push_back(frame);
    }
  }
  // auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  auto cylinder = make_uvarrow({4, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }
  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);

  added_path->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);

  return added_path;
}
Added_Path *add_generic_vector_field_shape(
    App &app, const vector<vec3f> &vector_field, const vector<vec3f> &instances,
    const vector<vec3f> &normals, const float &scale, const float &radius,
    const vec3f &color, const float &offset) {
  auto added_path = app.added_paths.emplace_back(new Added_Path{});
  added_path->color = color;
  added_path->radius = radius;
  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  auto frames = vector<frame3f>();
  frames.reserve(vector_field.size());
  for (auto i = 0; i < instances.size(); ++i) {
    auto from = instances[i] + offset * normals[i];
    auto to = from + scale * vector_field[i];
    if (from == to)
      continue;
    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }
  // auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  auto cylinder = make_uvarrow({4, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }
  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);

  added_path->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);

  return added_path;
}
Added_Path *add_glyph_shape(App &app, const vector<vec3f> &alpha,
                            const vector<vec3f> &instances,
                            const vector<vec3f> &normals, const float &scale,
                            const float &radius, const vec3f &color,
                            const float &offset) {
  auto added_path = app.added_paths.emplace_back(new Added_Path{});
  added_path->color = color;
  added_path->radius = radius;
  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  auto frames = vector<frame3f>();
  frames.reserve(alpha.size());
  for (auto i = 0; i < instances.size(); ++i) {
    auto from = instances[i] + offset * normals[i] - alpha[i] * scale;
    auto to = instances[i] + offset * normals[i] + alpha[i] * scale;
    if (from == to)
      continue;
    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }
  auto cylinder = make_rounded_uvcylinder({4, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }
  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);

  added_path->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);

  return added_path;
}
void update_vector_field_shape(shade_shape *shape, shape_data &mesh,
                               const vector<vec3f> &vector_field,
                               const float &scale, const float &radius,
                               const vec3f &color) {
  auto frames = vector<frame3f>();

  frames.reserve(vector_field.size());

  if (vector_field.size() == mesh.positions.size()) {
    for (int i = 0; i < mesh.positions.size(); i++) {
      auto from = mesh.positions[i];
      auto to = from + scale * vector_field[i];
      if (from == to)
        continue;
      auto frame = frame_fromz(from, normalize(to - from));
      frame.z *= length(to - from);
      frames.push_back(frame);
    }
  } else {
    for (int i = 0; i < mesh.triangles.size(); i++) {
      auto from = tid_centroid(mesh.triangles, mesh.positions, i);
      auto to = from + scale * vector_field[i];
      if (from == to)
        continue;
      auto frame = frame_fromz(from, normalize(to - from));
      frame.z *= length(to - from);
      frames.push_back(frame);
    }
  }
  // auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  auto cylinder = make_uvarrow({4, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}
void update_generic_vector_field_shape(shade_shape *shape,
                                       const vector<vec3f> &vector_field,
                                       const vector<vec3f> &instances,
                                       const vector<vec3f> &normals,
                                       const float &scale, const float &radius,
                                       const vec3f &color,
                                       const float &offset) {
  auto frames = vector<frame3f>();

  frames.reserve(vector_field.size());
  for (auto i = 0; i < instances.size(); ++i) {
    auto from = instances[i] + offset * normals[i];
    auto to = from + scale * vector_field[i];
    if (from == to)
      continue;
    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }

  // auto cylinder = make_uvcylinder({4, 1, 1}, {radius, 1});
  auto cylinder = make_uvarrow({4, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}
void update_glyph_shape(shade_shape *shape, const vector<vec3f> &alpha,
                        const vector<vec3f> &instances,
                        const vector<vec3f> &normals, const float &scale,
                        const float &radius, const vec3f &color,
                        const float &offset) {
  auto frames = vector<frame3f>();

  frames.reserve(alpha.size());
  for (auto i = 0; i < instances.size(); ++i) {
    auto from = instances[i] + offset * normals[i] - alpha[i] * scale;
    auto to = instances[i] + offset * normals[i] + alpha[i] * scale;
    if (from == to)
      continue;
    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }

  auto cylinder = make_rounded_uvcylinder({4, 1, 1}, {radius, 1});

  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius) {
  // if (thin) {
  //   set_positions(shape, positions);
  //   shape->shape->elements = ogl_element_type::line_strip;
  //   set_instances(shape, {});
  //   return;
  // }

  auto frames = vector<frame3f>();
  frames.reserve(positions.size() - 1);
  //  froms.reserve(positions.size() - 1);
  //  tos.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    auto from = positions[i];
    auto to = positions[i + 1];
    if (from == to)
      continue;

    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }

  auto cylinder = make_uvcylinder({8, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}

void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3f> &positions, float radius,
                       const float &treshold) {
  // if (thin) {
  //   set_positions(shape, positions);
  //   shape->shape->elements = ogl_element_type::line_strip;
  //   set_instances(shape, {});
  //   return;
  // }

  auto frames = vector<frame3f>();
  frames.reserve(positions.size() - 1);
  //  froms.reserve(positions.size() - 1);
  //  tos.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    auto from = positions[i];
    auto to = positions[i + 1];
    if (from == to || length(to - from) >= treshold)
      continue;

    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }

  auto cylinder = make_uvcylinder({8, 1, 1}, {radius, 1});
  for (auto &p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}
void update_path_shape(shade_shape *shape, const shape_data &mesh,
                       const vector<vec3i> &adjacencies,
                       const geodesic_path &path, float radius) {
  auto positions =
      path_positions(path, mesh.triangles, mesh.positions, adjacencies);
  update_path_shape(shape, mesh, positions, radius);
}

Added_Path *add_path_shape(App &app, const geodesic_path &path, float radius,
                           const vec3f &color) {
  auto added_path = app.added_paths.emplace_back(new Added_Path{});
  added_path->path = path;
  added_path->color = color;
  added_path->radius = radius;

  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_path_shape(shape, app.mesh, app.topology.adjacencies, path, radius);
  added_path->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);

  return added_path;
}
Added_Path *add_path_shape(App &app, const vector<vec3f> &positions,
                           float radius, const vec3f &color,
                           const float &threshold) {
  auto added_path = app.added_paths.emplace_back(new Added_Path{});
  added_path->color = color;
  added_path->radius = radius;
  added_path->positions = positions;
  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_path_shape(shape, app.mesh, positions, radius, threshold);
  added_path->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);

  return added_path;
}
void update_points_shape(shade_shape *shape, const shape_data &mesh,
                         const vector<mesh_point> &points, float radius) {
  auto positions = vector<vec3f>(points.size());
  for (int i = 0; i < positions.size(); i++) {
    positions[i] = eval_position(mesh.triangles, mesh.positions, points[i]);
  }
  update_points_shape(shape, positions, radius);
}
void update_points_shape(const vector<Added_Points *> &points,
                         const shape_data &data,
                         const vector<mesh_point> &new_points,
                         const int entry) {

  auto &shape = points[entry]->instance->shape;
  auto frames = vector<frame3f>(new_points.size(), identity3x4f);
  frames.reserve(new_points.size() - 1);
  for (int i = 0; i < new_points.size(); i++) {
    frames[i].o = eval_position(data.triangles, data.positions, new_points[i]);
  }

  auto sphere = make_sphere(32, points[entry]->radius);

  set_quads(shape, sphere.quads);
  set_positions(shape, sphere.positions);
  set_normals(shape, sphere.normals);
  set_texcoords(shape, sphere.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
  points[entry]->points = new_points;
}
Added_Points *add_points_shape(App &app, const vector<mesh_point> &points,
                               float radius, const vec3f &color) {
  auto added_points = app.added_points.emplace_back(new Added_Points{});
  added_points->points = points;
  added_points->color = color;
  added_points->radius = radius;

  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_points_shape(shape, app.mesh, points, radius);
  added_points->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);
  return added_points;
}
Added_Points *add_points_shape(App &app, const vector<vec3f> &points,
                               float radius, const vec3f &color) {
  auto added_points = app.added_points.emplace_back(new Added_Points{});
  // added_points->points = points; TODO(giacomo): fix this
  added_points->color = color;
  added_points->radius = radius;

  auto shape = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_points_shape(shape, points, radius);
  added_points->instance =
      add_instance(app.scene, identity3x4f, shape, material, false);
  return added_points;
}

void update_points_shape(shade_shape *shape, const vector<vec3f> &positions,
                         float radius) {
  auto frames = vector<frame3f>(positions.size(), identity3x4f);
  frames.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size(); i++) {
    frames[i].o = positions[i];
  }

  auto sphere = make_sphere(32, radius);

  set_quads(shape, sphere.quads);
  set_positions(shape, sphere.positions);
  set_normals(shape, sphere.normals);
  set_texcoords(shape, sphere.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}

void update_glpoints(App &app, const vector<vec3f> &positions,
                     const string &name) {
  auto &curr_shape = app.gpu_shapes[name];
  delete_shape(curr_shape);
  app.gpu_shapes[name] = gpu::make_points_shape(positions);
}
void update_glpoints(App &app, const vector<mesh_point> &points,
                     const string &name) {
  delete_shape(app.gpu_shapes[name]);
  if (points.empty())
    return;
  auto positions = vector<vec3f>(points.size());
  auto colors = vector<vec3f>(points.size(), {1, 1, 1});
  for (int i = 0; i < positions.size(); i++) {
    positions[i] =
        eval_position(app.mesh.triangles, app.mesh.positions, points[i]);
  }
  colors[0] = {0, 0.7, 0};
  if (points.size() >= 2)
    colors[1] = {0.7, 0, 0};
  app.gpu_shapes[name] = gpu::make_points_shape(positions);
  gpu::add_vertex_attribute(app.gpu_shapes[name], colors);
}

void update_glvector_field(App &app, const vector<vec3f> &vector_field,
                           float &scale, const string &name) {
  delete_shape(app.gpu_shapes[name]);
  auto &triangles = app.mesh.triangles;
  auto &positions = app.mesh.positions;
  if (vector_field.size() == triangles.size()) {
    app.gpu_shapes[name] =
        gpu::make_vector_field_shape(vector_field, triangles, positions, scale);
  } else {
    app.gpu_shapes[name] =
        gpu::make_vector_field_shape(vector_field, positions, scale);
  }
}
void save_curve(const App &app) {
  std::ofstream outfile;
  auto name = app.models_name + "_curve";
  outfile.open(name);
  if (outfile.fail())
    std::cout << "didn't mangae to open the file" << std::endl;

  for (auto i = 0; i < app.polyline.size(); ++i) {
    for (auto j = 0; j < app.polyline[i].size(); ++j)
      outfile << app.polyline[i][j].face << " " << app.polyline[i][j].uv.x
              << " " << app.polyline[i][j].uv.y << "\n";
  }

  outfile.close();
}
vector<vec3f> load_curve(App &app) {

  std::ifstream f;
  auto name = app.models_name;
  f.open(name);
  if (!f.is_open()) {
    std::cerr << "File loading the landamrks has not been found";
  }

  // f.seekg(0, f.end);
  // int len = f.tellg();
  // f.seekg(0, f.beg);
  auto path = vector<mesh_point>{};
  while (!f.eof()) {
    int tid;
    float uvx, uvy;
    f >> tid >> uvx >> uvy;
    path.push_back({tid, {uvx, uvy}});
  }
  f.close();
  // auto result = vector<vec3f>{};
  // for (auto i = 0; i < path.size() - 1; ++i) {
  //   auto curve = path_positions(
  //       compute_geodesic_path(app.mesh, app.topology,
  //                             app.dual_solvers[0], path[i], path[i + 1]),
  //       app.mesh.triangles, app.mesh.positions,
  //       app.topology.adjacencies);
  //   result.insert(result.end(), curve.begin(), curve.end());
  // }

  // return result;
  return polyline_pos(app.mesh.triangles, app.mesh.positions, path);
}

void save_control_points(const App &app) {

  std::ofstream outfile;
  outfile.open(app.models_name);
  if (outfile.fail())
    std::cout << "didn't mangae to open the file" << std::endl;
  for (auto i = 0; i < app.control_points.size(); ++i) {
    outfile << app.control_points[i].face << " " << app.control_points[i].uv.x
            << " " << app.control_points[i].uv.y << "\n";
  }

  outfile.close();
}

void load_control_points(App &app) {

  std::ifstream f;
  f.open(app.models_name);
  if (!f.is_open()) {
    std::cerr << "File loading the landamrks has not been found";

  } else {

    app.control_points.clear();

    while (!f.eof()) {
      int tid;
      float uvx, uvy;
      f >> tid >> uvx >> uvy;
      app.control_points.push_back({tid, {uvx, uvy}});
    }
    app.control_points.pop_back();
    f.close();
    for (auto i = 0; i < app.control_points.size(); ++i) {
      app.point_to_shape.insert(app.point_to_shape.end(),
                                {i, (int)app.added_points.size()});
      app.w.push_back(1);
      app.point_moved.push_back(true);
      app.control_points_shape.push_back(
          add_points_shape(app, {app.control_points[i]}, 0.0030, {0, 0, 0}));
    }

    std::cout << "Successfully loaded"
              << " " << app.control_points.size() << std::endl;
  }
}
