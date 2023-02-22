#include "mesh_io.h"

#include "ext/json.hpp"
#include <logging.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_sceneio.h>

using json = nlohmann::json;

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------

inline bool save_json(const string &filename, const json &js, string &error) {
  return save_text(filename, js.dump(2), error);
}

inline bool load_json(const string &filename, json &js, string &error) {
  // error helpers
  auto parse_error = [filename, &error]() -> json {
    error = filename + ": parse error in json";
    printf("error loading json %s\n:  %s\n", filename.c_str(), error.c_str());
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error))
    return false;
  try {
    js = json::parse(text);
    return true;
  } catch (std::exception &e) {
    return parse_error();
  }

  return true;
}

// support for json conversions

#define HEAVY 1
#if HEAVY
#include <diff_geo/diff_geo.h>

void init_mesh(shape_data &mesh, shape_geometry &geometry, shape_op &operators,
               geodesic_solver &solver, dual_geodesic_solver &dual_solver) {
  geometry.v2t =
      vertex_to_triangles(mesh.triangles, mesh.positions, geometry.adjacencies);

  dual_solver = make_dual_geodesic_solver(mesh.triangles, mesh.positions,
                                          geometry.adjacencies);
  geometry.angles = compute_angles_wo_opposite(
      mesh.triangles, mesh.positions, geometry.adjacencies, geometry.v2t,
      geometry.total_angles);
  solver = extended_solver(mesh, dual_solver, geometry, 6);

  operators.AGS_Grad = AGS_matrix(mesh, geometry);
}
#endif

#include "stl_reader.h"

bool load_mesh(const string &filename, shape_data &mesh,
               shape_geometry &geometry, shape_op &operators,
               geodesic_solver &solver, dual_geodesic_solver &dual_solver,
               string &error) {
  mesh = load_shape(filename);

  auto bbox = invalidb3f;
  for (auto &p : mesh.positions)
    bbox = merge(bbox, p);
  auto center = (bbox.max + bbox.min) / 2;
  auto scale = 1.0f / max(bbox.max - bbox.min);
  for (auto &p : mesh.positions)
    p = (p - center) * scale;

  geometry.adjacencies = face_adjacencies(mesh.triangles);

  mesh.normals = compute_normals(mesh.triangles, mesh.positions);

  init_mesh(mesh, geometry, operators, solver, dual_solver);

  return true;
}

#define NANOSVG_ALL_COLOR_KEYWORDS
#define NANOSVG_IMPLEMENTATION
#include "../nanosvg/src/nanosvg.h"

Svg load_svg(const string &filename) {
  struct NSVGimage *image;
  image = nsvgParseFromFile(filename.c_str(), "px", 96);
  printf("size: %f x %f\n", image->width, image->height);
  auto size = vec2f{image->width, image->height};

  // Use...
  auto svg = Svg();
  for (auto shape = image->shapes; shape != NULL; shape = shape->next) {
    auto &svg_shape = svg.emplace_back();
    unsigned int c;
    if (shape->fill.type == NSVG_PAINT_COLOR) {
      c = shape->fill.color;
    } else if (shape->fill.type >= NSVG_PAINT_LINEAR_GRADIENT) {
      c = shape->fill.gradient->stops[0].color;
    } else {
      c = 0;
    }
    float r = ((c >> 16) & 0xFF) / 255.0; // Extract the RR byte
    float g = ((c >> 8) & 0xFF) / 255.0;  // Extract the GG byte
    float b = ((c)&0xFF) / 255.0;
    svg_shape.color = yocto::pow(vec3f{b, g, r}, 2.2f);

    for (auto path = shape->paths; path != NULL; path = path->next) {
      auto &svg_path = svg_shape.paths.emplace_back();
      for (int i = 0; i < path->npts - 1; i += 3) {
        float *p = &path->pts[i * 2];
        auto &curve = svg_path.emplace_back();
        curve[0] = vec2f{p[0], size.y - p[1]} / size.y;
        curve[1] = vec2f{p[2], size.y - p[3]} / size.y;
        curve[2] = vec2f{p[4], size.y - p[5]} / size.y;
        curve[3] = vec2f{p[6], size.y - p[7]} / size.y;
        // printf("(%f %f) (%f %f) (%f %f) (%f %f)\n", curve[0].x, curve[0].y,
        //     curve[1].x, curve[1].y, curve[2].x, curve[2].y, curve[3].x,
        //     curve[3].y);
      }
    }
  }
  // Delete
  nsvgDelete(image);
  return svg;
}
