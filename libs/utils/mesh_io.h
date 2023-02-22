#pragma once
#include <diff_geo/diff_geo.h>

bool load_mesh(const string &filename, shape_data &mesh,
               shape_geometry &topology, shape_op &operators,
               geodesic_solver &solver, dual_geodesic_solver &dual_solver,
               string &error);

using Svg_Path = vector<array<vec2f, 4>>;
struct Svg_Shape {
  vec3f color = {};
  vector<Svg_Path> paths = {};
};
using Svg = vector<Svg_Shape>;

Svg load_svg(const string &filename);

#include "ext/json.hpp"

namespace yocto {

using json = nlohmann::json;
using std::array;

inline pair<Eigen::MatrixXd, Eigen::MatrixXi>
libigl_wrapper(const vector<vec3f> &positions, const vector<vec3i> &triangles) {
  Eigen::MatrixXd V(positions.size(), 3);
  Eigen::MatrixXi F(triangles.size(), 3);

  for (int i = 0; i < positions.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      V(i, j) = positions[i][j];
    }
  }
  for (int i = 0; i < triangles.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      F(i, j) = triangles[i][j];
    }
  }

  return {V, F};
}

// support for json conversions
inline void to_json(json &j, const vec2f &value) {
  nlohmann::to_json(j, (const array<float, 2> &)value);
}

inline void from_json(const json &j, vec2f &value) {
  nlohmann::from_json(j, (array<float, 2> &)value);
}

// support for json conversions
inline void to_json(json &js, const mesh_point &value) {
  js["face"] = value.face;
  js["uv"] = value.uv;
}

inline void from_json(const json &js, mesh_point &value) {
  js.at("face").get_to(value.face);
  js.at("uv").get_to(value.uv);
}

} // namespace yocto
