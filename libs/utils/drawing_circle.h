#pragma once

#include "struct.h"
#include <realtime/gpu.h>
#include <vector>
#include <yocto/yocto_mesh.h>
using gpu::Shape, gpu::Shader;
using std::vector;
using yocto::mesh_point, yocto::vec3i, yocto::vec3f, yocto::mat4f, yocto::vec2i;

struct closed_curve {
  vector<int> strip = {};
  vector<float> lerps = {};
};
struct circle_tids {
  float lerp = -1;
  int offset = -1;
};
typedef vector<closed_curve> Isoline;

struct Circle {
  mesh_point center;
  float radius;
  Isoline isoline;
  vector<int> tids = {};
  vector<vec3f> pos = {};
  vector<float> distances = {};
  float theta = 0.f;
  int primitive = -1;
  vec2i crop_range = vec2i{-1, -1};
  vector<mesh_point> vertices = {};
  vector<vector<mesh_point>> spider_vertices = {};
  vector<vector<vec3f>> spider_arcs = {};
  int garland_vertices = -1;
  int levels = -1;
};
struct Ellipse {
  mesh_point midpoint;
  float scale_factor;
  Isoline isoline;
  vector<float> distances_from_midpoint = {};
};

Circle create_circle(const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies, const mesh_point &center,
                     const float &radius, const vector<float> &distances);

vector<vec3f> circle_positions(const vector<vec3i> &triangles,
                               const vector<vec3f> &positions,
                               const vector<vec3i> &adjacencies,
                               const Circle &c0);
template <typename T> inline int find_in_vec(const T &vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x)
      return i;
  return -1;
}