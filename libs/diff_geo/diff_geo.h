#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

} // namespace yocto
#include <VTP/geodesic_algorithm_exact.h>
#include <VTP/geodesic_mesh.h>
#include <iostream>
#include <stdio.h>
#include <unordered_set>
#include <utils/logging.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
using namespace yocto;
enum transport_mode { V2V, V2T, T2V, T2T };
enum simplex { vert, edge };
enum field { exact, graph };

struct shape_geometry {
  vector<vec3i> adjacencies = {};
  vector<vector<int>> v2t = {};
  vector<vector<float>> angles = {};
  vector<float> total_angles = {};
  float avg_edge_length = 0.f;
};

struct shape_op {
  Eigen::SparseMatrix<double> PCE_Grad;
  Eigen::SparseMatrix<double> AGS_Grad;
  Eigen::SparseMatrix<double> Grad;
  Eigen::SparseMatrix<double> Lap;
  Eigen::SparseMatrix<double> Hess;
  vector<Eigen::MatrixXd> C;
  vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> CMat;
  vector<Eigen::MatrixXd> Rhs;
  vector<Eigen::MatrixXd> quadrics;
};
int forced_vert_from_point(const vector<vec3i> &triangles, const mesh_point &p);
void clean_bary(mesh_point &sample);
vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int pid);

vector<int> one_ring(const vector<vec3i> &triangles,
                     const vector<vector<int>> &v2t, const int vid);

vector<vec3f> polyline_pos(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<mesh_point> &poly);

std::pair<float, float> intersect(const vec2f &direction, const vec2f &left,
                                  const vec2f &right);

vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p);
float tid_area(const vector<vec3i> &triangles, const vector<vec3f> &positions,
               const int tid);

Eigen::Matrix3d tranformation_matrix(const vec2f &x1, const vec2f &y1,
                                     const vec2f &O1);

Eigen::Matrix4d tranformation_matrix(const vec3f &x1, const vec3f &y1,
                                     const vec3f &z1, const vec3f &O1);

Eigen::Matrix3d switch_reference_frame(const shape_data &data,
                                       const shape_geometry &geometry,
                                       const int tid0, const int tid1);

vec3f switch_reference_frame(const Eigen::Matrix4d &T, const vec3f &p);

vec3f switch_reference_frame_vector(const Eigen::Matrix4d &T, const vec3f &p);

vec2f switch_reference_frame(const Eigen::Matrix3d &T, const vec2f &p);

vec2f switch_reference_frame_vector(const Eigen::Matrix3d &T, const vec2f &p);

std::tuple<int, float> straightest_path_in_tri(const shape_data &data,
                                               const mesh_point &p,
                                               const vec2f &dir,
                                               const int k_start = -1);

vector<mesh_point> straightest_path_from_vert_in_triangle(
    const shape_data &data, const shape_geometry &geometry, const int tid,
    const int vid, const float &len, const vec3f &dir);

vector<vec3f> get_basis_at_p(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid);
std::tuple<vec2f, mesh_point, int>
polthier_condition_at_vert(const shape_data &data,
                           const shape_geometry &geometry, const int vid,
                           const int tid, const vec3f &dir);
std::tuple<vector<vec3f>, mesh_point> polthier_straightest_geodesic(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3i> &adjacencies, const vector<vector<int>> &v2t,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &p, const vec2f &dir, const float &len,
    const vector<bool> &boundaries = {});
inline std::tuple<vector<vec3f>, mesh_point> polthier_straightest_geodesic(
    const shape_data &data, const shape_geometry &geometry, const mesh_point &p,
    const vec2f &dir, const float &len, const vector<bool> &boundaries = {}) {
  return polthier_straightest_geodesic(
      data.triangles, data.positions, geometry.adjacencies, geometry.v2t,
      geometry.angles, geometry.total_angles, p, dir, len, boundaries);
}

mesh_point polthier_straightest_geodesic_for_Newton_method(
    const shape_data &data, const shape_geometry &geometry, const mesh_point &p,
    const vector<vector<float>> &f, const vector<float> &w, const vec2f &dir,
    const float &len, const vector<bool> &boundaries = {});

vec3f project_vec(const vec3f &v, const vec3f &n);

mesh_point point_from_vert(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int j = 0);
mesh_point point_from_vert(const vector<vec3i> &triangles, const int vid,
                           const int tid);
int forced_vert_from_point(const vector<vec3i> &triangles, const mesh_point &p);
vector<float> subdivide_angles(const int number_of_subdivision,
                               const vec2f &range = vec2f{0, 2 * pif});

vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle);
vec2f rot_vect(const vec2f &p, const float theta);
mesh_point make_mesh_point(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int tid = -1);

float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vector<int>> &v2t,
                             const vector<vec3f> &normals, const int vid,
                             const vec3f &v);

vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions,
                  const vector<vector<int>> &v2t, const vector<vec3f> &normals,
                  int vid);
vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int tid);

std::pair<vector<double>, vector<int>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions);

std::pair<vector<double>, vector<int>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions,
                       const mesh_point &source);

Eigen::VectorXd wrapper(const vector<float> &f);

float nbr_avg_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, const int vid);
float nbr_min_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, int vid);
float avg_edge_length(const shape_data &mesh, const shape_geometry &topology);

vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source);
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const mesh_point &source);

vec3f compute_PCE(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, const vector<float> &f,
                  const mesh_point &p);
Eigen::SparseMatrix<double> PCE_matrix(const vector<vec3i> &triangles,
                                       const vector<vec3f> &positions);

Eigen::SparseMatrix<double> AGS_matrix(const shape_data &data,
                                       const shape_geometry &geometry);

Eigen::VectorXd compute_hessian(const shape_op &operators,
                                const vector<float> &field);
Eigen::Matrix2d compute_hessian(const Eigen::VectorXd &Hessian, const int vid);

void parallel_transp(const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies,
                     const vector<vector<int>> &v2t, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode);

geodesic_path compute_geodesic_path(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const dual_geodesic_solver &solver,
                                    const mesh_point &start,
                                    const mesh_point &end);

geodesic_solver extended_solver(const shape_data &data,
                                const dual_geodesic_solver &dual_solver,
                                shape_geometry &geometry, const int k);

inline vec3f tri_bary_coords(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid,
                             const vec3f &p) {
  auto px = positions[triangles[tid].x];
  auto py = positions[triangles[tid].y];
  auto pz = positions[triangles[tid].z];
  return tri_bary_coords(px, py, pz, p);
}

inline vec3f tid_normal(const vector<vec3i> &triangles,
                        const vector<vec3f> &positions, const int tid) {
  auto p0 = positions[triangles[tid].x];
  auto p1 = positions[triangles[tid].y];
  auto p2 = positions[triangles[tid].z];

  return normalize(cross(p1 - p0, p2 - p0));
}

inline float tid_area(const shape_data &mesh, const int tid) {
  return tid_area(mesh.triangles, mesh.positions, tid);
}

inline int node_is_adjacent(const geodesic_solver &solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}
template <typename T> inline int find(const T &vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x)
      return i;
  return -1;
}
inline vec3f get_bary(const vec2f &uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}
