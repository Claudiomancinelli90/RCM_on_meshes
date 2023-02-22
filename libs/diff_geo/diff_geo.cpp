//
//  karcher.cpp
//  glfw
//
//  Created by Claudio Mancinelli on 21/07/2020.
//

#include "diff_geo.h"
using namespace logging;
using namespace yocto;

#include <deque>

// utility (geometry)
// project a vector v onto a plane having n as normal
vec3f project_vec(const vec3f &v, const vec3f &n) {
  auto proj = n * dot(v, n);

  return v - proj;
}
vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle) {
  auto M = rotation_frame(axis, angle);
  auto v = p;
  return transform_vector(M, v);
}
vec2f rot_vect(const vec2f &p, const float theta) {
  auto M = mat2f{{yocto::cos(theta), -yocto::sin(theta)},
                 {yocto::sin(theta), yocto::cos(theta)}};
  auto v = M * p;
  return v;
}

vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int tid) {
  vec3f p0 = positions[triangles[tid].x];
  vec3f p1 = positions[triangles[tid].y];
  vec3f p2 = positions[triangles[tid].z];

  return (p0 + p1 + p2) / 3.0;
}
float tid_area(const vector<vec3i> &triangles, const vector<vec3f> &positions,
               const int tid) {
  return triangle_area(positions[triangles[tid].x], positions[triangles[tid].y],
                       positions[triangles[tid].z]);
}

vector<vec3f> polyline_pos(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<mesh_point> &poly) {
  auto result = vector<vec3f>(poly.size());
  for (auto i = 0; i < poly.size(); ++i) {
    result[i] = eval_position(triangles, positions, poly[i]);
  }
  return result;
}
int forced_vert_from_point(const vector<vec3i> &triangles,
                           const mesh_point &p) {
  auto bary = vector<pair<float, int>>{
      {1 - p.uv.x - p.uv.y, 0}, {p.uv.x, 1}, {p.uv.y, 2}};

  sort(bary.begin(), bary.end());
  return triangles[p.face][bary.back().second];
}

void clean_bary(mesh_point &sample) {

  auto bary = sample.uv;
  auto coords = vector<pair<float, int>>{
      {1 - bary.x - bary.y, 0}, {bary.x, 1}, {bary.y, 2}};
  sort(coords.begin(), coords.end());
  vec3f bary3d = get_bary(bary);
  if (coords[0].first < 0 && coords[1].first > 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[2].second] =
        1 - bary3d[coords[1].second] - bary3d[coords[0].second];
  } else if (coords[0].first < 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[1].second] = 0;
    bary3d[coords[2].second] = 1;
  }

  sample = {sample.face, {bary3d.y, bary3d.z}};
}

float avg_edge_length(const shape_data &data, const shape_geometry &topology) {
  auto avg = 0.f;
  auto count = 0.f;
  for (auto i = 0; i < data.positions.size(); ++i) {
    auto star = topology.v2t[i];
    for (auto tid : star) {
      auto k = find(data.triangles[tid], i);
      auto vid0 = data.triangles[tid][(k + 1) % 3];
      if (vid0 < i)
        continue;

      avg += length(data.positions[vid0] - data.positions[i]);
      ++count;
    }
  }
  return avg / count;
}
vec3f polar_basis(const geodesic_solver &solver, const vector<vec3f> &positions,
                  const vector<vec3f> &normals, int vid) {
  int vid0 = solver.graph[vid][0].node;
  vec3f v = positions[vid0] - positions[vid];
  vec3f e = project_vec(v, normals[vid]);
  return normalize(e);
}
vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions,
                  const vector<vector<int>> &v2t, const vector<vec3f> &normals,
                  int vid) {
  auto tid = v2t[vid][0];
  auto k = find(triangles[tid], vid);
  vec3f v = positions[triangles[tid][(k + 1) % 3]] - positions[vid];
  vec3f e = normalize(project_vec(v, normals[vid]));
  return e;
}

vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int tid) {
  auto c = tid_centroid(triangles, positions, tid);
  vec3f v = positions[triangles[tid].x];
  return normalize(v - c);
}
// Transformation matrix s.t [x,y,1]=T[x1,y1,1]
Eigen::Matrix3d tranformation_matrix(const vec2f &x1, const vec2f &y1,
                                     const vec2f &O1) {
  Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
  T << x1.x, y1.x, O1.x, x1.y, y1.y, O1.y, 0, 0, 1;

  return T;
}

Eigen::Matrix4d tranformation_matrix(const vec3f &x1, const vec3f &y1,
                                     const vec3f &z1, const vec3f &O1) {
  Eigen::Matrix4d T;

  T << x1.x, x1.y, x1.z, -dot(x1, O1), y1.x, y1.y, y1.z, -dot(y1, O1), z1.x,
      z1.y, z1.z, -dot(z1, O1), 0, 0, 0, 1;

  return T;
}
Eigen::Matrix3d switch_reference_frame(const vector<vec3i> &triangles,
                                       const vector<vec3f> &positions,
                                       const int tid0, const int tid1) {
  auto flat_t0 = init_flat_triangle(positions, triangles[tid0], 0);
  auto flat_t1 = unfold_face(triangles, positions, flat_t0, tid0, tid1);
  auto e0 = normalize(flat_t1[1] - flat_t1[0]);
  auto e1 = vec2f{-e0.y, e0.x};
  return tranformation_matrix(e0, e1, flat_t1[0]);
}
Eigen::Matrix3d switch_reference_frame(const shape_data &data,
                                       const shape_geometry &geometry,
                                       const int tid0, const int tid1) {
  auto flat_t0 = init_flat_triangle(data.positions, data.triangles[tid0], 0);
  auto flat_t1 =
      unfold_face(data.triangles, data.positions, flat_t0, tid0, tid1);
  auto e0 = normalize(flat_t1[1] - flat_t1[0]);
  auto e1 = vec2f{-e0.y, e0.x};
  return tranformation_matrix(e0, e1, flat_t1[0]);
}
vec3f switch_reference_frame(const Eigen::Matrix4d &T, const vec3f &p) {
  Eigen::Vector4d V;
  V << p.x, p.y, p.z, 1;
  Eigen::Vector4d t = T * V;
  return vec3f{(float)t(0), (float)t(1), (float)t(2)};
}

vec3f switch_reference_frame_vector(const Eigen::Matrix4d &T, const vec3f &p) {
  Eigen::Vector4d V;
  V << p.x, p.y, p.z, 0;
  Eigen::Vector4d t = T * V;
  return vec3f{(float)t(0), (float)t(1), (float)t(2)};
}

vec2f switch_reference_frame(const Eigen::Matrix3d &T, const vec2f &p) {
  Eigen::Vector3d V;
  V << p.x, p.y, 1;
  Eigen::Vector3d t = T * V;
  return vec2f{(float)t(0), (float)t(1)};
}
vec2f switch_reference_frame_vector(const Eigen::Matrix3d &T, const vec2f &p) {
  Eigen::Vector3d V;
  V << p.x, p.y, 0;
  Eigen::Vector3d t = T * V;
  return vec2f{(float)t(0), (float)t(1)};
}
// Compute polar coordinates
float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vector<int>> &v2t,
                             const vector<vec3f> &normals, const int vid,
                             const vec3f &v) {
  auto e = polar_basis(triangles, positions, v2t, normals, vid);

  auto teta = angle(v, e);

  if (dot(cross(e, v), normals[vid]) < 0)
    teta *= -1;

  return teta;
}
float angle_in_tangent_space(const geodesic_solver &solver,
                             const vector<vec3f> &positions, const vec3f &v,
                             const int vid, const vec3f &n) {
  float teta;
  int vid0 = solver.graph[vid][0].node;
  vec3f e0 = positions[vid0] - positions[vid];
  vec3f e = normalize(project_vec(e0, n));

  teta = angle(v, e);
  if (dot(cross(e, v), n) < 0)
    teta *= -1;

  return teta;
}
float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid,
                             const vec3f &v) {
  auto teta = 0.f;

  auto e = polar_basis(triangles, positions, tid);
  auto n = tid_normal(triangles, positions, tid);
  teta = angle(v, e);
  if (dot(cross(e, v), n) < 0)
    teta = 2 * M_PI - teta;

  return teta;
}
mesh_point point_from_vert(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int j) {
  auto tid = v2t[vid][j];
  auto k = find(triangles[tid], vid);
  auto bary = zero3f;
  bary[k] = 1;
  return {tid, {bary.y, bary.z}};
}
mesh_point point_from_vert(const vector<vec3i> &triangles, const int vid,
                           const int tid) {
  auto k = find(triangles[tid], vid);
  auto bary = zero3f;
  bary[k] = 1;
  return {tid, {bary.y, bary.z}};
}
int vert_from_point(const vector<vec3i> &triangles, const mesh_point &p) {
  auto b = vector<pair<float, int>>(3);
  auto bary = get_bary(p.uv);
  for (auto i = 0; i < 3; ++i) {
    b[i] = std::make_pair(bary[i], i);
  }
  sort(b.begin(), b.end());

  return triangles[p.face][b.back().second];
}
vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p) {
  vec3f wgts = vec3f{0.0, 0.0, 0.0};
  vec3f u = v1 - v0, v = v2 - v0, w = p - v0;
  float d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
        d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0)
    return zero3f;

  wgts[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(wgts[2]));
  wgts[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(wgts[1]));
  wgts[0] = 1.0 - wgts[1] - wgts[2];
  assert(!isnan(wgts[0]));

  return wgts;
}
vector<int> one_ring(const vector<vec3i> &triangles,
                     const vector<vector<int>> &v2t, const int vid) {
  auto &star = v2t[vid];
  auto nbr = vector<int>(star.size());
  for (auto i = 0; i < star.size(); ++i) {
    auto tid = star[i];
    auto tr = triangles[tid];
    auto k = find(tr, vid);
    nbr[i] = tr[(k + 1) % 3];
  }
  return nbr;
}
vector<int> k_ring(const vector<vec3i> &triangles,
                   const vector<vector<int>> &v2t, const int vid, const int k) {
  vector<int> ring;
  vector<int> active_set = {vid};
  for (int i = 0; i < k; ++i) {
    vector<int> next_active_set;
    for (int j = 0; j < active_set.size(); ++j) {
      auto nbr = one_ring(triangles, v2t, active_set[j]);
      for (int h = 0; h < nbr.size(); ++h) {
        int curr = nbr[h];
        if (find(ring, curr) == -1 && curr != vid) {
          next_active_set.push_back(curr);
          ring.push_back(curr);
        }
      }
    }
    active_set = next_active_set;
  }

  return ring;
}
int opposite_vert(const vector<vec3i> &triangles,
                  const vector<vec3i> &adjacencies, const int tid,
                  const int vid) {
  auto nei = opposite_face(triangles, adjacencies, tid, vid);
  auto h = find(adjacencies[nei], tid);
  return triangles[nei][(h + 2) % 3];
}

void parallel_transp(const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies,
                     const vector<vector<int>> &v2t, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode) {
  switch (mode) {
  case V2V: {
    auto mag = length(v);
    float teta =
        angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
    auto nbr_from = one_ring(triangles, v2t, from);
    auto nbr_to = one_ring(triangles, v2t, to);
    float phi_ij = -1;
    float phi_ji = -1;

    for (int i = 0; i < nbr_from.size(); ++i) {
      if (nbr_from[i] == to) {
        phi_ij = angles[from][i];

        break;
      }
    }

    for (int j = 0; j < nbr_to.size(); ++j) {
      if (nbr_to[j] == from) {
        phi_ji = angles[to][j];
        break;
      }
    }
    assert(phi_ij != -1);
    assert(phi_ji != -1);

    vec3f e0 = polar_basis(triangles, positions, v2t, normals, to);
    float rotation = teta + phi_ji + M_PI - phi_ij;

    v = rot_vect(e0, normals[to], rotation);
    v *= mag;

  } break;

  case V2T: {
    float teta =
        angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
    vec3i tri = triangles[to];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f normal = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);

    vec3f coords = positions[from] - centroid;
    float phi_ji = angle(e, coords);
    if (dot(cross(e, coords), normal) < 0)
      phi_ji = 2 * M_PI - phi_ji;
    int offset = find(tri, from);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];
    int vid2 = tri[(offset + 2) % 3];
    float factor = 2 * M_PI / total_angles[from];
    auto nbr_from = one_ring(triangles, v2t, from);
    float phi_ij = -1;
    coords *= -1;
    if (nbr_from[0] == vid2) {
      vec3f edge = positions[vid2] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      curr_angle = 2 * M_PI - curr_angle;
      phi_ij = curr_angle;
    } else {
      for (int i = 0; i < nbr_from.size(); ++i) {
        if (nbr_from[i] == vid1) {
          phi_ij = angles[from][i];
          break;
        }
      }

      vec3f edge = positions[vid1] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      phi_ij += curr_angle;
    }

    float rot = teta + phi_ji + M_PI - phi_ij;

    e *= length(v);
    v = rot_vect(e, normal, rot);

  }

  break;

  case T2V: {
    vec3i tri = triangles[from];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f n = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);
    float teta = angle(e, v);

    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;
    int offset = find(tri, to);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];

    vec3f vert = positions[tri[offset]];
    vec3f v1 = positions[vid1] - vert;

    vec3f coords = vert - centroid;
    float phi_ij = angle(e, coords);
    if (dot(cross(e, coords), n) < 0)
      phi_ij = 2 * M_PI - phi_ij;

    coords *= -1;
    float phi_ji = angle(v1, coords);
    float factor = 2 * M_PI / total_angles[to];
    phi_ji *= factor;
    auto nbr = one_ring(triangles, v2t, to);
    for (int i = 0; i < nbr.size(); ++i) {
      if (nbr[i] == vid1) {
        float phi = angles[to][i];
        phi_ji += phi;
        break;
      }
    }

    float rot = teta + phi_ji + M_PI - phi_ij;
    vec3f e0 = polar_basis(triangles, positions, v2t, normals, to);
    e0 *= length(v);
    v = rot_vect(e0, normals[to], rot);

  }

  break;

  case T2T: {
    auto flat_from = init_flat_triangle(positions, triangles[from]);
    auto flat_to = unfold_face(triangles, positions, flat_from, from, to);
    auto bary = vec2f{0.333, 0.333};
    auto c0 =
        interpolate_triangle(flat_from[0], flat_from[1], flat_from[2], bary);
    auto c1 = interpolate_triangle(flat_to[0], flat_to[1], flat_to[2], bary);
    auto e0 = flat_from[0] - c0;
    auto e1 = flat_to[0] - c1;

    auto w = c1 - c0;
    auto phi_ij = angle(e0, w);
    if (cross(e0, w) < 0)
      phi_ij = 2 * M_PI - phi_ij;
    w *= -1;
    auto phi_ji = angle(e1, w);

    if (cross(e1, w) < 0)
      phi_ji = 2 * M_PI - phi_ji;

    auto n = tid_normal(triangles, positions, from);
    auto e = polar_basis(triangles, positions, from);
    float teta = angle(e, v);
    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;

    auto e_to = polar_basis(triangles, positions, to);
    auto n_to = tid_normal(triangles, positions, to);
    float rot = teta + phi_ji + M_PI - phi_ij;
    e_to *= length(v);
    v = rot_vect(e_to, n_to, rot);

  }

  break;
  }
}
// VTP
std::pair<vector<double>, vector<int>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> points(3 * V);
  vector<int> faces(3 * F);
  vector<double> f(V);

  for (int i = 0; i < V; ++i) {
    points[3 * i] = positions[i].x;
    points[3 * i + 1] = positions[i].y;
    points[3 * i + 2] = positions[i].z;
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i] = triangles[i].x;
    faces[3 * i + 1] = triangles[i].y;
    faces[3 * i + 2] = triangles[i].z;
  }

  return std::make_pair(points, faces);
}
std::pair<vector<double>, vector<int>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions,
                       const mesh_point &source) {
  int V = (int)positions.size() + 1;
  int F = (int)triangles.size();
  auto tid = source.face;
  vector<double> points(3 * V);
  vector<int> faces(3 * (F + 2));
  vector<float> f(V);
  auto pos = eval_position(triangles, positions, source);

  for (int i = 0; i < V; ++i)
    // if (v2t[i].size() == 0) continue;
    if (i != V - 1) {
      points[3 * i] = positions[i].x;
      points[3 * i + 1] = positions[i].y;
      points[3 * i + 2] = positions[i].z;
    } else {
      points[3 * i] = pos.x;
      points[3 * i + 1] = pos.y;
      points[3 * i + 2] = pos.z;
    }

  for (int i = 0; i < F; ++i) {
    if (i != tid) {
      faces[3 * i] = triangles[i].x;
      faces[3 * i + 1] = triangles[i].y;
      faces[3 * i + 2] = triangles[i].z;
    } else {
      faces[3 * i] = triangles[i].x;
      faces[3 * i + 1] = triangles[i].y;
      faces[3 * i + 2] = V - 1;

      faces[3 * F] = triangles[i].y;
      faces[3 * F + 1] = triangles[i].z;
      faces[3 * F + 2] = V - 1;

      faces[3 * (F + 1)] = triangles[i].z;
      faces[3 * (F + 1) + 1] = triangles[i].x;
      faces[3 * (F + 1) + 2] = V - 1;
    }
  }
  return {points, faces};
}
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source) {
  int V = (int)positions.size();
  vector<float> f(V);
  auto [points, faces] = exact_geodesic_wrapper(triangles, positions);
  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v = verts[j];
    float value = (float)v.geodesic_distance();
    f[j] = value;
  }

  return f;
}

// utility (gradient matrix)
Eigen::VectorXd wrapper(const vector<float> &f) {
  Eigen::VectorXd F(f.size());
  for (int i = 0; i < f.size(); ++i) {
    F(i) = f[i];
  }
  return F;
}
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const mesh_point &source) {
  time_function();
  auto tid = source.face;
  auto [is_vert, offset] = point_is_vert(source);
  if (is_vert)
    return exact_geodesic_distance(triangles, positions,
                                   triangles[tid][offset]);
  else {
    int V = (int)positions.size() + 1;

    vector<float> f(V);
    auto [points, faces] = exact_geodesic_wrapper(triangles, positions, source);

    geodesic_VTP::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
    algorithm.propagate(V - 1);
    vector<geodesic_VTP::Vertex> verts = mesh.vertices();
    for (int j = 0; j < V; ++j) {
      geodesic_VTP::Vertex v = verts[j];
      float value = (float)v.geodesic_distance();
      f[j] = value;
    }
    f.pop_back();
    return f;
  }
}
Eigen::MatrixXd rhs(int s) {
  Eigen::MatrixXd E(s, s + 1);
  Eigen::MatrixXd X = Eigen::MatrixXd::Constant(s, 1, -1);
  E.topLeftCorner(s, 1) = X;
  Eigen::MatrixXd I(s, s);
  I.setIdentity();
  E.topRightCorner(s, s) = I;
  return E;
}
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
double angle_at_vertex(const shape_data &data, const int tid, const int vid) {
  auto k = find(data.triangles[tid], vid);
  if (k == -1) {
    std::cout << "Error: the vertex must belongs to the triangles. Check the "
                 "construction of the v2t adjacency"
              << std::endl;
    return 0.f;
  }
  return angle(data.positions[data.triangles[tid][(k + 1) % 3]] -
                   data.positions[data.triangles[tid][k]],
               data.positions[data.triangles[tid][(k + 2) % 3]] -
                   data.positions[data.triangles[tid][k]]);
}
inline double cot(const float &theta) {
  return yocto::cos(theta) / yocto::sin(theta);
}
Eigen::VectorXd compute_hessian(const shape_op &operators,
                                const vector<float> &field) {
  auto f = wrapper(field);
  Eigen::VectorXd Hessian = operators.Hess * f;
  return Hessian;
}
Eigen::Matrix2d compute_hessian(const Eigen::VectorXd &Hessian, const int vid) {
  Eigen::Matrix2d H;
  auto n = Hessian.rows() / 4;
  H(0, 0) = Hessian(vid);
  H(0, 1) = Hessian(n + vid);
  H(1, 0) = Hessian(2 * n + vid);
  H(1, 1) = Hessian(3 * n + vid);

  return H;
}
vec3f compute_PCE(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, const vector<float> &f,
                  const mesh_point &p) {
  auto tid = p.face;
  vec3f n = cross(positions[triangles[tid].y] - positions[triangles[tid].x],
                  positions[triangles[tid].z] - positions[triangles[tid].x]);
  double area = length(n);
  n = normalize(n);
  auto vid0 = triangles[tid][0];
  auto vid1 = triangles[tid][1];
  auto vid2 = triangles[tid][2];
  auto v0 = positions[vid0];
  auto v1 = positions[vid1];
  auto v2 = positions[vid2];

  return 1 / area *
         ((f[vid1] - f[vid0]) * cross(n, v0 - v2) +
          (f[vid2] - f[vid0]) * cross(n, v1 - v0));
}
vec3f PCE_contribute(const vector<vec3i> &triangles,
                     const vector<vec3f> &positions, const int tid,
                     const vec3f &n, const float &double_area, const int off) {
  int prev = triangles[tid][off];
  int curr = triangles[tid][(off + 1) % 3];
  int next = triangles[tid][(off + 2) % 3];
  vec3f u = positions[next] - positions[curr];
  vec3f v = positions[curr] - positions[prev];
  vec3f u_90 = normalize(cross(u, n));
  vec3f v_90 = normalize(cross(v, n));

  vec3f contribute = u_90 * length(u) + v_90 * length(v);
  contribute /= double_area;

  return contribute;
}
Eigen::SparseMatrix<double> PCE_matrix(const vector<vec3i> &triangles,
                                       const vector<vec3f> &positions) {
  time_function();
  Eigen::SparseMatrix<double> G(triangles.size() * 3, positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T> entries;

  for (int i = 0; i < triangles.size(); ++i) {
    vec3f n = cross(positions[triangles[i].y] - positions[triangles[i].x],
                    positions[triangles[i].z] - positions[triangles[i].x]);
    double area = length(n); // division by 2 is missing because
                             // there is a *2 in the computation of
                             // the gradient
    n = normalize(n);

    for (int off = 0; off < 3; ++off) {
      int prev = triangles[i][off];
      int curr = triangles[i][(off + 1) % 3];
      int next = triangles[i][(off + 2) % 3];
      vec3f u = positions[next] - positions[curr];
      vec3f v = positions[curr] - positions[prev];
      vec3f u_90 = normalize(cross(u, n));
      vec3f v_90 = normalize(cross(v, n));

      vec3f contribute = u_90 * length(u) + v_90 * length(v);
      contribute /= area;

      int row = 3 * i;
      entries.push_back(T(row, curr, contribute.x));
      ++row;
      entries.push_back(T(row, curr, contribute.y));
      ++row;
      entries.push_back(T(row, curr, contribute.z));
    }
  }
  G.setFromTriplets(entries.begin(), entries.end());
  return G;
}

std::pair<float, float> intersect(const vec2f &direction, const vec2f &left,
                                  const vec2f &right) {
  auto v1 = -left;
  auto v2 = right - left;
  auto v3 = vec2f{-direction.y, direction.x};
  auto t0 = cross(v2, v1) / dot(v2, v3);
  auto t1 = -dot(left, v3) / dot(v2, v3);
  return std::make_pair(t0, t1);
};

vector<vec3f> get_basis_at_p(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid) {
  auto v0 = positions[triangles[tid][0]];
  auto v1 = positions[triangles[tid][1]];
  auto v2 = positions[triangles[tid][2]];
  auto e = normalize(v1 - v0);
  auto n = triangle_normal(v0, v1, v2);
  auto e_perp = cross(n, e);
  return {e, e_perp, n};
}
std::tuple<int, float> straightest_path_in_tri(const vector<vec3i> &triangles,
                                               const vector<vec3f> &positions,
                                               const mesh_point &p,
                                               const vec2f &dir,
                                               const int k_start) {
  auto tid = p.face;
  auto ep = get_basis_at_p(triangles, positions, tid);
  auto vid0 = -1;
  auto vid1 = -1;
  auto vid2 = -1;

  if (k_start == -1) {
    vid0 = triangles[tid][0];
    vid1 = triangles[tid][1];
    vid2 = triangles[tid][2];
  } else {
    vid0 = triangles[tid][(k_start + 1) % 3];
    vid1 = triangles[tid][(k_start + 2) % 3];
    vid2 = triangles[tid][k_start];
  }
  auto Tp = tranformation_matrix(ep[0], ep[1], ep[2],
                                 eval_position(triangles, positions, p));

  auto p0p = switch_reference_frame(Tp, positions[vid0]);
  auto p1p = switch_reference_frame(Tp, positions[vid1]);
  auto p2p = switch_reference_frame(Tp, positions[vid2]);
  auto verts = vector<vec2f>{{p0p.x, p0p.y}, {p1p.x, p1p.y}, {p2p.x, p2p.y}};
  for (auto k = 0; k < 3; ++k) {

    auto [t0, t1] = intersect(dir, verts[k], verts[(k + 1) % 3]);

    if (t0 > 0 && t1 >= -1e-4 && t1 <= 1 + 1e-4) {
      auto result_k = (k_start == -1) ? k : (k_start + 1 + k) % 3;
      return {result_k, clamp(t1, 0.f, 1.f)};
    }
  }
  return {-1, -1};
}
std::tuple<int, float> straightest_path_in_tri(const shape_data &data,
                                               const mesh_point &p,
                                               const vec2f &dir,
                                               const int k_start) {
  auto tid = p.face;
  auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  auto vid0 = -1;
  auto vid1 = -1;
  auto vid2 = -1;

  if (k_start == -1) {
    vid0 = data.triangles[tid][0];
    vid1 = data.triangles[tid][1];
    vid2 = data.triangles[tid][2];
  } else {
    vid0 = data.triangles[tid][(k_start + 1) % 3];
    vid1 = data.triangles[tid][(k_start + 2) % 3];
    vid2 = data.triangles[tid][k_start];
  }
  auto Tp = tranformation_matrix(
      ep[0], ep[1], ep[2], eval_position(data.triangles, data.positions, p));

  auto p0p = switch_reference_frame(Tp, data.positions[vid0]);
  auto p1p = switch_reference_frame(Tp, data.positions[vid1]);
  auto p2p = switch_reference_frame(Tp, data.positions[vid2]);
  auto verts = vector<vec2f>{{p0p.x, p0p.y}, {p1p.x, p1p.y}, {p2p.x, p2p.y}};
  for (auto k = 0; k < 3; ++k) {

    auto [t0, t1] = intersect(dir, verts[k], verts[(k + 1) % 3]);

    if (t0 > 0 && t1 >= -1e-4 && t1 <= 1 + 1e-4) {
      auto result_k = (k_start == -1) ? k : (k_start + 1 + k) % 3;
      return {result_k, clamp(t1, 0.f, 1.f)};
    }
    // if (k == 2)
    //   std::cerr << "Error in tracing within tri";
  }
  return {-1, -1};
}
std::tuple<vec2f, mesh_point, int> polthier_condition_at_vert(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3i> &adjacencies, const vector<float> &total_angles,
    const int vid, const int tid, const vec3f &dir) {
  auto k = find(triangles[tid], vid);
  auto theta = 0.5 * total_angles[vid];
  auto vert = positions[vid];
  auto v = positions[triangles[tid][(k + 2) % 3]] - vert;
  auto acc = angle(v, dir);
  auto prev = acc;
  auto curr_tid = tid;
  auto count = 0;
  while (acc < theta) {
    ++count;
    curr_tid = adjacencies[curr_tid][(k + 2) % 3];
    if (curr_tid == -1)
      return {zero2f, {-1, zero2f}, -1};
    k = find(triangles[curr_tid], vid);
    prev = acc;
    acc += angle(positions[triangles[curr_tid][(k + 1) % 3]] - vert,
                 positions[triangles[curr_tid][(k + 2) % 3]] - vert);
  }
  auto vid1 = triangles[curr_tid][(k + 1) % 3];
  auto vid2 = triangles[curr_tid][(k + 2) % 3];
  auto offset = theta - prev;
  auto l = length(positions[vid1] - vert);
  auto phi = angle(vert - positions[vid1], positions[vid2] - positions[vid1]);
  auto x = l * std::sin(offset) / std::sin(pif - phi - offset);
  auto alpha = x / length(positions[vid2] - positions[vid1]);
  auto flat_tid = init_flat_triangle(positions, triangles[curr_tid], 0);

  auto q = (1 - alpha) * flat_tid[(k + 1) % 3] + alpha * flat_tid[(k + 2) % 3];

  auto new_dir = q - flat_tid[k];

  return {new_dir, point_from_vert(triangles, vid, curr_tid), k};
}
std::tuple<vec2f, mesh_point, int>
polthier_condition_at_vert(const shape_data &data,
                           const shape_geometry &geometry, const int vid,
                           const int tid, const vec3f &dir) {
  return polthier_condition_at_vert(data.triangles, data.positions,
                                    geometry.adjacencies, geometry.total_angles,
                                    vid, tid, dir);
}
std::tuple<vector<vec3f>, mesh_point> polthier_straightest_geodesic(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3i> &adjacencies, const vector<vector<int>> &v2t,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &p, const vec2f &dir, const float &len,
    const vector<bool> &boundaries) {
  auto result = vector<vec3f>{};
  auto directions = vector<vec3f>{};
  auto accumulated = 0.f;
  auto curr_tid = p.face;
  auto curr_p = p;
  auto prev_p = mesh_point{};
  auto curr_dir = dir;
  auto ep = get_basis_at_p(triangles, positions, curr_tid);

  result.push_back(eval_position(triangles, positions, p));
  directions.push_back(ep[0] * dir.x + ep[1] * dir.y);

  Eigen::Matrix4d T_l = tranformation_matrix(ep[0], ep[1], ep[2],
                                             positions[triangles[curr_tid][0]]);
  auto k_start = -1;
  auto [is_vert, kv] = point_is_vert(p);
  auto [is_edge, ke] = point_is_edge(p);
  if (is_vert) {
    k_start = kv;
  } else if (is_edge)
    k_start = ke;

  auto count = 0;
  while (accumulated < len) {
    ++count;
    auto [k, t1] = straightest_path_in_tri(triangles, positions, curr_p,
                                           curr_dir, k_start);

    auto new_bary = zero3f;
    auto point_on_edge = mesh_point{};
    if (k != -1) {
      new_bary[k] = 1 - t1;
      new_bary[(k + 1) % 3] = t1;
      point_on_edge = mesh_point{curr_tid, vec2f{new_bary.y, new_bary.z}};
    } else {
      std::tie(is_edge, ke) = point_is_edge(curr_p, 5e-3);
      std::tie(is_vert, kv) = point_is_vert(curr_p, 5e-3);
      auto bary = get_bary(curr_p.uv);
      if (is_edge) {
        k = ke;
        t1 = bary[(k + 1) % 3];
        point_on_edge = curr_p;
      } else if (is_vert) {
        auto bary3 = zero3f;
        bary3[kv] = 1;
        point_on_edge = {curr_p.face, {bary3.y, bary3.z}};
      } else {
        std::cout << "Error!This should not happen" << std::endl;
        return {result, curr_p};
      }
    }

    std::tie(is_vert, kv) = point_is_vert(point_on_edge);

    if (is_vert) {
      auto vid = triangles[curr_tid][kv];
      if (angles[vid].size() == 0)
        return {result, curr_p};

      accumulated +=
          length(eval_position(triangles, positions, curr_p) - positions[vid]);
      auto dir3d = normalize(eval_position(triangles, positions, curr_p) -
                             positions[vid]);
      prev_p = curr_p;
      std::tie(curr_dir, curr_p, k_start) =
          polthier_condition_at_vert(triangles, positions, adjacencies,
                                     total_angles, vid, curr_tid, dir3d);
      curr_tid = curr_p.face;
      if (boundaries.size() > 0 && boundaries[curr_tid]) {
        return {result, curr_p};
      }
      if (curr_tid == -1)
        return {result, curr_p};

      ep = get_basis_at_p(triangles, positions, curr_tid);
      T_l = tranformation_matrix(ep[0], ep[1], ep[2],
                                 positions[triangles[curr_tid][0]]);

    } else {
      auto adj = adjacencies[curr_tid][k];
      if (adj == -1)
        return {result, curr_p};
      auto h = find(adjacencies[adj], curr_tid);

      new_bary = zero3f;
      new_bary[h] = t1;
      new_bary[(h + 1) % 3] = 1 - t1;

      prev_p = curr_p;
      curr_p = mesh_point{adj, vec2f{new_bary.y, new_bary.z}};
      if (boundaries.size() > 0 && boundaries[adj]) {
        return {result, curr_p};
      }
      accumulated += length(eval_position(triangles, positions, curr_p) -
                            eval_position(triangles, positions, prev_p));

      auto T = switch_reference_frame(triangles, positions, adj, curr_tid);
      curr_dir = switch_reference_frame_vector(T, curr_dir);

      curr_tid = adj;
      k_start = h;
    }

    result.push_back(eval_position(triangles, positions, curr_p));
    if (boundaries.size() > 0 && boundaries[curr_tid]) {
      return {result, curr_p};
    }
  }

  auto excess = accumulated - len;
  auto prev = result.rbegin()[1];
  auto last = result.back();
  auto alpha = excess / length(last - prev);
  auto pos = alpha * prev + (1 - alpha) * last;

  auto [inside, bary] =
      point_in_triangle(triangles, positions, prev_p.face, pos);
  if (!inside)
    std::cout << "error!This point should be in the triangle" << std::endl;

  result.pop_back();
  result.push_back(pos);

  return {result, {prev_p.face, bary}};
}

float energy_at_point(const vector<vec3i> &triangles,
                      const vector<vector<float>> &f, const vector<float> &w,
                      const mesh_point &p) {
  auto tid = p.face;
  auto en = zero3f;

  for (auto i = 0; i < 3; ++i) {
    auto vid = triangles[tid][i];
    for (auto j = 0; j < w.size(); ++j) {
      en[i] += w[j] * f[j][vid];
    }
  }

  return (1 - p.uv.x - p.uv.y) * en[0] + p.uv.x * en[1] + p.uv.y * en[2];
}
mesh_point polthier_straightest_geodesic_for_Newton_method(
    const shape_data &data, const shape_geometry &geometry, const mesh_point &p,
    const vector<vector<float>> &f, const vector<float> &w, const vec2f &dir,
    const float &len, const vector<bool> &boundaries) {
  auto line = vector<mesh_point>{};
  auto directions = vector<vec3f>{};
  auto accumulated = 0.f;
  auto curr_tid = p.face;
  auto curr_p = p;
  auto prev_p = mesh_point{};
  auto curr_dir = dir;
  auto ep = get_basis_at_p(data.triangles, data.positions, curr_tid);

  line.push_back(p);
  directions.push_back(ep[0] * dir.x + ep[1] * dir.y);
  auto en = energy_at_point(data.triangles, f, w, p);
  auto entry = 0;
  Eigen::Matrix4d T_l = tranformation_matrix(
      ep[0], ep[1], ep[2], data.positions[data.triangles[curr_tid][0]]);
  auto k_start = -1;
  auto [is_vert, kv] = point_is_vert(p);
  auto [is_edge, ke] = point_is_edge(p);
  if (is_vert) {
    k_start = kv;

  } else if (is_edge) {
    k_start = ke;
    auto dir3d = ep[0] * dir.x + ep[1] * dir.y;
    auto edge = data.positions[data.triangles[curr_tid][(ke + 1) % 3]] -
                data.positions[data.triangles[curr_tid][ke]];
    auto n = tid_normal(data.triangles, data.positions, curr_tid);
    if (dot(cross(edge, dir3d), n) > 0) {
      auto adj = geometry.adjacencies[curr_tid][ke];
      auto h = find(geometry.adjacencies[adj], curr_tid);
      auto bary = get_bary(p.uv);
      auto new_bary = zero3f;
      new_bary[h] = bary[(ke + 1) % 3];
      new_bary[(h + 1) % 3] = bary[ke];

      curr_p = mesh_point{adj, vec2f{new_bary.y, new_bary.z}};

      auto T = switch_reference_frame(data, geometry, adj, curr_tid);
      curr_dir = switch_reference_frame_vector(T, curr_dir);
      curr_tid = adj;
      k_start = h;
    }
  }

  auto count = 0;
  while (accumulated < len) {
    ++count;
    auto [k, t1] = straightest_path_in_tri(data, curr_p, curr_dir, k_start);

    auto new_bary = zero3f;
    auto point_on_edge = mesh_point{};
    if (k != -1) {
      new_bary[k] = 1 - t1;
      new_bary[(k + 1) % 3] = t1;
      point_on_edge = mesh_point{curr_tid, vec2f{new_bary.y, new_bary.z}};
    } else {
      std::tie(is_edge, ke) = point_is_edge(curr_p, 5e-3);
      std::tie(is_vert, kv) = point_is_vert(curr_p, 5e-3);
      auto bary = get_bary(curr_p.uv);
      if (is_edge) {
        k = ke;
        t1 = bary[(k + 1) % 3];
        point_on_edge = curr_p;
      } else if (is_vert) {
        auto bary3 = zero3f;
        bary3[kv] = 1;
        point_on_edge = {curr_p.face, {bary3.y, bary3.z}};
      } else {
        std::cout << "Error!This should not happen" << std::endl;
        return line[entry];
      }
    }

    std::tie(is_vert, kv) = point_is_vert(point_on_edge);

    if (is_vert) {
      auto vid = data.triangles[curr_tid][kv];
      if (geometry.angles[vid].size() == 0)
        return line[entry];

      accumulated +=
          length(eval_position(data.triangles, data.positions, curr_p) -
                 data.positions[vid]);
      auto dir3d =
          normalize(eval_position(data.triangles, data.positions, curr_p) -
                    data.positions[vid]);
      prev_p = curr_p;
      std::tie(curr_dir, curr_p, k_start) =
          polthier_condition_at_vert(data, geometry, vid, curr_tid, dir3d);

      curr_tid = curr_p.face;
      if (boundaries.size() > 0 && boundaries[curr_tid]) {
        return line[entry];
      }
      if (curr_tid == -1)
        return line[entry];

      ep = get_basis_at_p(data.triangles, data.positions, curr_tid);
      T_l = tranformation_matrix(ep[0], ep[1], ep[2],
                                 data.positions[data.triangles[curr_tid][0]]);

    } else {
      auto adj = geometry.adjacencies[curr_tid][k];
      if (adj == -1)
        return line[entry];
      auto h = find(geometry.adjacencies[adj], curr_tid);

      new_bary = zero3f;
      new_bary[h] = t1;
      new_bary[(h + 1) % 3] = 1 - t1;

      prev_p = curr_p;
      curr_p = mesh_point{adj, vec2f{new_bary.y, new_bary.z}};
      if (boundaries.size() > 0 && boundaries[adj]) {
        return line[entry];
      }

      accumulated +=
          length(eval_position(data.triangles, data.positions, curr_p) -
                 eval_position(data.triangles, data.positions, prev_p));

      auto T = switch_reference_frame(data, geometry, adj, curr_tid);
      curr_dir = switch_reference_frame_vector(T, curr_dir);

      curr_tid = adj;
      k_start = h;
    }
    auto curr_en = energy_at_point(data.triangles, f, w, curr_p);
    if (curr_en < en) {
      en = curr_en;
      entry = count;
    }
    line.push_back(curr_p);
    if (boundaries.size() > 0 && boundaries[curr_tid]) {
      return line[entry];
    }
  }

  return line[entry];
}

vec3f flip_bary_to_adjacent_tri(const vector<vec3i> &adjacencies,
                                const int tid0, const int tid1,
                                const vec3f &bary) {
  if (tid0 == tid1)
    return bary;
  auto new_bary = zero3f;
  auto k1 = find(adjacencies[tid1], tid0);
  auto k0 = find(adjacencies[tid0], tid1);
  if (k1 == -1) {
    std::cout << "Error, faces are not adjacent" << std::endl;
    return zero3f;
  }
  new_bary[k1] = bary[(k0 + 1) % 3];
  new_bary[(k1 + 1) % 3] = bary[k0];
  new_bary[(k1 + 2) % 3] = bary[(k0 + 2) % 3];

  return new_bary;
}
std::tuple<vector<int>, mesh_point, mesh_point>
handle_short_strips(const vector<vec3i> &triangles,
                    const vector<vec3f> &positions, const vector<int> &strip,
                    const mesh_point &start, const mesh_point &end) {
  if (strip.size() == 1) {
    return {strip, start, end};
  } else if (strip.size() == 2) {
    auto [inside, b2f] =
        point_in_triangle(triangles, positions, start.face,
                          eval_position(triangles, positions, end));
    if (inside) {
      auto new_end = mesh_point{start.face, b2f};
      return {{start.face}, start, new_end};
    }
    std::tie(inside, b2f) =
        point_in_triangle(triangles, positions, end.face,
                          eval_position(triangles, positions, start));

    if (inside) {
      auto new_start = mesh_point{end.face, b2f};
      return {{end.face}, new_start, end};
    }

    return {strip, start, end};
  }
  return {{-1}, {}, {}};
}
std::tuple<vector<int>, mesh_point, mesh_point>
cleaned_strip(const vector<vec3i> &triangles, const vector<vec3f> &positions,
              const vector<vec3i> &adjacencies, const vector<int> &strip,
              const mesh_point &start, const mesh_point &end) {
  vector<int> cleaned = strip;

  auto start_entry = 0, end_entry = (int)strip.size() - 1;
  auto b3f = zero3f;
  auto new_start = start;
  auto new_end = end;
  auto [is_vert, kv] = point_is_vert(end);
  auto [is_edge, ke] = point_is_edge(end);
  if (strip.size() <= 2)
    return handle_short_strips(triangles, positions, strip, start, end);
  // Erasing from the bottom
  if (is_vert) {
    auto vid = triangles[end.face][kv];
    auto curr_tid = strip[end_entry - 1];
    kv = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.pop_back();
      --end_entry;
      if (end_entry == 1)
        break;
      // see comment below
      auto curr_tid = strip[end_entry - 1];
      kv = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned.back()], vid);
    assert(kv != -1);
    b3f[kv] = 1;
    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}}; // updating end
  } else if (is_edge) {
    if (end.face != strip.back()) {
      assert(adjacencies[end.face][ke] == strip.back());

      if (end.face == strip[end_entry - 1])
        cleaned.pop_back();
    } else if (adjacencies[end.face][ke] == strip[end_entry - 1])
      cleaned.pop_back();

    b3f = flip_bary_to_adjacent_tri(adjacencies, end.face, cleaned.back(),
                                    get_bary(end.uv));

    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}}; // updating end
  }
  std::tie(is_vert, kv) = point_is_vert(start);
  std::tie(is_edge, ke) = point_is_vert(start);

  if (is_vert) {
    auto vid = triangles[start.face][kv];
    auto curr_tid = strip[start_entry + 1];
    kv = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.erase(cleaned.begin());
      ++start_entry;
      if (start_entry > end_entry - 1)
        break;
      auto curr_tid = strip[start_entry + 1];
      kv = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned[0]], vid);
    assert(kv != -1);
    b3f = zero3f;
    b3f[kv] = 1;
    new_start = mesh_point{cleaned[0], vec2f{b3f.y, b3f.z}}; // udpdating start

  } else if (is_edge) {
    if (start.face != strip[0]) {
      assert(adjacencies[start.face][ke] == strip[0]);
      if (start.face == strip[1])
        cleaned.erase(cleaned.begin());
    } else if (adjacencies[start.face][ke] == strip[1]) {
      cleaned.erase(cleaned.begin());
    }
    b3f = flip_bary_to_adjacent_tri(adjacencies, start.face, cleaned[0],
                                    get_bary(start.uv));
    new_start = {cleaned[0], vec2f{b3f.y, b3f.z}}; // updating start
  }
  return {cleaned, new_start, new_end};
}
template <typename Update, typename Stop, typename Exit>
void search_strip(vector<float> &field, vector<bool> &in_queue,
                  const dual_geodesic_solver &solver,
                  const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int start, int end,
                  Update &&update, Stop &&stop, Exit &&exit) {
  auto destination_pos =
      eval_position(triangles, positions, {end, {1.0f / 3, 1.0f / 3}});

  auto estimate_dist = [&](int face) {
    auto p = eval_position(triangles, positions, {face, {1.0f / 3, 1.0f / 3}});
    return length(p - destination_pos);
  };
  field[start] = estimate_dist(start);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

  while (!queue.empty()) {
    auto node = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight)
        break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node))
      break;
    if (stop(node))
      continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1)
        continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance)
        continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      if (update(node, neighbor, new_distance))
        return;
    }
  }
}
vector<int> compute_strip_tlv(const shape_data &data,
                              const dual_geodesic_solver &solver, int start,
                              int end) {
  if (start == end)
    return {start};

  thread_local static auto parents = vector<int>{};
  thread_local static auto field = vector<float>{};
  thread_local static auto in_queue = vector<bool>{};

  if (parents.size() != solver.graph.size()) {
    parents.assign(solver.graph.size(), -1);
    field.assign(solver.graph.size(), flt_max);
    in_queue.assign(solver.graph.size(), false);
  }

  // initialize once for all and sparsely cleanup at the end of every solve
  auto visited = vector<int>{start};
  auto sources = vector<int>{start};
  auto update = [&visited, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(field, in_queue, solver, data.triangles, data.positions, start,
               end, update, stop, exit);

  // extract_strip
  auto strip = vector<int>{};
  auto node = end;
  strip.reserve((int)sqrt(parents.size()));
  while (node != -1) {
    assert(find(strip, node) != 1);
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  for (auto &v : visited) {
    parents[v] = -1;
    field[v] = flt_max;
    in_queue[v] = false;
  }
  // assert(check_strip(geometry.adjacencies, strip));
  return strip;
}

geodesic_path compute_geodesic_path(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const dual_geodesic_solver &solver,
                                    const mesh_point &start,
                                    const mesh_point &end) {
  // profile_function();

  vector<int> parents;
  auto strip = compute_strip_tlv(data, solver, end.face, start.face);

  auto path = geodesic_path{};
  auto [cleaned, new_start, new_end] = cleaned_strip(
      data.triangles, data.positions, geometry.adjacencies, strip, start, end);

  if (new_start.face == new_end.face) {
    path.start = new_start;
    path.end = new_end;
    path.strip = {new_start.face};
    return path;
  }

  path = shortest_path(data.triangles, data.positions, geometry.adjacencies,
                       new_start, new_end, cleaned);
  return path;
}

Eigen::SparseMatrix<double> AGS_matrix(const shape_data &data,
                                       const shape_geometry &geometry) {
  Eigen::SparseMatrix<double> A(3 * data.positions.size(),
                                data.positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T> entries;
  vector<pair<int, vec3f>> pesi;
  for (int j = 0; j < data.positions.size(); ++j) {
    if (geometry.angles[j].size() == 0)
      continue;
    vector<int> poly_nbr = geometry.v2t[j];

    auto star = geometry.v2t[j];
    auto s = star.size();
    float wgt = 0.0;

    pesi.clear();

    if (s > 0) {
      for (int i = 0; i < s; ++i) {
        auto tid = star[i];
        auto n = tid_normal(data.triangles, data.positions, tid);
        auto area = tid_area(data, tid);
        wgt += area;
        for (auto off = 0; off < 3; ++off) {
          auto contribute = PCE_contribute(data.triangles, data.positions, tid,
                                           n, 2 * area, off);
          contribute *= area;
          parallel_transp(geometry.angles, geometry.total_angles,
                          data.triangles, data.positions, geometry.adjacencies,
                          geometry.v2t, contribute, data.normals, tid, j, T2V);
          pesi.push_back(
              make_pair(data.triangles[tid][(off + 1) % 3], contribute));
        }
        auto opp = opposite_face(data.triangles, geometry.adjacencies, tid, j);
        n = tid_normal(data.triangles, data.positions, opp);
        area = tid_area(data, opp);
        wgt += area;
        for (auto off = 0; off < 3; ++off) {
          auto contribute = PCE_contribute(data.triangles, data.positions, opp,
                                           n, 2 * area, off);
          contribute *= area;
          parallel_transp(geometry.angles, geometry.total_angles,
                          data.triangles, data.positions, geometry.adjacencies,
                          geometry.v2t, contribute, data.normals, opp, tid,
                          T2T);
          parallel_transp(geometry.angles, geometry.total_angles,
                          data.triangles, data.positions, geometry.adjacencies,
                          geometry.v2t, contribute, data.normals, tid, j, T2V);
          pesi.push_back(
              make_pair(data.triangles[opp][(off + 1) % 3], contribute));
        }
      }
      int row = j * 3;
      for (int h = 0; h < pesi.size(); ++h) {
        entries.push_back(T(row, pesi[h].first, pesi[h].second.x / wgt));
        entries.push_back(T(row + 1, pesi[h].first, pesi[h].second.y / wgt));
        entries.push_back(T(row + 2, pesi[h].first, pesi[h].second.z / wgt));
      }
    }
  }
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

pair<float, int> optimal_sample(const vector<vec3i> &triangles,
                                const geodesic_path &path, const int vid) {
  auto lerp = path.lerps[0];
  auto tid = path.strip[0];
  auto entry = 0;
  while (lerp == 0 || lerp == 1) {
    ++entry;
    if (path.lerps.size() == entry)
      return {lerp, tid};
    if (find(triangles[path.strip[entry]], vid) != -1) {
      lerp = path.lerps[entry];
      tid = path.strip[entry];
    } else
      return {lerp, tid};
  }

  return {lerp, tid};
}
geodesic_solver extended_solver(const shape_data &data,
                                const dual_geodesic_solver &dual_solver,
                                shape_geometry &geometry, const int k) {
  time_function();
  auto old_solver = make_geodesic_solver(data.triangles, data.positions,
                                         geometry.adjacencies, geometry.v2t);
  geometry.angles = compute_angles_wo_opposite(
      data.triangles, data.positions, geometry.adjacencies, geometry.v2t,
      geometry.total_angles);
  auto avg_valence = 0;
  geodesic_solver solver;
  solver.graph.resize(data.positions.size());

  for (auto i = 0; i < data.positions.size(); ++i) {
    if (geometry.v2t[i].size() == 0)
      continue;
    auto neighborhood = k_ring(data.triangles, geometry.v2t, i, k);
    auto valence = neighborhood.size();
    auto lengths = unordered_map<int, float>{};
    auto source = point_from_vert(data.triangles, geometry.v2t, i);
    for (auto j = 0; j < neighborhood.size(); ++j) {
      auto curr_point =
          point_from_vert(data.triangles, geometry.v2t, neighborhood[j]);
      auto path = compute_geodesic_path(data, geometry, dual_solver, source,
                                        curr_point);
      if (path.lerps.size() > 0) {

        lengths[neighborhood[j]] = path_length(
            path, data.triangles, data.positions, geometry.adjacencies);
      } else {
        auto entry = node_is_adjacent(old_solver, i, neighborhood[j]);
        if (entry == -1) {
          lengths[neighborhood[j]] = flt_max;
          continue;
        }

        lengths[neighborhood[j]] = old_solver.graph[i][entry].length;
      }
    }
    avg_valence += valence;

    solver.graph[i].resize(valence);
    for (auto j = 0; j < valence; ++j) {
      solver.graph[i][j] = {neighborhood[j], lengths.at(neighborhood[j])};
    }
  }
  return solver;
}
