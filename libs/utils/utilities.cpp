#include "utilities.h"
#include <deque>

std::tuple<mesh_point, int>
check_star_around_vert_pce(const shape_data &data,
                           const shape_geometry &geometry, const int vid,
                           const mesh_point &target, const vector<vec3f> &grad,
                           const int tid, const int next_tid) {
  auto star = geometry.v2t[vid];
  for (auto i = 0; i < star.size(); ++i) {
    auto curr_tid = star[i];
    if (curr_tid == target.face) {
      return {target, -1};
    }
    if (curr_tid == next_tid || curr_tid == tid)
      continue;
    auto k = find(data.triangles[curr_tid], vid);
    auto p = point_from_vert(data.triangles, vid, curr_tid);
    auto ep = get_basis_at_p(data.triangles, data.positions, p.face);
    auto T_l = tranformation_matrix(ep[0], ep[1], ep[2], data.positions[vid]);
    auto curr_dir = switch_reference_frame_vector(T_l, grad[p.face]);
    auto p0 = switch_reference_frame(
        T_l, data.positions[data.triangles[curr_tid][(k + 1) % 3]]);
    auto p1 = switch_reference_frame(
        T_l, data.positions[data.triangles[curr_tid][(k + 2) % 3]]);
    auto [t0, t1] =
        intersect({curr_dir.x, curr_dir.y}, {p0.x, p0.y}, {p1.x, p1.y});

    if (t0 > 0 && t1 >= 0 && t1 <= 1) {
      auto adj = geometry.adjacencies[curr_tid][(k + 1) % 3];
      auto h = find(geometry.adjacencies[adj], curr_tid);
      clamp(t1, 0.f, 1.f);
      auto bary = zero3f;
      bary[h] = t1;
      bary[(h + 1) % 3] = 1 - t1;
      return {{adj, {bary.y, bary.z}}, h};
    }
  }
  return {{-1, zero2f}, -1};
}
std::tuple<int, int, int> best_next_vid(const shape_data &data,
                                        const shape_geometry &geometry,
                                        const int vid, const int vid0) {
  auto nbr = one_ring(data.triangles, geometry.v2t, vid);
  float phi = -1;
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i] == vid0) {
      phi = geometry.angles[vid][i];
      break;
    }
  }

  if (phi == -1) {
    std::cout << "We messed up with PCE" << std::endl;
  }

  phi = std::fmod(phi + pif, 2 * pif);
  auto new_vid = -1;
  auto prev_tid = -1;
  auto next_tid = -1;
  auto offset = flt_max;
  auto star = geometry.v2t[vid];
  auto s = star.size();
  for (auto i = 0; i < nbr.size(); ++i) {
    auto diff = std::abs(phi - geometry.angles[vid][i]);
    if (diff < offset) {
      offset = diff;
      new_vid = nbr[i];
      prev_tid = star[(s + i - 1) % s];
      next_tid = star[i];
    }
  }

  return {new_vid, prev_tid, next_tid};
}
std::tuple<int, int> next_vid(const shape_data &data,
                              const shape_geometry &geometry,
                              const vector<vec3f> &grad, const int vid) {
  auto star = geometry.v2t[vid];
  auto vert = data.positions[vid];
  auto min_angle = flt_max;
  auto best_tid = -1;
  auto best_vid = -1;
  for (auto tid : star) {

    auto k = find(data.triangles[tid], vid);
    auto v1 = data.triangles[tid][(k + 1) % 3];
    auto e1 = data.positions[v1] - vert;
    auto phi = angle(e1, grad[tid]);
    if (phi < min_angle) {
      min_angle = phi;
      best_tid = tid;
      best_vid = v1;
    }
  }

  return {best_vid, best_tid};
}
std::tuple<mesh_point, int>
handle_vert_pce_tracing(const shape_data &data, const shape_geometry &geometry,
                        const vector<vec3f> &grad, const mesh_point &curr_p,
                        const int next_tid, const mesh_point &target,
                        vector<vec3f> &result) {
  auto curr_dir = grad[curr_p.face];
  auto k = find(geometry.adjacencies[curr_p.face], next_tid);
  auto edge = data.positions[data.triangles[curr_p.face][(k + 1) % 3]] -
              data.positions[data.triangles[curr_p.face][k]];
  auto vid0 = -1;
  auto vid1 = -1;
  if (dot(curr_dir, edge) > 0) {
    vid0 = data.triangles[curr_p.face][(k + 1) % 3];
    vid1 = data.triangles[curr_p.face][k];
  } else {
    vid0 = data.triangles[curr_p.face][k];
    vid1 = data.triangles[curr_p.face][(k + 1) % 3];
  }

  result.push_back(data.positions[vid0]);
  auto prev_tid = curr_p.face;
  auto curr_tid = next_tid;
  while (true) {
    auto [new_point, k] = check_star_around_vert_pce(
        data, geometry, vid0, target, grad, prev_tid, curr_tid);
    if (new_point.face != -1 && k == -1) {
      result.push_back(eval_position(data.triangles, data.positions, target));
      return {{-1, zero2f}, -1};
    } else if (new_point.face != -1)
      return {new_point, k};
    else {
      auto tmp = vid0;
      std::tie(vid0, prev_tid, curr_tid) =
          best_next_vid(data, geometry, vid0, vid1);
      result.push_back(data.positions[tmp]);
      vid1 = tmp;
    }
  }
}
std::tuple<mesh_point, int>
handle_vert_pce_tracing(const shape_data &data, const shape_geometry &geometry,
                        const vector<vec3f> &grad, const int prev_tid,
                        const int &vid, const mesh_point &target,
                        vector<vec3f> &result) {

  auto curr_tid = prev_tid;
  auto curr_vid = vid;
  while (true) {
    auto [new_point, k] = check_star_around_vert_pce(
        data, geometry, curr_vid, target, grad, curr_tid, -1);
    if (new_point.face != -1 && k == -1) {
      result.push_back(eval_position(data.triangles, data.positions, target));
      return {{-1, zero2f}, -1};
    } else if (new_point.face != -1) {
      result.push_back(
          eval_position(data.triangles, data.positions, new_point));
      return {new_point, k};
    }

    else {
      std::tie(curr_vid, curr_tid) = next_vid(data, geometry, grad, curr_vid);
    }
  }
}
bool check_next_tid_pce(const shape_data &data, const vector<vec3f> &grad,
                        const int tid, const int next_tid, const int k) {
  auto vid0 = data.triangles[tid][k];
  auto vid1 = data.triangles[tid][(k + 1) % 3];
  auto edge = normalize(data.positions[vid0] - data.positions[vid1]);
  auto n = tid_normal(data.triangles, data.positions, next_tid);
  if (dot(cross(edge, grad[next_tid]), n) < 0)
    return false;

  return true;
}

vector<vec3f>
shortest_path_using_pce(const shape_data &data, const shape_geometry &geometry,
                        const vector<vec3f> &grad, const mesh_point &source,
                        const mesh_point &target, const int max_it = INT_MAX) {
  auto curr_p = source;
  auto k_start = -1;
  auto result = vector<vec3f>{};
  auto count = 0;
  result.push_back(eval_position(data.triangles, data.positions, curr_p));
  while (curr_p.face != target.face) {
    // while (count < max_it) {
    ++count;

    auto ep = get_basis_at_p(data.triangles, data.positions, curr_p.face);
    auto T_l = tranformation_matrix(
        ep[0], ep[1], ep[2], data.positions[data.triangles[curr_p.face][0]]);
    auto curr_dir =
        switch_reference_frame_vector(T_l, normalize(grad[curr_p.face]));
    auto [k, t] = straightest_path_in_tri(data, curr_p,
                                          {curr_dir.x, curr_dir.y}, k_start);
    if (k == -1) {
      auto [is_vert, kv] = point_is_vert(curr_p);
      auto [is_edge, ke] = point_is_edge(curr_p);

      if (!is_vert && !is_edge) {
        std::cout << "PCE ti odio" << std::endl;
        return result;
      }

      if (is_vert) {
        std::tie(curr_p, k_start) = handle_vert_pce_tracing(
            data, geometry, grad, curr_p.face, data.triangles[curr_p.face][kv],
            target, result);
      } else {
        auto curr_adj = geometry.adjacencies[curr_p.face][ke];
        std::tie(curr_p, k_start) = handle_vert_pce_tracing(
            data, geometry, grad, curr_p, curr_adj, target, result);
      }

    } else {
      auto adj = geometry.adjacencies[curr_p.face][k];
      if (count > 1 && !check_next_tid_pce(data, grad, curr_p.face, adj, k)) {
        std::tie(curr_p, k_start) = handle_vert_pce_tracing(
            data, geometry, grad, curr_p, adj, target, result);
        if (curr_p.face == -1)
          return result;
      } else {
        auto h = find(geometry.adjacencies[adj], curr_p.face);
        auto bary = zero3f;
        bary[h] = t;
        bary[(h + 1) % 3] = 1 - t;
        curr_p = {adj, {bary.y, bary.z}};
        k_start = h;
      }
      if (curr_p.face == -1) {
        std::cout << "we messed up" << std::endl;

        return result;
      }

      result.push_back(eval_position(data.triangles, data.positions, curr_p));
    }
  }
  result.push_back(eval_position(data.triangles, data.positions, target));
  return result;
}

int bin_coefficient(const int &n, const int &k) {
  int bin = 1;

  for (int i = 0; i < k; ++i) {
    bin *= (n - i);
    bin /= (i + 1);
  }

  return bin;
}

vec3f gradient_blending(const vector<vector<vec3f>> &gradients,
                        const vector<float> &weights, const int vid) {
  auto v = zero3f;

  for (int i = 0; i < weights.size(); ++i) {
    auto w = gradients[i][vid] * weights[i];
    v += w;
  }

  return v;
}
pair<bool, int> bary_for_minimum_hessian(const vec3f &i, const vec3f &j,
                                         const vec3f &k, float &alpha,
                                         float &beta, float &gamma) {
  auto g0_norm = length_squared(i);
  auto g1_norm = length_squared(j);
  auto g2_norm = length_squared(k);
  auto g01 = dot(i, j);
  auto g02 = dot(i, k);
  auto g12 = dot(j, k);
  auto l01 = length_squared(i - j);
  auto l02 = length_squared(i - k);

  auto A = g0_norm + g12 - g02 - g01;
  Eigen::Matrix3d A3;
  A3 << g0_norm, g01, g02, g01, g1_norm, g12, g02, g12, g2_norm;

  auto det = l02 * l01 - pow(A, 2);

  auto det_A_3 = g0_norm * length_squared(cross(j, k)) +
                 g1_norm * length_squared(cross(i, k)) -
                 g2_norm * length_squared(cross(i, j));
  auto contribute = dot(cross(i, k), cross(j, i)) +
                    dot(cross(j, k), cross(i, j)) +
                    dot(cross(k, j), cross(i, k));

  auto cond_A = (l01 + l02 + sqrt(pow(l01 - l02, 2) - 4 * pow(A, 2))) /
                (l01 + l02 - sqrt(pow(l01 - l02, 2) - 4 * pow(A, 2)));
  cond_A = std::abs(cond_A);

  if (det > 1e-10) {

    gamma = ((g0_norm - g02) * l01 + A * (g01 - g0_norm)) / det;
    beta = (g0_norm - g01 - gamma * A) / l01;
    alpha = 1 - beta - gamma;
    if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 &&
        gamma <= 1) {
      return {true, -1};

    } else
      return {false, -1};
  } else

    return {false, -1};
}
pair<bool, int> tri_contains_min_hessian(const vec3f &g0, const vec3f &g1,
                                         const vec3f &g2, vec3f &bary) {

  return bary_for_minimum_hessian(g0, g1, g2, bary.x, bary.y, bary.z);
}

std::pair<vec2f, bool>
maximum_descent_direction(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions, const int tid,
                          const Eigen::Matrix4d &T, const vec3f &grad,
                          const vec3f &gx, const vec3f &gy, const vec3f &gz) {

  auto vid1 = triangles[tid].y;
  auto vid2 = triangles[tid].z;
  auto v0 = switch_reference_frame_vector(T, gx);
  auto v1 = switch_reference_frame_vector(T, gy);
  auto v2 = switch_reference_frame_vector(T, gz);
  auto v = switch_reference_frame_vector(T, grad);
  auto p1p = switch_reference_frame(T, positions[vid1]);
  auto p2p = switch_reference_frame(T, positions[vid2]);

  auto ax = (v1.x - v0.x) / p1p.x;
  auto ay = (v1.y - v0.y) / p1p.x;
  auto bx = (v2.x - v0.x - ax * p2p.x) / p2p.y;
  auto by = (v2.y - v0.y - ay * p2p.x) / p2p.y;
  Eigen::Matrix2d H;
  H(0, 0) = ax;
  H(0, 1) = (bx + ay) / 2;
  H(1, 0) = (ay + bx) / 2;
  H(1, 1) = by;
  Eigen::LLT<Eigen::Matrix2d> lltOfH(H);

  if (lltOfH.info() == Eigen::NumericalIssue) {
    std::cout << "is not a descent direction" << std::endl;
    return {zero2f, false};
  }

  Eigen::Vector2d b;
  b << -v.x, -v.y;
  auto x = lltOfH.solve(b);

  return {vec2f{(float)x(0), (float)x(1)}, true};
}
std::tuple<vec3f, mesh_point>
Riemannian_Newton_method(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const vector<vector<vec3f>> &grads,
                         const vector<vector<float>> &fields,
                         const vector<float> &w, const mesh_point &start) {
  auto tid = start.face;
  auto tags = vector<int>(data.triangles.size(), 0);
  auto curr_sample = start;
  mesh_point next_sample = {-1, zero2f};
  auto gx = gradient_blending(grads, w, data.triangles[tid].x);
  auto gy = gradient_blending(grads, w, data.triangles[tid].y);
  auto gz = gradient_blending(grads, w, data.triangles[tid].z);
  auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  auto T = tranformation_matrix(ep[0], ep[1], ep[2],
                                data.positions[data.triangles[tid].x]);
  // tags[tid] += 1;
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, gx,
                  data.normals, data.triangles[tid].x, tid, V2T);
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, gy,
                  data.normals, data.triangles[tid].y, tid, V2T);
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, gz,
                  data.normals, data.triangles[tid].z, tid, V2T);
  auto grad = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
              curr_sample.uv.x * gy + curr_sample.uv.y * gz;
  if (length(grad) < 1e-14)
    return {eval_position(data.triangles, data.positions, curr_sample),
            curr_sample};
  auto [dir, is_descent_direction] = maximum_descent_direction(
      data.triangles, data.positions, tid, T, -grad, -gx, -gy, -gz);
  if (!is_descent_direction)
    return {eval_position(data.triangles, data.positions, curr_sample),
            curr_sample};
  auto bary = zero3f;
  auto [arrived, type] = tri_contains_min_hessian(gx, gy, gz, bary);
  if (arrived) {
    curr_sample = {tid, {bary.y, bary.z}};
    return {eval_position(data.triangles, data.positions, curr_sample),
            curr_sample};
  }

  auto count = 0;
  while (!arrived) {

    ++count;

    next_sample = polthier_straightest_geodesic_for_Newton_method(
        data, geometry, curr_sample, fields, w, normalize(dir), length(dir));

    if (length(eval_position(data.triangles, data.positions, next_sample) -
               eval_position(data.triangles, data.positions, curr_sample)) <
        1e-16) {

      return {zero3f, {-1, zero2f}};
    }
    auto [is_vert, kv] = point_is_vert(next_sample);
    auto [is_edge, ke] = point_is_edge(next_sample);
    if (is_vert) {
      auto vid = data.triangles[next_sample.face][kv];

      auto [tmp_dir, tmp_point, k] =
          polthier_condition_at_vert(data, geometry, vid, next_sample.face,
                                     gradient_blending(grads, w, vid));
      next_sample = tmp_point;
      tid = next_sample.face;
      gx = gradient_blending(grads, w, data.triangles[tid].x);
      gy = gradient_blending(grads, w, data.triangles[tid].y);
      gz = gradient_blending(grads, w, data.triangles[tid].z);
      ep = get_basis_at_p(data.triangles, data.positions, tid);
      T = tranformation_matrix(ep[0], ep[1], ep[2],
                               data.positions[data.triangles[tid].x]);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gx,
                      data.normals, data.triangles[tid].x, tid, V2T);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gy,
                      data.normals, data.triangles[tid].y, tid, V2T);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gz,
                      data.normals, data.triangles[tid].z, tid, V2T);
      grad = (1 - tmp_point.uv.x - tmp_point.uv.y) * gx + tmp_point.uv.x * gy +
             tmp_point.uv.y * gz;

      std::tie(arrived, type) = tri_contains_min_hessian(gx, gy, gz, bary);
      if (arrived) {
        curr_sample = {tid, {bary.y, bary.z}};
        return {eval_position(data.triangles, data.positions, curr_sample),
                curr_sample};
      }

      std::tie(dir, is_descent_direction) = maximum_descent_direction(
          data.triangles, data.positions, tid, T, -grad, -gx, -gy, -gz);

      if (!is_descent_direction)
        return {eval_position(data.triangles, data.positions, curr_sample),
                curr_sample};
    } else if (is_edge) {
      tid = next_sample.face;
      gx = gradient_blending(grads, w, data.triangles[tid].x);
      gy = gradient_blending(grads, w, data.triangles[tid].y);
      gz = gradient_blending(grads, w, data.triangles[tid].z);
      ep = get_basis_at_p(data.triangles, data.positions, tid);
      T = tranformation_matrix(ep[0], ep[1], ep[2],
                               data.positions[data.triangles[tid].x]);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gx,
                      data.normals, data.triangles[tid].x, tid, V2T);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gy,
                      data.normals, data.triangles[tid].y, tid, V2T);
      parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                      data.positions, geometry.adjacencies, geometry.v2t, gz,
                      data.normals, data.triangles[tid].z, tid, V2T);
      grad = (1 - next_sample.uv.x - next_sample.uv.y) * gx +
             next_sample.uv.x * gy + next_sample.uv.y * gz;
      std::tie(arrived, type) = tri_contains_min_hessian(gx, gy, gz, bary);
      if (arrived) {
        curr_sample = {tid, {bary.y, bary.z}};
        return {eval_position(data.triangles, data.positions, curr_sample),
                curr_sample};
      }

      std::tie(dir, is_descent_direction) = maximum_descent_direction(
          data.triangles, data.positions, tid, T, -grad, -gx, -gy, -gz);

      if (!is_descent_direction)
        return {zero3f, {-1, zero2f}};
    } else
      std::cout << "This should not happen" << std::endl;

    curr_sample = next_sample;
    if (length(grad) < 1e-14)
      return {eval_position(data.triangles, data.positions, curr_sample),
              curr_sample};
  }

  return {eval_position(data.triangles, data.positions, curr_sample),
          curr_sample};
}

// Bernstein polynomials
void bernstein_polynomials(const int &n, const float &t, vector<float> &w) {
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    float lambda =
        bin_coefficient(n - 1, i) * pow(t, i) * pow(1 - t, n - i - 1);
    w[i] = lambda;
  }
}
vector<float> compute_distance_field(const shape_data &data,
                                     const shape_geometry &geometry,
                                     const geodesic_solver &geo_solver,
                                     const mesh_point &source,
                                     const int solver) {
  if (solver == exact)
    return exact_geodesic_distance(data.triangles, data.positions, source);
  else
    return compute_geodesic_distances(geo_solver, data.triangles,
                                      data.positions, geometry.adjacencies,
                                      {source});
}
vector<vector<float>> fields_from_lndmrks(const shape_data &data,
                                          const shape_geometry &geometry,
                                          const geodesic_solver &geo_solver,
                                          const vector<mesh_point> &lndmrks,
                                          const int solver) {
  auto n = lndmrks.size();
  auto result = vector<vector<float>>(n);
  for (auto i = 0; i < n; ++i) {
    result[i] =
        compute_distance_field(data, geometry, geo_solver, lndmrks[i], solver);
    std::transform(result[i].begin(), result[i].end(), result[i].begin(),
                   [](float lambda) { return lambda * lambda; });
  }
  return result;
}
void update_fields_from_lndmrks(const shape_data &data,
                                const shape_geometry &geometry,
                                const geodesic_solver &geo_solver,
                                const vector<mesh_point> &lndmrks,
                                const int solver, const vector<bool> &moved,
                                vector<vector<float>> &f) {
  if (f.size() == 0) {
    f = fields_from_lndmrks(data, geometry, geo_solver, lndmrks, solver);
  } else if (moved.size() != f.size()) {
    for (auto i = 0; i < moved.size(); ++i) {
      if (i < f.size() && moved[i]) {
        f[i] = compute_distance_field(data, geometry, geo_solver, lndmrks[i],
                                      solver);
        std::transform(f[i].begin(), f[i].end(), f[i].begin(),
                       [](float lambda) { return lambda * lambda; });
      } else if (i >= f.size()) {
        f.push_back({});
        f[i] = compute_distance_field(data, geometry, geo_solver, lndmrks[i],
                                      solver);
        std::transform(f[i].begin(), f[i].end(), f[i].begin(),
                       [](float lambda) { return lambda * lambda; });
      }
    }

  } else {
    for (auto i = 0; i < moved.size(); ++i) {
      if (moved[i]) {
        f[i] = compute_distance_field(data, geometry, geo_solver, lndmrks[i],
                                      solver);
        std::transform(f[i].begin(), f[i].end(), f[i].begin(),
                       [](float lambda) { return lambda * lambda; });
      }
    }
  }
}
void update_grads_from_lndmrks(const shape_op &op,
                               const vector<vector<float>> &f,
                               const vector<bool> &moved,
                               vector<vector<vec3f>> &grds) {
  auto V = (int)f[0].size();

  if (grds.size() == 0) {
    grds = grads_from_lndmrks(op, f);
  } else if (moved.size() != grds.size()) {

    for (auto i = 0; i < moved.size(); ++i) {
      if (i < grds.size() && moved[i]) {
        grds[i] = AGS_grad(op.AGS_Grad, wrapper(f[i]));
      } else if (i >= grds.size()) {
        grds.push_back({});
        grds[i] = AGS_grad(op.AGS_Grad, wrapper(f[i]));
      }
    }
  } else {
    for (auto i = 0; i < moved.size(); ++i) {
      if (moved[i]) {
        grds[i] = AGS_grad(op.AGS_Grad, wrapper(f[i]));
      }
    }
  }
}

vector<vector<vec3f>> grads_from_lndmrks(const shape_op &op,
                                         const vector<vector<float>> &fields) {
  auto n = fields.size();
  auto result = vector<vector<vec3f>>(n);
  for (auto i = 0; i < n; ++i) {
    result[i] = AGS_grad(op.AGS_Grad, wrapper(fields[i]));
  }
  return result;
}

vector<vec3f> bézier_curve(const shape_data &data,
                           const shape_geometry &geometry, const shape_op &op,
                           const geodesic_solver &geo_solver,
                           const dual_geodesic_solver &solver,
                           const vector<mesh_point> &control_points,
                           const int k, const int type_of_solver) {
  int n = (int)control_points.size();

  vector<vector<float>> f = fields_from_lndmrks(data, geometry, geo_solver,
                                                control_points, type_of_solver);
  auto grds = grads_from_lndmrks(op, f);
  vector<float> w(n);
  bernstein_polynomials(n, 0.0, w);

  auto result = vector<mesh_point>((int)pow(2, k));

  auto seed = control_points[0];

  double step = 1 / pow(2, k);
  auto t = step;
  result[0] = seed;
  auto pos = zero3f;
  auto point = mesh_point{};
  bernstein_polynomials(n, t, w);
  auto count = 0;

  for (auto i = 1; i < pow(2, k); ++i) {

    std::tie(pos, point) =
        Riemannian_Newton_method(data, geometry, op, grds, f, w, seed);
    if (point.face != -1) {
      result[i] = point;
      seed = point;
    } else {
      ++count;
      result[i] = result[i - 1];
    }

    t += step;
    bernstein_polynomials(n, t, w);
  }

  result.push_back(control_points.back());

  auto polyline = vector<vec3f>{};
  for (auto i = 0; i < result.size() - 1; ++i) {
    auto path = path_positions(
        compute_geodesic_path(data, geometry, solver, result[i], result[i + 1]),
        data.triangles, data.positions, geometry.adjacencies);
    polyline.insert(polyline.end(), path.begin(), path.end());
  }
  return polyline;
}

vector<float> bspline_wgts(const float &t) {
  auto result = vector<float>(4);
  result[0] = 1.f / 6 * (-pow(t, 3) + 3 * pow(t, 2) - 3 * t + 1);
  result[1] = 1.f / 6 * (3 * pow(t, 3) - 6 * pow(t, 2) + 4);
  result[2] = 1.f / 6 * (-3 * pow(t, 3) + 3 * pow(t, 2) + 3 * t + 1);
  result[3] = 1.f / 6 * pow(t, 3);
  return result;
}
vector<float> rational_bspline_wgts(const float &t, const vector<float> &w) {
  auto result = bspline_wgts(t);
  auto den = 0.f;
  for (auto i = 0; i < w.size(); ++i) {
    result[i] *= w[i];
    den += result[i];
  }
  for (auto j = 0; j < w.size(); ++j) {
    result[j] /= den;
  }
  return result;
}

inline mesh_point geodesic_lerp(const shape_data &data,
                                const shape_geometry &geometry,
                                const dual_geodesic_solver &solver,
                                const mesh_point &start, const mesh_point &end,
                                float t) {

  if (start.face == end.face) {
    return mesh_point{start.face, lerp(start.uv, end.uv, t)};
  }

  auto path = compute_geodesic_path(data, geometry, solver, start, end);
  auto point = eval_path_point(path, data.triangles, data.positions,
                               geometry.adjacencies, t);

  return point;
}
mesh_point warm_start_bspline(const shape_data &data,
                              const shape_geometry &geometry,
                              const dual_geodesic_solver &solver,
                              const mesh_point &p0, const mesh_point &p1,
                              const mesh_point &p2) {
  auto q0 = geodesic_lerp(data, geometry, solver, p0, p1, 2.f / 3);
  auto q1 = geodesic_lerp(data, geometry, solver, p1, p2, 1.f / 3);
  return geodesic_lerp(data, geometry, solver, q0, q1, 0.5);
}
std::tuple<vector<vec3f>, vector<mesh_point>> trace_B_spline_segment(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const vector<vector<float>> &fields,
    const vector<vector<vec3f>> &grds, const vector<mesh_point> &control_points,
    mesh_point &seed, const int i, const int k) {

  auto curr_control_points = vector<mesh_point>(4);
  auto curr_grads = vector<vector<vec3f>>(4);
  auto curr_fields = vector<vector<float>>(4);
  for (auto j = -3; j <= 0; ++j) {
    curr_control_points[j + 3] = control_points[i + j];
    curr_grads[j + 3] = grds[i + j];
    curr_fields[j + 3] = fields[i + j];
  }
  auto step = 1.f / pow(2, k);
  auto t = 0.f;
  if (seed.face == -1)
    seed = warm_start_bspline(data, geometry, solver, curr_control_points[0],
                              curr_control_points[1], curr_control_points[2]);
  auto wgts = bspline_wgts(t);
  auto pos = zero3f;
  auto point = mesh_point{};
  auto count = 0;
  auto points = vector<mesh_point>((int)pow(2, k) + 1);
  for (auto j = 0; j <= pow(2, k); ++j) {
    std::tie(pos, point) = Riemannian_Newton_method(
        data, geometry, op, curr_grads, curr_fields, wgts, seed);
    if (point.face != -1) {
      points[j] = point;
      seed = point;
    } else if (j > 0) {
      ++count;
      points[j] = points[j - 1];
    } else
      points[0] = seed;

    t += step;
    wgts = bspline_wgts(t);
  }

  auto polyline = vector<vec3f>{};
  seed = points.back();
  for (auto i = 0; i < points.size() - 1; ++i) {
    auto path = path_positions(
        compute_geodesic_path(data, geometry, solver, points[i], points[i + 1]),
        data.triangles, data.positions, geometry.adjacencies);
    polyline.insert(polyline.end(), path.begin(), path.end());
  }

  return {polyline, points};
}
std::tuple<vector<vec3f>, vector<mesh_point>> trace_rational_bspline_segment(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const vector<vector<float>> &fields,
    const vector<vector<vec3f>> &grds, const vector<mesh_point> &control_points,
    const vector<float> &weights, mesh_point &seed, const int i, const int k) {

  auto curr_control_points = vector<mesh_point>(4);
  auto curr_grads = vector<vector<vec3f>>(4);
  auto curr_fields = vector<vector<float>>(4);
  auto curr_weights = vector<float>(4);

  for (auto j = -3; j <= 0; ++j) {
    curr_control_points[j + 3] = control_points[i + j];
    curr_grads[j + 3] = grds[i + j];
    curr_fields[j + 3] = fields[i + j];
    curr_weights[j + 3] = weights[i + j];
  }
  auto step = 1.f / pow(2, k);
  auto t = 0.f;
  if (seed.face == -1)
    seed = warm_start_bspline(data, geometry, solver, curr_control_points[0],
                              curr_control_points[1], curr_control_points[2]);
  auto wgts = rational_bspline_wgts(t, curr_weights);
  auto pos = zero3f;
  auto point = mesh_point{};
  auto count = 0;
  auto points = vector<mesh_point>((int)pow(2, k) + 1);
  for (auto j = 0; j <= pow(2, k); ++j) {
    std::tie(pos, point) = Riemannian_Newton_method(
        data, geometry, op, curr_grads, curr_fields, wgts, seed);
    if (point.face != -1) {
      points[j] = point;
      seed = point;
    } else if (j > 0) {
      ++count;
      points[j] = points[j - 1];
    } else
      points[0] = seed;

    t += step;
    wgts = rational_bspline_wgts(t, curr_weights);
  }

  auto polyline = vector<vec3f>{};
  seed = points.back();
  for (auto i = 0; i < points.size() - 1; ++i) {
    auto path = path_positions(
        compute_geodesic_path(data, geometry, solver, points[i], points[i + 1]),
        data.triangles, data.positions, geometry.adjacencies);
    polyline.insert(polyline.end(), path.begin(), path.end());
  }

  return {polyline, points};
}

std::tuple<vector<vector<vec3f>>, vector<vector<mesh_point>>>
trace_bspline(const shape_data &data, const shape_geometry &geometry,
              const shape_op &op, const dual_geodesic_solver &solver,
              const vector<vector<float>> &fields,
              const vector<vector<vec3f>> &grds,
              const vector<mesh_point> &control_points, const int k) {
  auto n = control_points.size();
  auto result_pos = vector<vector<vec3f>>(n - 3);
  auto result_point = vector<vector<mesh_point>>(n - 3);
  auto seed = mesh_point{-1, zero2f};

  for (auto i = 3; i < n; ++i) {
    std::tie(result_pos[i - 3], result_point[i - 3]) = trace_B_spline_segment(
        data, geometry, op, solver, fields, grds, control_points, seed, i, k);
  }
  return {result_pos, result_point};
}

std::tuple<vector<vector<vec3f>>, vector<vector<mesh_point>>>
trace_rational_bspline(const shape_data &data, const shape_geometry &geometry,
                       const shape_op &op, const dual_geodesic_solver &solver,
                       const vector<vector<float>> &fields,
                       const vector<vector<vec3f>> &grds,
                       const vector<mesh_point> &control_points,
                       const vector<float> &weights, const int k) {
  auto n = control_points.size();
  auto result_pos = vector<vector<vec3f>>(n - 3);
  auto result_point = vector<vector<mesh_point>>(n - 3);
  auto seed = mesh_point{-1, zero2f};

  for (auto i = 3; i < n; ++i)
    std::tie(result_pos[i - 3], result_point[i - 3]) =
        trace_rational_bspline_segment(data, geometry, op, solver, fields, grds,
                                       control_points, weights, seed, i, k);

  return {result_pos, result_point};
}

void update_rational_spline(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const vector<vector<float>> &fields,
    const vector<vector<vec3f>> &grds, const vector<mesh_point> &control_points,
    const vector<float> &weights, const int k, const int moved_point,
    vector<vector<vec3f>> &pos, vector<vector<mesh_point>> &points) {
  auto segments_to_update = vector<int>{};
  segments_to_update.reserve(4);
  auto n = control_points.size();
  for (auto i = 3; i < n; ++i) {
    if (moved_point >= i - 3 && moved_point <= i)
      segments_to_update.push_back(i);
  }
  mesh_point seed = {};
  if (segments_to_update[0] > 3)
    seed = points[segments_to_update[0] - 4].back();

  for (auto i : segments_to_update) {
    std::tie(pos[i - 3], points[i - 3]) =
        trace_rational_bspline_segment(data, geometry, op, solver, fields, grds,
                                       control_points, weights, seed, i, k);
  }
}
vector<float> rational_weights(const vector<float> &weights, const float &t) {
  auto n = weights.size();
  vector<float> w(n);
  bernstein_polynomials(n, t, w);
  auto result = vector<float>(n);
  auto den = 0.f;
  for (auto i = 0; i < n; ++i) {
    result[i] = weights[i] * w[i];
    den += result[i];
  }
  for (auto i = 0; i < n; ++i) {
    result[i] /= den;
  }

  return result;
}

vector<mesh_point> circle_control_points(const shape_data &data,
                                         const shape_geometry &geometry,
                                         const mesh_point &center,
                                         const float &r, const int n) {
  auto result = vector<mesh_point>(2 * n + 1);
  auto alpha = 2 * pif / (2 * n);
  auto R = r / std::cos(alpha);

  for (auto i = 0; i < n; ++i) {
    auto odd_v = vec2f{R * std::sin((2 * i + 1) * alpha),
                       -R * std::cos((2 * i + 1) * alpha)};
    auto even_v =
        vec2f{r * std::sin(2 * i * alpha), -r * std::cos(2 * i * alpha)};

    auto [pos, end] = polthier_straightest_geodesic(data, geometry, center,
                                                    even_v, length(even_v));
    result[2 * i] = end;
    std::tie(pos, end) = polthier_straightest_geodesic(data, geometry, center,
                                                       odd_v, length(odd_v));
    result[2 * i + 1] = end;
  }

  result[2 * n] = result[0];
  return result;
}
vector<vec3f> rational_bézier_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const geodesic_solver &geo_solver, const dual_geodesic_solver &solver,
    const vector<mesh_point> &control_points, const vector<float> &weights,
    const int k, const int method) {
  vector<vector<float>> f =
      fields_from_lndmrks(data, geometry, geo_solver, control_points, method);
  auto grds = grads_from_lndmrks(op, f);
  auto seed = control_points[0];
  double step = 1 / pow(2, k);
  auto t = step;
  vector<float> w = rational_weights(weights, t);
  auto curve = vector<mesh_point>((int)pow(2, k));
  curve[0] = seed;
  auto pos = zero3f;
  auto point = mesh_point{};
  auto count = 0;
  for (auto i = 1; i < pow(2, k); ++i) {

    std::tie(pos, point) =
        Riemannian_Newton_method(data, geometry, op, grds, f, w, seed);

    if (point.face != -1) {
      curve[i] = point;
      seed = point;
    } else {
      ++count;
      curve[i] = curve[i - 1];
    }

    t += step;
    w = rational_weights(weights, t);
  }
  curve.push_back(control_points.back());
  auto polyline = vector<vec3f>{};
  for (auto i = 0; i < curve.size() - 1; ++i) {
    auto path = path_positions(
        compute_geodesic_path(data, geometry, solver, curve[i], curve[i + 1]),
        data.triangles, data.positions, geometry.adjacencies);
    polyline.insert(polyline.end(), path.begin(), path.end());
  }
  return polyline;
}
vector<vec3f> rational_bézier_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const geodesic_solver &geo_solver,
    const vector<vector<float>> &f, const vector<vector<vec3f>> &grds,
    const vector<mesh_point> &control_points, const vector<float> &weights,
    const int k) {
  auto seed = control_points[0];
  double step = 1 / pow(2, k);
  auto t = step;
  vector<float> w = rational_weights(weights, t);
  auto curve = vector<mesh_point>((int)pow(2, k));
  curve[0] = seed;
  auto pos = zero3f;
  auto point = mesh_point{};
  auto count = 0;
  for (auto i = 1; i < pow(2, k); ++i) {

    std::tie(pos, point) =
        Riemannian_Newton_method(data, geometry, op, grds, f, w, seed);

    if (point.face != -1) {
      curve[i] = point;
      seed = point;
    } else {
      ++count;
      curve[i] = curve[i - 1];
    }

    t += step;
    w = rational_weights(weights, t);
  }
  curve.push_back(control_points.back());

  auto polyline = vector<vec3f>{};
  for (auto i = 0; i < curve.size() - 1; ++i) {
    auto path = path_positions(
        compute_geodesic_path(data, geometry, solver, curve[i], curve[i + 1]),
        data.triangles, data.positions, geometry.adjacencies);
    polyline.insert(polyline.end(), path.begin(), path.end());
  }
  return polyline;
}

vector<vector<vec3f>>
trace_circle(const shape_data &data, const shape_geometry &geometry,
             const shape_op &op, const geodesic_solver &geo_solver,
             const dual_geodesic_solver &solver, const mesh_point &center,
             const float &r, const int num_segments, const int k,
             const int type_of_solver) {
  auto control_points =
      circle_control_points(data, geometry, center, r, num_segments);

  auto f = fields_from_lndmrks(data, geometry, geo_solver, control_points,
                               type_of_solver);

  auto grds = grads_from_lndmrks(op, f);
  auto n_p = control_points.size();
  auto result = vector<vector<vec3f>>(num_segments);
  auto weights = vector<float>(3, 1);
  weights[1] = std::cos(2 * pif / (2 * num_segments));
  auto curr_points = vector<mesh_point>(3);
  auto curr_fields = vector<vector<float>>(3);
  auto curr_grads = vector<vector<vec3f>>(3);
  for (auto i = 0; i < num_segments; ++i) {

    curr_points[0] = control_points[2 * i];
    curr_points[1] = control_points[2 * i + 1];
    curr_points[2] = control_points[2 * i + 2];

    curr_fields[0] = f[2 * i];
    curr_fields[1] = f[2 * i + 1];
    curr_fields[2] = f[2 * i + 2];

    curr_grads[0] = grds[2 * i];
    curr_grads[1] = grds[2 * i + 1];
    curr_grads[2] = grds[2 * i + 2];
    result[i] =
        rational_bézier_curve(data, geometry, op, solver, geo_solver,
                              curr_fields, curr_grads, curr_points, weights, k);
  }

  return result;
}
vector<float> béz_interp_weights(const vector<float> &betas,
                                 const vector<float> &ti, const float &t) {
  auto n = betas.size();
  auto w = vector<float>(n);
  auto sum = 0.f;
  for (auto i = 0; i < n; ++i) {
    if (ti[i] == t) {
      w = vector<float>(n, 0.f);
      w[i] = 1;
      return w;
    } else {
      w[i] = pow(-1, i) * betas[i] / (t - ti[i]);
      sum += w[i];
    }
  }
  for (auto i = 0; i < n; ++i) {
    w[i] /= sum;
  }

  return w;
}
vector<vec3f> béz_interp_rational_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const geodesic_solver &geo_solver, const dual_geodesic_solver &solver,
    const vector<vector<float>> &f, const vector<vector<vec3f>> &grds,
    const vector<mesh_point> &control_points, const vector<float> &betas,
    const int k) {

  auto tis = vector<float>{0.f, 1.f / 2, 1};
  auto seed = control_points[0];
  double step = 1 / pow(2, k);
  auto t = step;
  vector<float> w = béz_interp_weights(betas, tis, t);
  auto curve = vector<mesh_point>((int)pow(2, k));
  curve[0] = seed;
  auto pos = zero3f;
  auto point = mesh_point{};
  auto count = 0;
  for (auto i = 1; i < pow(2, k); ++i) {

    std::tie(pos, point) =
        Riemannian_Newton_method(data, geometry, op, grds, f, w, seed);

    if (point.face != -1) {
      curve[i] = point;
      seed = point;
    } else {
      ++count;
      curve[i] = curve[i - 1];
    }

    t += step;
    w = béz_interp_weights(betas, tis, t);
  }
  curve.push_back(control_points.back());

  return polyline_pos(data.triangles, data.positions, curve);
}
vector<vec3f> PCE_grad(const Eigen::SparseMatrix<double> &G,
                       const Eigen::VectorXd &f, const int F) {
  auto result = vector<vec3f>(F);
  Eigen::VectorXd grad = G * f;
  for (auto i = 0; i < F; ++i) {
    result[i].x = -grad(3 * i);
    result[i].y = -grad(3 * i + 1);
    result[i].z = -grad(3 * i + 2);
    // result[i] = normalize(result[i]);
  }

  return result;
}
vector<vec3f> AGS_grad(const Eigen::SparseMatrix<double> &G,
                       const Eigen::VectorXd &f) {
  int V = (int)f.rows();
  auto result = vector<vec3f>(V);
  Eigen::VectorXd grad = G * f;
  for (auto i = 0; i < V; ++i) {
    result[i].x = -grad(3 * i);
    result[i].y = -grad(3 * i + 1);
    result[i].z = -grad(3 * i + 2);
  }

  return result;
}
