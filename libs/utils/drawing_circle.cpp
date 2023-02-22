#include "drawing_circle.h"
using namespace yocto;
vector<vec3f> closed_curve_positions(const closed_curve &curve,
                                     const vector<vec3i> &triangles,
                                     const vector<vec3f> &positions,
                                     const vector<vec3i> &adjacencies,
                                     const vec2i &range) {

  auto pos = vector<vec3f>(curve.lerps.size() + 1);
  auto s = curve.strip.size();
  for (auto i = 0; i < s; i++) {
    auto e = get_edge(triangles, positions, adjacencies, curve.strip[i],
                      curve.strip[(i + 1) % s]);
    if (e.x == -1)
      continue;
    auto x = curve.lerps[i];
    auto p = lerp(positions[e.x], positions[e.y], x);
    pos[i] = p;
  }
  pos.back() = pos.front();
  auto result = vector<vec3f>{};
  if (range.x < range.y) {
    result.insert(result.begin(), pos.begin() + range.x,
                  pos.begin() + range.y + 1);
  } else if (range.x > range.y) {
    result.insert(result.begin(), pos.begin() + range.y, pos.end());
    result.insert(result.end(), pos.begin(), pos.begin() + range.x);
  } else
    result = pos;

  return result;
}

vector<vec3f> circle_positions(const vector<vec3i> &triangles,
                               const vector<vec3f> &positions,
                               const vector<vec3i> &adjacencies,
                               const Circle &c0) {
  auto pos = vector<vec3f>{};
  for (auto curve : c0.isoline) {
    auto curve_pos = closed_curve_positions(curve, triangles, positions,
                                            adjacencies, vec2i{-1, -1});
    pos.insert(pos.end(), curve_pos.begin(), curve_pos.end());
  }

  return pos;
}

std::tuple<int, float, int>
find_seed_in_circle(const vector<vec3i> &triangles,
                    const vector<vec3i> &adjacencies,
                    const vector<float> &distances, const float &radius,
                    const vector<bool> &parsed) {

  for (auto i = 0; i < triangles.size(); ++i) {
    if (parsed[i])
      continue;
    auto &t = triangles[i];
    auto d = vec3f{distances[t.x] - radius, distances[t.y] - radius,
                   distances[t.z] - radius};
    if (d.x * d.y <= 0 && !parsed[adjacencies[i][0]]) {
      return {i, -d.x / (d.y - d.x), 0};
    } else if (d.x * d.z <= 0 && !parsed[adjacencies[i][2]]) {
      return {i, -d.x / (d.z - d.x), 2};
    } else if (d.z * d.y <= 0 && !parsed[adjacencies[i][1]]) {
      return {i, -d.z / (d.y - d.z), 1};
    }
  }
  return {-1, 0, -1};
}
Isoline create_isoline(const vector<vec3i> &triangles,
                       const vector<vec3i> &adjacencies,
                       const vector<float> &distances, const float &radius) {
  Isoline iso = {};
  auto parsed = vector<bool>(triangles.size(), false);
  auto [seed, lerp, offset] =
      find_seed_in_circle(triangles, adjacencies, distances, radius, parsed);
  while (seed != -1) {
    auto curve = closed_curve{};
    while (!parsed[seed]) {

      curve.strip.push_back(seed);
      curve.lerps.push_back(lerp);
      parsed[seed] = true;

      auto next = adjacencies[seed][offset];
      auto t = triangles[next];
      auto dist = vec3f{distances[t.x] - radius, distances[t.y] - radius,
                        distances[t.z] - radius};
      auto h = find_in_vec(adjacencies[next], seed);
      if (dist[h] * dist[(h + 2) % 3] <= 0) {
        seed = next;
        offset = (h + 2) % 3;
        lerp = -dist[(h + 2) % 3] / (dist[h] - dist[(h + 2) % 3]);
      } else if (dist[(h + 1) % 3] * dist[(h + 2) % 3] <= 0) {
        seed = next;
        offset = (h + 1) % 3;
        lerp = -dist[(h + 1) % 3] / (dist[(h + 2) % 3] - dist[(h + 1) % 3]);
      } else if (dist[h] * dist[(h + 1) % 3] <= 0) { // closing circle
        seed = next;
        offset = h;
        lerp = -dist[h] / (dist[(h + 1) % 3] - dist[h]);
        break;
      } else {
        break;
      }
      lerp = yocto::clamp(lerp, 0.f, 1.f);
    }

    if (curve.strip.size() > 0) {
      iso.push_back(curve);
    }
    std::tie(seed, lerp, offset) =
        find_seed_in_circle(triangles, adjacencies, distances, radius, parsed);
  }
  return iso;
}

bool set_radius(const vector<vec3i> &triangles,
                const vector<vec3i> &adjacencies, Circle *circle,
                const float &radius) {

  if (radius == 0.f)
    return true;

  circle->isoline =
      create_isoline(triangles, adjacencies, circle->distances, radius);

  if (circle->isoline.size() == 0) {
    return false;
  }

  circle->radius = radius;

  return true;
}

Circle create_circle(const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies, const mesh_point &center,
                     const float &radius, const vector<float> &distances) {
  auto circle = Circle();
  circle.center = center;
  circle.distances = distances;

  if (!set_radius(triangles, adjacencies, &circle, radius)) {
    return Circle();
  }

  return circle;
}

// void map_tids(const bezier_mesh &mesh, Circle &circle) {

//   auto initialized = false;
//   for (auto &curve : circle.isoline) {
//     for (auto &tid : curve.strip) {
//       auto path =
//           compute_geodesic_path(mesh, circle.center, {tid, vec2f{0.33,
//           0.33}});
//       auto v = tangent_path_direction(mesh, path);
//       if (!initialized) {
//         circle.e = v;
//         circle.tids_mapping[tid] = 0.f;
//         initialized = true;

//       } else {
//         auto teta = angle(circle.e, v);
//         if (cross(circle.e, v) < 0)
//           teta = 2 * pif - teta;

//         circle.tids_mapping[tid] = teta;
//       }
//     }
//   }
// }
