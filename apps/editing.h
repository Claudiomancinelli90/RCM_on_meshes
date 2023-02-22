#pragma once

using namespace yocto;

inline bool is_anchor_control_point(int i) {
  if (i == -1)
    return false;
  return i % 3 == 0;
}
inline bool is_handle_control_point(int i) {
  if (i == -1)
    return false;
  return !is_anchor_control_point(i);
}

// inline mat2f parallel_transport_rotation(const bezier_mesh &mesh,
//                                          const mesh_point &start,
//                                          const mesh_point &end) {
//   auto path = compute_geodesic_path(mesh, start, end);
//   return parallel_transport_rotation(mesh.triangles, mesh.positions,
//                                      mesh.adjacencies, path);
// }
