#ifndef UTILITIES_H
#define UTILITIES_H

#include <diff_geo/diff_geo.h>
#include <fstream>
#include <utils/mesh_io.h>
#include <yocto/yocto_sceneio.h>

vector<vec3f> PCE_grad(const Eigen::SparseMatrix<double> &G,
                       const Eigen::VectorXd &f, const int F);
vector<vec3f> AGS_grad(const Eigen::SparseMatrix<double> &G,
                       const Eigen::VectorXd &f);

vector<vector<float>> fields_from_lndmrks(const shape_data &data,
                                          const shape_geometry &geometry,
                                          const geodesic_solver &geo_solver,
                                          const vector<mesh_point> &lndmrks,
                                          const int solver);
void update_fields_from_lndmrks(const shape_data &data,
                                const shape_geometry &geometry,
                                const geodesic_solver &geo_solver,
                                const vector<mesh_point> &lndmrks,
                                const int solver, const vector<bool> &moved,
                                vector<vector<float>> &f);

vector<vector<vec3f>> grads_from_lndmrks(const shape_op &op,
                                         const vector<vector<float>> &fields);

void update_grads_from_lndmrks(const shape_op &op,
                               const vector<vector<float>> &f,
                               const vector<bool> &moved,
                               vector<vector<vec3f>> &grds);

vector<vec3f> bézier_curve(const shape_data &data,
                           const shape_geometry &geometry, const shape_op &op,
                           const geodesic_solver &geo_solver,
                           const dual_geodesic_solver &solver,
                           const vector<mesh_point> &control_points,
                           const int k, const int type_of_solver);

vector<mesh_point> circle_control_points(const shape_data &data,
                                         const shape_geometry &geometry,
                                         const mesh_point &center,
                                         const float &radius, const int n);
vector<vector<vec3f>>
trace_circle(const shape_data &data, const shape_geometry &geometry,
             const shape_op &op, const geodesic_solver &geo_solver,
             const dual_geodesic_solver &solver, const mesh_point &center,
             const float &r, const int num_segments, const int k,
             const int type_of_solver);
vector<vec3f> rational_bézier_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const geodesic_solver &geo_solver, const dual_geodesic_solver &solver,
    const vector<mesh_point> &control_points, const vector<float> &weights,
    const int k, const int method);
vector<vec3f> rational_bézier_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const geodesic_solver &geo_solver,
    const vector<vector<float>> &f, const vector<vector<vec3f>> &grds,
    const vector<mesh_point> &control_points, const vector<float> &weights,
    const int k);
vector<vec3f> béz_interp_rational_curve(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const geodesic_solver &geo_solver, const dual_geodesic_solver &solver,
    const vector<vector<float>> &f, const vector<vector<vec3f>> &grds,
    const vector<mesh_point> &control_points, const vector<float> &betas,
    const int k);

std::tuple<vector<vector<vec3f>>, vector<vector<mesh_point>>>
trace_bspline(const shape_data &data, const shape_geometry &geometry,
              const shape_op &op, const dual_geodesic_solver &solver,
              const vector<vector<float>> &fields,
              const vector<vector<vec3f>> &grds,
              const vector<mesh_point> &control_points, const int k);

std::tuple<vector<vector<vec3f>>, vector<vector<mesh_point>>>
trace_rational_bspline(const shape_data &data, const shape_geometry &geometry,
                       const shape_op &op, const dual_geodesic_solver &solver,
                       const vector<vector<float>> &fields,
                       const vector<vector<vec3f>> &grds,
                       const vector<mesh_point> &control_points,
                       const vector<float> &weights, const int k);

void update_rational_spline(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const dual_geodesic_solver &solver, const vector<vector<float>> &fields,
    const vector<vector<vec3f>> &grds, const vector<mesh_point> &control_points,
    const vector<float> &weights, const int k, const int moved_point,
    vector<vector<vec3f>> &pos, vector<vector<mesh_point>> &points);

#endif