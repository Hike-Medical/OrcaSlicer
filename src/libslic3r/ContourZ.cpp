// ZAA: Z Anti-Aliasing / Non-planar contouring
// Based on BambuStudio-ZAA by adob (https://github.com/adob/BambuStudio-ZAA)
// Ported to OrcaSlicer with Polyline3/Point3 approach.
//
// For each extrusion point on top surfaces, external perimeters, and ironing:
//   - Subdivide the path into small segments (~0.1mm)
//   - Cast ray up/down from each point to find the mesh surface distance
//   - Store the Z offset in the 3D polyline point

#include "Exception.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Layer.hpp"
#include "Print.hpp"
#include "Point.hpp"
#include "Line.hpp"
#include "libslic3r.h"
#include <cfloat>
#include <cmath>
#include <initializer_list>
#include <string>
#include <atomic>
#include <boost/log/trivial.hpp>

namespace Slic3r {

static std::atomic<int> s_zaa_debug_counter{0};

static void contour_extrusion_entity(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionEntity *extr);

static double follow_slope_down(double angle_rad, double dist) {
	return -dist * std::sin(angle_rad);
}

static double slope_from_normal(const Eigen::Vector3d& normal) {
    Eigen::Vector3d n = normal.normalized();
    double angle_rad = std::acos(std::abs(n.z()));
	return angle_rad;
}

static bool contour_extrusion_path(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionPath &path) {
	if (region->region().config().zaa_region_disable) {
		return false;
	}

	if (path.role() != erTopSolidInfill && path.role() != erIroning && path.role() != erExternalPerimeter && path.role() != erPerimeter) {
		return false;
	}

	Layer *layer = region->layer();
	coordf_t mesh_z = layer->print_z + mesh.ground_level();
	coordf_t min_z = layer->object()->config().zaa_min_z;

	const Points3 &points = path.polyline.points;
	double resolution_mm = 0.1;

	// Debug: log first few outer wall paths to /tmp/zaa_contour_debug.txt
	int cnt = s_zaa_debug_counter.fetch_add(1);
	if (cnt < 20 && path.role() == erExternalPerimeter) {
		FILE *df = fopen("/tmp/zaa_contour_debug.txt", "a");
		if (df) {
			fprintf(df, "--- path %d: role=ExternalPerimeter layer_id=%zu layer_z=%.4f mesh_z=%.4f height=%.4f points=%zu\n",
				cnt, layer->id(), layer->print_z, mesh_z, layer->height, points.size());
			if (!points.empty()) {
				Vec2d p0(unscale_(points.front().x()), unscale_(points.front().y()));
				sla::IndexedMesh::hit_result h_up = mesh.query_ray_hit({p0.x(), p0.y(), mesh_z}, {0.0, 0.0, 1.0});
				sla::IndexedMesh::hit_result h_dn = mesh.query_ray_hit({p0.x(), p0.y(), mesh_z}, {0.0, 0.0, -1.0});
				double up = h_up.distance(), dn = h_dn.distance();
				double d = up < dn ? up : -dn;
				fprintf(df, "  first_pt=(%.4f, %.4f) up=%.4f down=%.4f d=%.4f\n", p0.x(), p0.y(), up, dn, d);
			}
			fclose(df);
		}
	}

	coordf_t height = layer->height;
	double minimize_perimeter_height_angle = region->region().config().zaa_minimize_perimeter_height;

	Pointf3s contoured_points;
	bool was_contoured = false;

	// Helper lambda: compute clamped Z offset d for a given XY point
	auto compute_d = [&](double x, double y) -> double {
		sla::IndexedMesh::hit_result hit_up = mesh.query_ray_hit({x, y, mesh_z}, {0.0, 0.0, 1.0});
		sla::IndexedMesh::hit_result hit_down = mesh.query_ray_hit({x, y, mesh_z}, {0.0, 0.0, -1.0});

		double up = hit_up.distance();
		double down = hit_down.distance();
		double d = up < down ? up : -down;
		const Vec3d &normal = (up < down ? hit_up : hit_down).normal();

		double max_up = min_z;
		double min_down_val = -(height - min_z);
		double half_width = path.width / 2.0;
		if (path.role() == erIroning) {
			max_up = height;
			min_down_val = -(height + 0.1);
		}

		double slope_rad = slope_from_normal(normal);
		double slope_degrees = slope_rad * 180.0 / M_PI;

		if (d > min_down_val && minimize_perimeter_height_angle > 0 && minimize_perimeter_height_angle < slope_degrees && path.role() == erExternalPerimeter) {
			double adjustment = follow_slope_down(slope_rad, half_width);
			if (adjustment > 0) {
				throw RuntimeError("ContourZ: got positive adjustment");
			}
			d += adjustment;
			if (d < min_down_val) {
				d = min_down_val;
			}
		}

		if (d > max_up + 0.03 || d < min_down_val) {
			d = 0;
		} else {
			if (d > max_up) {
				d = max_up;
			}
		}

		if (path.role() == erExternalPerimeter && d > 0) {
			d = 0;
		}

		// Skip contouring for points outside original mesh boundary
		// (e.g. overhang-expanded areas from make_overhang_printable).
		// If no mesh above (up=inf) and closest surface below is steep
		// (normal mostly horizontal), the point is on an expanded overhang
		// shelf, not a genuine top surface — force d=0.
		if (std::abs(d) > EPSILON && std::isinf(up) && !std::isinf(down) && slope_degrees > 50.0) {
			d = 0;
		}

		return d;
	};

	// Add the first original point
	{
		Vec2d p0(unscale_(points.front().x()), unscale_(points.front().y()));
		double d0 = compute_d(p0.x(), p0.y());
		if (std::abs(d0) > EPSILON) was_contoured = true;
		contoured_points.push_back({p0.x(), p0.y(), d0});
	}

	for (Points3::const_iterator it = points.begin(); it != points.end()-1; ++it) {
		Vec2d p1d(unscale_(it->x()), unscale_(it->y()));
		Vec2d p2d(unscale_((it+1)->x()), unscale_((it+1)->y()));

		Vec2d delta = p2d - p1d;
		double length_mm = delta.norm();
		int num_segments = int(std::ceil(length_mm / resolution_mm));

		// Compute d for all subdivided points in this segment (single pass)
		std::vector<std::pair<Vec2d, double>> sub_points;
		bool segment_has_contour = false;
		for (int i = 1; i <= num_segments; i++) {
			Vec2d p = p1d + delta * i / num_segments;
			double d = compute_d(p.x(), p.y());
			sub_points.push_back({p, d});
			if (std::abs(d) > EPSILON)
				segment_has_contour = true;
		}

		if (!segment_has_contour) {
			// All d=0 in this segment: emit only the original endpoint.
			// This preserves the original curve geometry from the slicer.
			contoured_points.push_back({p2d.x(), p2d.y(), 0.0});
		} else {
			// This segment has Z variation: emit all subdivided points
			was_contoured = true;
			for (auto &[p, d] : sub_points) {
				Vec3d new_point = {p.x(), p.y(), d};

				// Point simplification: merge collinear subdivided points
				if (contoured_points.size() > 2) {
					double dist = line_alg::distance_to_infinite_squared(
						Linef3(contoured_points[contoured_points.size() - 2],
						        contoured_points[contoured_points.size() - 1]),
						new_point);
					if (dist < EPSILON) {
						contoured_points[contoured_points.size() - 1] = new_point;
						continue;
					}
				}

				contoured_points.push_back(new_point);
			}
		}
	}

	if (!was_contoured) {
		return false;
	}

	// Debug: log when outer wall gets contoured (this causes jagged walls)
	if (path.role() == erExternalPerimeter && cnt < 20) {
		FILE *df = fopen("/tmp/zaa_contour_debug.txt", "a");
		if (df) {
			// Find the point with max |d| to understand what triggered was_contoured
			double max_abs_d = 0;
			size_t max_idx = 0;
			for (size_t i = 0; i < contoured_points.size(); i++) {
				if (std::abs(contoured_points[i].z()) > max_abs_d) {
					max_abs_d = std::abs(contoured_points[i].z());
					max_idx = i;
				}
			}
			fprintf(df, "  WAS_CONTOURED! points=%zu max_d=%.6f at idx=%zu (%.4f, %.4f)\n",
				contoured_points.size(), contoured_points[max_idx].z(), max_idx,
				contoured_points[max_idx].x(), contoured_points[max_idx].y());
			fclose(df);
		}
	}

	Polyline3 polyline;
	for (const Vec3d &point : contoured_points) {
		polyline.append(Point3(scale_(point.x()), scale_(point.y()), scale_(point.z())));
	}

	path.polyline = std::move(polyline);
	path.z_contoured = true;
	return true;
}

static void contour_extrusion_multipath(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionMultiPath &multipath)
{
	for (ExtrusionPath &path : multipath.paths) {
		contour_extrusion_path(region, mesh, path);
	}
}

static void contour_extrusion_loop(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionLoop &loop)
{
	for (ExtrusionPath &path : loop.paths) {
		contour_extrusion_path(region, mesh, path);
	}
}

static void contour_extrusion_entity_collection(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionEntityCollection &collection) {
	for (ExtrusionEntity *entity : collection.entities) {
		contour_extrusion_entity(region, mesh, entity);
	}
}

static void contour_extrusion_entity(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionEntity *extr) {
	ExtrusionMultiPath *multipath = dynamic_cast<ExtrusionMultiPath*>(extr);
	if (multipath != nullptr) {
		contour_extrusion_multipath(region, mesh, *multipath);
		return;
	}

	ExtrusionPath *path = dynamic_cast<ExtrusionPath*>(extr);
	if (path != nullptr) {
		contour_extrusion_path(region, mesh, *path);
		return;
	}

	ExtrusionLoop *loop = dynamic_cast<ExtrusionLoop*>(extr);
	if (loop != nullptr) {
		contour_extrusion_loop(region, mesh, *loop);
		return;
	}

	ExtrusionEntityCollection *collection = dynamic_cast<ExtrusionEntityCollection*>(extr);
	if (collection != nullptr) {
		contour_extrusion_entity_collection(region, mesh, *collection);
		return;
	}

	// Other types (ExtrusionPathSloped, ExtrusionLoopSloped) — skip silently
}

static void handle_extrusion_collection(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionEntityCollection &collection, std::initializer_list<ExtrusionRole> roles) {
	for (ExtrusionEntity *extr : collection.entities) {
		if (!contains(roles, extr->role())) {
			continue;
		}
		contour_extrusion_entity(region, mesh, extr);
	}
}

void Layer::make_contour_z(const sla::IndexedMesh &mesh)
{
	for (LayerRegion *region : this->regions()) {
		handle_extrusion_collection(region, mesh, region->fills, {erTopSolidInfill, erIroning, erExternalPerimeter, erMixed});
		handle_extrusion_collection(region, mesh, region->perimeters, {erExternalPerimeter, erMixed});
	}
}

} // namespace Slic3r
