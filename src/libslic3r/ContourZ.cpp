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

	// Debug: log ALL external perimeter paths
	int cnt = s_zaa_debug_counter.fetch_add(1);

	coordf_t height = layer->height;
	double minimize_perimeter_height_angle = region->region().config().zaa_minimize_perimeter_height;

	Pointf3s contoured_points;
	bool was_contoured = false;

	// If lslices_original is populated (make_overhang_printable active),
	// use it to filter out points in overhang-expanded areas.
	const ExPolygons &orig_slices = layer->lslices_original;
	bool has_orig_slices = !orig_slices.empty();

	// Helper lambda: compute clamped Z offset d for a given XY point
	auto compute_d = [&](double x, double y) -> double {
		// Skip contouring for points outside original mesh footprint
		// (overhang-expanded areas where raycast against original mesh is unreliable)
		if (has_orig_slices) {
			Point pt(scale_(x), scale_(y));
			bool inside = false;
			for (const ExPolygon &ep : orig_slices) {
				if (ep.contains(pt)) {
					inside = true;
					break;
				}
			}
			if (!inside)
				return 0.0;
		}

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

		// For minimize_perimeter_height, use DOWN hit slope (side wall)
		// instead of closest hit slope. At the mesh edge, the UP hit
		// sees the gentle top surface while the DOWN hit sees the steep
		// side wall that determines the actual perimeter height.
		double down_slope_rad = slope_from_normal(hit_down.normal());
		double down_slope_degrees = down_slope_rad * 180.0 / M_PI;

		if (d > min_down_val && minimize_perimeter_height_angle > 0 && minimize_perimeter_height_angle < down_slope_degrees && path.role() == erExternalPerimeter) {
			double adjustment = follow_slope_down(down_slope_rad, half_width);
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

		// BambuStudio-ZAA zeros d>0 for external perimeters to avoid seam appearance.
		// However, this relies on ground_level being large negative (from 3MF project
		// transforms), making d naturally negative for walls. In OrcaSlicer CLI mode,
		// objects are placed on bed (ground_level≈0), so d is near-zero or slightly
		// positive for walls — zeroing kills all wall contouring.
		// Fix: only zero d that exceeds max_up (already clamped above). Small positive
		// d values (0 to max_up=0.05mm) are physically acceptable and needed for
		// wall contouring when ground_level≈0.
		// The clamping block above already handles d > max_up + 0.03 → d = 0.

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

	// Debug: log every external perimeter path result
	if (path.role() == erExternalPerimeter && cnt < 200) {
		// Find max |d| in contoured_points
		double max_abs_d = 0;
		int nonzero_count = 0;
		for (const auto &cp : contoured_points) {
			if (std::abs(cp.z()) > EPSILON) nonzero_count++;
			if (std::abs(cp.z()) > max_abs_d) max_abs_d = std::abs(cp.z());
		}
		FILE *df = fopen("/tmp/zaa_contour_debug.txt", "a");
		if (df) {
			fprintf(df, "path_%d layer=%zu pts_in=%zu pts_out=%zu was_contoured=%d nonzero_d=%d max_d=%.6f\n",
				cnt, layer->id(), points.size(), contoured_points.size(),
				was_contoured ? 1 : 0, nonzero_count, max_abs_d);
			// Log first few points with d values
			for (size_t i = 0; i < std::min(contoured_points.size(), (size_t)5); i++) {
				fprintf(df, "  [%zu] x=%.3f y=%.3f d=%.6f\n", i,
					contoured_points[i].x(), contoured_points[i].y(), contoured_points[i].z());
			}
			fclose(df);
		}
	}

	if (!was_contoured) {
		return false;
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
	static std::atomic<int> s_filter_debug{0};
	for (ExtrusionEntity *extr : collection.entities) {
		if (!contains(roles, extr->role())) {
			int fc = s_filter_debug.fetch_add(1);
			if (fc < 50) {
				FILE *df = fopen("/tmp/zaa_filter_debug.txt", "a");
				if (df) {
					fprintf(df, "FILTERED layer=%zu role=%s type=%s\n",
						region->layer()->id(),
						ExtrusionEntity::role_to_string(extr->role()).c_str(),
						dynamic_cast<ExtrusionLoop*>(extr) ? "Loop" :
						dynamic_cast<ExtrusionEntityCollection*>(extr) ? "Collection" :
						dynamic_cast<ExtrusionPath*>(extr) ? "Path" : "Other");
					fclose(df);
				}
			}
			continue;
		}
		contour_extrusion_entity(region, mesh, extr);
	}
}

void Layer::make_contour_z(const sla::IndexedMesh &mesh)
{
	for (LayerRegion *region : this->regions()) {
		handle_extrusion_collection(region, mesh, region->fills, {erTopSolidInfill, erIroning, erExternalPerimeter, erMixed});
		// Process ALL perimeter entities — OrcaSlicer wraps them in collections
		// with various top-level roles. The path-level filter in contour_extrusion_path
		// already checks for erExternalPerimeter/erPerimeter.
		for (ExtrusionEntity *extr : region->perimeters.entities) {
			contour_extrusion_entity(region, mesh, extr);
		}
	}
}

} // namespace Slic3r
