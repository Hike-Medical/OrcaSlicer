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

	// Debug: log first few paths
	int cnt = s_zaa_debug_counter.fetch_add(1);
	if (cnt < 10) {
		BOOST_LOG_TRIVIAL(warning) << "ZAA contour_extrusion_path: role=" << path.role()
			<< " layer_z=" << layer->print_z
			<< " mesh_z=" << mesh_z
			<< " ground_level=" << mesh.ground_level()
			<< " height=" << layer->height
			<< " min_z_cfg=" << min_z
			<< " points=" << points.size()
			<< " width=" << path.width;
		if (!points.empty()) {
			Vec2d p0(unscale_(points.front().x()), unscale_(points.front().y()));
			sla::IndexedMesh::hit_result h_up = mesh.query_ray_hit({p0.x(), p0.y(), mesh_z}, {0.0, 0.0, 1.0});
			sla::IndexedMesh::hit_result h_dn = mesh.query_ray_hit({p0.x(), p0.y(), mesh_z}, {0.0, 0.0, -1.0});
			BOOST_LOG_TRIVIAL(warning) << "ZAA raycast at (" << p0.x() << "," << p0.y() << "," << mesh_z
				<< "): up=" << h_up.distance() << " down=" << h_dn.distance()
				<< " is_inside=" << h_up.is_inside();
		}
	}

	coordf_t height = layer->height;
	double minimize_perimeter_height_angle = region->region().config().zaa_minimize_perimeter_height;

	Pointf3s contoured_points;
	bool was_contoured = false;

	for (Points3::const_iterator it = points.begin(); it != points.end()-1; ++it) {
		Vec2d p1d(unscale_(it->x()), unscale_(it->y()));
		Vec2d p2d(unscale_((it+1)->x()), unscale_((it+1)->y()));
		Linef line(p1d, p2d);

		Vec2d delta = line.b - line.a;
		double length_mm = delta.norm();
		int num_segments = int(std::ceil(length_mm / resolution_mm));

		for (int i = 0; i < num_segments+1; i++) {
			Vec2d p = p1d + delta*i/num_segments;

			coordf_t x = p.x();
			coordf_t y = p.y();

			sla::IndexedMesh::hit_result hit_up = mesh.query_ray_hit({x, y, mesh_z}, {0.0, 0.0, 1.0});
			sla::IndexedMesh::hit_result hit_down = mesh.query_ray_hit({x, y, mesh_z}, {0.0, 0.0, -1.0});

			double up = hit_up.distance();
			double down = hit_down.distance();
			double d = up < down ? up : -down;
			const Vec3d &normal = (up < down ? hit_up : hit_down).normal();

			double max_up = min_z;
			double min_down = -(height - min_z);
			double half_width = path.width / 2.0;
			if (path.role() == erIroning) {
				max_up = height;
				min_down = -(height + 0.1);
			}

			double slope_rad = slope_from_normal(normal);
			double slope_degrees = slope_rad * 180.0 / M_PI;

			if (d > min_down && minimize_perimeter_height_angle > 0 && minimize_perimeter_height_angle < slope_degrees && path.role() == erExternalPerimeter) {
				double adjustment = follow_slope_down(slope_rad, half_width);
				if (adjustment > 0) {
					throw RuntimeError("ContourZ: got positive adjustment");
				}
				d += adjustment;
				if (d < min_down) {
					d = min_down;
				}
			}

			if (d > max_up + 0.03 || d < min_down) {
				d = 0;
			} else {
				if (d > max_up) {
					d = max_up;
				}
			}

			if (path.role() == erExternalPerimeter && d > 0) {
				// do not increase height of external perimeters as this may create an appearance of a seam
				d = 0;
			}

			if (std::abs(d) > EPSILON) {
				was_contoured = true;
			}

			Vec3d new_point = {p.x(), p.y(), d};

			// Point simplification: if the new point is collinear with the previous two, replace the last
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

static std::atomic<int> s_zaa_collection_debug{0};

static void handle_extrusion_collection(LayerRegion *region, const sla::IndexedMesh &mesh, ExtrusionEntityCollection &collection, std::initializer_list<ExtrusionRole> roles) {
	int dbg = s_zaa_collection_debug.fetch_add(1);
	for (ExtrusionEntity *extr : collection.entities) {
		if (dbg < 5) {
			BOOST_LOG_TRIVIAL(warning) << "ZAA handle_collection: entity role=" << extr->role()
				<< " type=" << typeid(*extr).name()
				<< " accepted=" << (contains(roles, extr->role()) ? "YES" : "NO");
		}
		if (!contains(roles, extr->role())) {
			continue;
		}
		contour_extrusion_entity(region, mesh, extr);
	}
}

static std::atomic<int> s_zaa_layer_debug{0};

void Layer::make_contour_z(const sla::IndexedMesh &mesh)
{
	int dbg = s_zaa_layer_debug.fetch_add(1);
	if (dbg < 3) {
		BOOST_LOG_TRIVIAL(warning) << "ZAA make_contour_z: layer print_z=" << this->print_z
			<< " height=" << this->height
			<< " regions=" << this->regions().size()
			<< " ground_level=" << mesh.ground_level();
		for (size_t i = 0; i < this->regions().size(); i++) {
			LayerRegion *r = this->regions()[i];
			BOOST_LOG_TRIVIAL(warning) << "ZAA   region[" << i << "]: fills=" << r->fills.entities.size()
				<< " perimeters=" << r->perimeters.entities.size();
		}
	}
	for (LayerRegion *region : this->regions()) {
		handle_extrusion_collection(region, mesh, region->fills, {erTopSolidInfill, erIroning, erExternalPerimeter, erMixed});
		handle_extrusion_collection(region, mesh, region->perimeters, {erExternalPerimeter, erMixed});
	}
}

} // namespace Slic3r
