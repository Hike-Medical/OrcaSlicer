// ZAA: Z Anti-Aliasing / Non-planar contouring
// Based on BambuStudio-ZAA by adob (https://github.com/adob/BambuStudio-ZAA)
// Adapted for OrcaSlicer with z_offsets approach.
//
// For each extrusion point on top surfaces, external perimeters, and ironing:
//   - Subdivide the path into small segments (~0.1mm)
//   - Cast ray up/down from each point to find the mesh surface distance
//   - Store the Z offset so the nozzle follows the actual surface

#include "Exception.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Layer.hpp"
#include "Print.hpp"
#include "Point.hpp"
#include "libslic3r.h"
#include <cfloat>
#include <cmath>
#include <string>

namespace Slic3r {

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

    const Points &points = path.polyline.points;
    if (points.size() < 2)
        return false;

    double resolution_mm = 0.1;
    coordf_t height = layer->height;
    double minimize_perimeter_height_angle = region->region().config().zaa_minimize_perimeter_height;

    // Build new polyline with subdivided points and corresponding z_offsets
    Points new_points;
    std::vector<coordf_t> new_z_offsets;
    bool was_contoured = false;

    for (Points::const_iterator it = points.begin(); it != points.end() - 1; ++it) {
        Vec2d p1d(unscale_(it->x()), unscale_(it->y()));
        Vec2d p2d(unscale_((it + 1)->x()), unscale_((it + 1)->y()));
        Vec2d delta = p2d - p1d;

        double length_mm = delta.norm();
        int num_segments = int(std::ceil(length_mm / resolution_mm));

        for (int i = 0; i < num_segments + 1; i++) {
            Vec2d p = p1d + delta * i / num_segments;

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

            // Point simplification: if the new point is collinear with the previous two, replace the last
            if (new_points.size() > 1) {
                Vec2d prev2(unscale_(new_points[new_points.size() - 2].x()), unscale_(new_points[new_points.size() - 2].y()));
                Vec2d prev1(unscale_(new_points[new_points.size() - 1].x()), unscale_(new_points[new_points.size() - 1].y()));
                Vec2d curr(p.x(), p.y());
                // Check if prev1 is on the line from prev2 to curr and z_offsets are similar
                Vec2d v1 = prev1 - prev2;
                Vec2d v2 = curr - prev2;
                double cross = v1.x() * v2.y() - v1.y() * v2.x();
                double z_diff = std::abs(new_z_offsets.back() - d);
                if (std::abs(cross) < EPSILON * v2.norm() && z_diff < EPSILON) {
                    new_points.back() = Point(scale_(p.x()), scale_(p.y()));
                    new_z_offsets.back() = d;
                    continue;
                }
            }

            new_points.emplace_back(scale_(p.x()), scale_(p.y()));
            new_z_offsets.push_back(d);
        }
    }

    if (!was_contoured) {
        return false;
    }

    path.polyline = Polyline(std::move(new_points));
    path.polyline.fitting_result.clear();
    path.z_offsets = std::move(new_z_offsets);
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
