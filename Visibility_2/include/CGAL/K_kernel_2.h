// Copyright (c) 2026 Torotono Metropolitan University (Canada).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Teodor Cirstoiu <tcirstoiu@torontomu.com>
//

#ifndef VISIBILITY_2_DEV_K_KERNEL_2_H
#define VISIBILITY_2_DEV_K_KERNEL_2_H

#include <CGAL/tags.h>
#include <CGAL/assertions.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

namespace CGAL {
    template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true>
    class K_kernel_2 {
    public:
        typedef Arrangement_2_ Arrangement_2;
        typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;

        typedef typename Geometry_traits_2::Point_2 Point_2;
        typedef typename Geometry_traits_2::Vector_2 Vector_2;
        typedef typename Geometry_traits_2::Direction_2 Direction_2;
        typedef typename Geometry_traits_2::Ray_2 Ray_2;
        typedef typename Geometry_traits_2::Line_2 Line_2;

        typedef RegularizationCategory Regularization_category;

        struct Wedge_ray {
            Direction_2 direction;
            std::size_t vertex; // index of the polygon vertex defining the line
            bool reversed; // true iff the ray points away from `vertex`

            Wedge_ray() : vertex(0), reversed(false) {
            }

            Wedge_ray(const Direction_2 &d, const std::size_t v, const bool r)
                : direction(d), vertex(v), reversed(r) {
            }
        };

        struct Wedge {
            Point_2 apex;
            std::size_t apex_vertex; // index of the apex (wedge bounding intersection) in the polygon
            Wedge_ray start; // clockwise-most bounding ray
            Wedge_ray end; // counterclockwise-most bounding ray

            Wedge() : apex_vertex(0) {
            }

            Wedge(const Point_2 &a, const std::size_t ai,
                  const Wedge_ray &s, const Wedge_ray &e)
                : apex(a), apex_vertex(ai), start(s), end(e) {
            }
        };

        typedef std::vector<Wedge> Wedge_container;

        template<class PointRange, class OutputIterator>
        OutputIterator wedges_at_vertex(const PointRange &polygon,
                                        std::size_t apex_index,
                                        OutputIterator oi) const {
            const std::vector<Point_2> points(std::begin(polygon), std::end(polygon));
            CGAL_precondition(apex_index < points.size());

            const Point_2 &apex = points[apex_index];

            // Both rays of the line through the apex and every other vertex.
            std::vector<Wedge_ray> rays;
            rays.reserve(2 * points.size());

            for (std::size_t i = 0; i < points.size(); ++i) {
                if (i == apex_index || points[i] == apex)
                    continue;
                const Direction_2 d(Vector_2(apex, points[i]));

                rays.push_back(Wedge_ray(d, i, false));
                rays.push_back(Wedge_ray(-d, i, true));
            }

            // Sort counterclockwise, starting at the positive x-axis, and merge the rays supported by collinear vertices.
            std::sort(rays.begin(), rays.end(), less_wedge_ray_comp());
            rays.erase(std::unique(rays.begin(), rays.end(), equal_wedge_ray_comp()),
                       rays.end());

            if (rays.size() < 2) {
                throw std::invalid_argument(
                    "Cannot construct wedges: Polygon has fewer than 2 distinct ray directions from apex."
                );
            }

            for (std::size_t i = 0; i < rays.size(); ++i) {
                const std::size_t j = (i + 1) % rays.size();
                *oi++ = Wedge(apex, apex_index, rays[i], rays[j]);
            }

            return oi;
        }

        template<class PointRange>
        void wedges(const PointRange &polygon,
                    std::vector<Wedge_container> &out) {
            const std::vector<Point_2> points(std::begin(polygon), std::end(polygon));

            out.clear();
            out.resize(points.size());

            // For each vertex build its wedges
            for (std::size_t i = 0; i < points.size(); ++i)
                wedges_at_vertex(points, i, std::back_inserter(out[i]));
        }

    private:
        struct less_wedge_ray_comp {
            bool operator()(const Wedge_ray &r1, const Wedge_ray &r2) const { return r1.direction < r2.direction; }
        };

        struct equal_wedge_ray_comp {
            bool operator()(const Wedge_ray &r1, const Wedge_ray &r2) const { return r1.direction == r2.direction; }
        };
    };
}

#endif //VISIBILITY_2_DEV_K_KERNEL_2_H
