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
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions_2.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace CGAL {
    template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true>
    class K_kernel_2 {
    public:
        typedef Arrangement_2_ Arrangement_2;
        typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;

        typedef typename Geometry_traits_2::FT FT;
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

        enum Wedge_side {
            WEDGE_RIGHT, // closed sector between -end and start, clockwise from A
            WEDGE_LEFT, // closed sector between end and -start, counterclockwise from A
            WEDGE_INTERIOR, // strictly inside A or A~
            WEDGE_APEX // the apex itself
        };

        // Vertices are tagged as lying in AL (left region) or AR (right region)
        // Tests a point against wedge.start.direction and wedge.end.direction (the two lines through the apex)
        Wedge_side side_of_wedge(const Wedge &wedge, const Point_2 &p) const {
            const Vector_2 v(wedge.apex, p);

            if (v == NULL_VECTOR)
                return WEDGE_APEX;

            const Orientation o_start = orientation(wedge.start.direction.vector(), v);
            const Orientation o_end = orientation(wedge.end.direction.vector(), v);

            // The tests are negated so that a point collinear with a bounding
            // ray (i.e. a vertex defining it) is assigned to the sector it bounds, not to the interior.
            if (o_start != LEFT_TURN && o_end != LEFT_TURN) return WEDGE_RIGHT;
            if (o_start != RIGHT_TURN && o_end != RIGHT_TURN) return WEDGE_LEFT;
            return WEDGE_INTERIOR; // left of one line and right of the other: inside A or A~
        }

        struct Clipping_edge {
            Point_2 source; // endpoint on the right of the wedge
            Point_2 target; // endpoint on the left of the wedge
            std::size_t edge; // index i of the polygon edge (p_i, p_{i+1})

            // the parameter at which the edge's supporting line crosses the wedge axis
            FT param_num;
            FT param_den;

            Clipping_edge() : edge(0), param_num(0), param_den(1) {
            }

            Clipping_edge(const Point_2 &s, const Point_2 &t, const std::size_t e,
                          const FT &num, const FT &den)
                : source(s), target(t), edge(e), param_num(num), param_den(den) {
            }

            Sign signed_distance_sign() const { return CGAL::sign(param_num); }
        };

        // E(A), with `list[i - 1]` denoting E(A)_i.
        typedef std::vector<Clipping_edge> Clipping_list;

        template<class PointRange>
        Clipping_list create_clipping_list(const PointRange &polygon, const Wedge &wedge) const {
            Clipping_list list;

            const std::vector<Point_2> points(std::begin(polygon), std::end(polygon));
            const std::size_t n = points.size();

            CGAL_precondition(n > 2);
            CGAL_precondition(wedge.apex_vertex < n);
            CGAL_precondition(points[wedge.apex_vertex] == wedge.apex);

            const std::size_t apex = wedge.apex_vertex;
            const std::size_t prev = (apex + n - 1) % n; // s-
            const std::size_t next = (apex + 1) % n; // s+

            // Used to order the elements of E(A) according to their signed distance from s
            const Vector_2 axis = wedge_axis(wedge);

            // Map every vertex i to its relative location with respect to the wedge
            std::vector<Wedge_side> sides(n);
            for (std::size_t i = 0; i < n; ++i)
                sides[i] = side_of_wedge(wedge, points[i]);

            CGAL_precondition(sides[prev] != WEDGE_APEX && sides[next] != WEDGE_APEX);

            // Always append (s-,s) and conditionally append (s,s+)
            append_clipping_edge(points, sides, prev, apex, wedge.apex, axis, list);
            if (sides[prev] == sides[next])
                append_clipping_edge(points, sides, apex, next, wedge.apex, axis, list);

            // Build the rest of the clipping edges
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t j = (i + 1) % n;

                if (i == apex || j == apex)
                    continue; // incident to the apex, already handled above

                CGAL_precondition_msg(sides[i] != WEDGE_INTERIOR && sides[j] != WEDGE_INTERIOR,
                                      "A vertex lies strictly inside the wedge: the wedge was not "
                                      "built from the lines through its apex and the vertices of "
                                      "this polygon.");

                // Check if we should append the edge (if it clips), and if so append to output list
                append_clipping_edge(points, sides, i, j, wedge.apex, axis, list);
            }

            // Sort E(A) by their signed distance relative to the apex s
            std::stable_sort(list.begin(), list.end(), less_clipping_edge_comp());

            return list;
        }

        template<class PointRange, class OutputIterator>
        OutputIterator create_clipping_list(const PointRange &polygon, const Wedge &wedge,
                                            OutputIterator oi) const {
            const Clipping_list list = create_clipping_list(polygon, wedge);
            return std::copy(list.begin(), list.end(), oi);
        }

    private:
        struct less_wedge_ray_comp {
            bool operator()(const Wedge_ray &r1, const Wedge_ray &r2) const { return r1.direction < r2.direction; }
        };

        struct equal_wedge_ray_comp {
            bool operator()(const Wedge_ray &r1, const Wedge_ray &r2) const { return r1.direction == r2.direction; }
        };

        struct less_clipping_edge_comp {
            bool operator()(const Clipping_edge &e1, const Clipping_edge &e2) const {
                return e1.param_num * e2.param_den < e2.param_num * e1.param_den;
            }
        };

        // A ray strictly inside A, used as a ruler: every edge of E(A) crosses it exactly once,
        // and the crossing position (`param_num`/`param_den`) orders E(A) by distance from the apex.
        Vector_2 wedge_axis(const Wedge &wedge) const {
            const Vector_2 start = wedge.start.direction.vector();
            const Vector_2 end = wedge.end.direction.vector();

            CGAL_precondition_msg(orientation(start, end) == LEFT_TURN,
                                  "Degenerate wedge: its bounding directions span a half-plane.");

            return start + end;
        }

        void append_clipping_edge(const std::vector<Point_2> &points,
                                  const std::vector<Wedge_side> &sides,
                                  const std::size_t a, const std::size_t b,
                                  const Point_2 &apex, const Vector_2 &axis,
                                  Clipping_list &list) const {
            const bool a_right = (sides[a] == WEDGE_RIGHT || sides[a] == WEDGE_APEX);
            const bool a_left = (sides[a] == WEDGE_LEFT || sides[a] == WEDGE_APEX);
            const bool b_right = (sides[b] == WEDGE_RIGHT || sides[b] == WEDGE_APEX);
            const bool b_left = (sides[b] == WEDGE_LEFT || sides[b] == WEDGE_APEX);

            if (a_right && b_left)
                list.push_back(make_clipping_edge(points[a], points[b], a, apex, axis));
            else if (a_left && b_right)
                list.push_back(make_clipping_edge(points[b], points[a], a, apex, axis));
        }

        // Should never need to run `make_clipping_edge` unless we are going to append it
        Clipping_edge make_clipping_edge(const Point_2 &source, const Point_2 &target,
                                         const std::size_t edge, const Point_2 &apex,
                                         const Vector_2 &axis) const {
            const Vector_2 d(source, target);

            const FT num = determinant(Vector_2(apex, source), d);
            const FT den = determinant(axis, d);

            CGAL_precondition_msg(orientation(d, axis) == RIGHT_TURN,
                                  "Clipping edge does not run from the right sector to the left one.");

            return Clipping_edge(source, target, edge, num, den);
        }
    };
}

#endif //VISIBILITY_2_DEV_K_KERNEL_2_H
