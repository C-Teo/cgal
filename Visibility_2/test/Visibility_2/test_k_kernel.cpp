// Copyright (c) 2026 Toronto Metropolitan University (Canada).
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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/K_kernel_2.h>

#include <cstddef>
#include <iostream>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::K_kernel_2<Arrangement_2> K_kernel;

typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Direction_2 Direction_2;

typedef K_kernel::Wedge Wedge;
typedef K_kernel::Wedge_container Wedge_container;
typedef K_kernel::Clipping_edge Clipping_edge;
typedef K_kernel::Clipping_list Clipping_list;

static int failures = 0;

// Kept independent of NDEBUG so the checks also run in a release build.
#define CHECK(cond)                                                            \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::cerr << "  FAILED (line " << __LINE__ << "): " #cond << '\n';   \
            ++failures;                                                          \
        }                                                                        \
    } while (false)

// A simple counterclockwise polygon mixing spikes with reflex chains, and with no three vertices
// collinear, so every wedge is bounded by exactly one vertex direction on each side.
static std::vector<Point_2> test_polygon() {
    static const int coords[][2] = {
        {47, 447}, {172, 114}, {330, 227}, {476, 227}, {468, 97},
        {897, 255}, {711, 453}, {694, 502}, {628, 290}, {458, 395},
        {600, 572}, {222, 508}, {170, 288}
    };
    constexpr std::size_t n = std::size(coords);

    std::vector<Point_2> polygon;
    polygon.reserve(n);
    for (const auto coord: coords)
        polygon.emplace_back(coord[0], coord[1]);

    return polygon;
}

// Mirrors the private `wedge_axis`, which is just the sum of the two bounding directions.
static Vector_2 axis_of(const Wedge &wedge) {
    return wedge.start.direction.vector() + wedge.end.direction.vector();
}

// The crossing position on the axis, stored as a fraction to keep it division-free.
static FT param_of(const Clipping_edge &e) {
    return e.param_num / e.param_den;
}

// The wedges at each vertex partition the directions around it into consecutive sub-half-plane cones.
static void test_wedges(const std::vector<Point_2> &polygon,
                        const std::vector<Wedge_container> &all) {
    constexpr K_kernel kk;
    const std::size_t n = polygon.size();

    // One wedge container per vertex.
    CHECK(all.size() == n);

    for (std::size_t i = 0; i < n; ++i) {
        const Wedge_container &wedges = all[i];

        // 12 other vertices and no three collinear: 12 distinct lines, hence 24 distinct rays.
        CHECK(wedges.size() == 2 * (n - 1));

        for (std::size_t j = 0; j < wedges.size(); ++j) {
            const Wedge &w = wedges[j];

            // Every wedge is anchored at the vertex it was built for.
            CHECK(w.apex == polygon[i]);
            CHECK(w.apex_vertex == i);

            // Bounding rays are supported by real polygon vertices, never by the apex itself.
            CHECK(w.start.vertex < n && w.start.vertex != i);
            CHECK(w.end.vertex < n && w.end.vertex != i);

            // Each wedge spans strictly less than a half-plane, which `wedge_axis` requires.
            CHECK(CGAL::orientation(w.start.direction.vector(),
                w.end.direction.vector()) == CGAL::LEFT_TURN);

            // Consecutive wedges share a bounding ray, so together they cover the full turn.
            CHECK(w.end.direction == wedges[(j + 1) % wedges.size()].start.direction);

            // `reversed` says whether the ray points at its defining vertex or away from it.
            const Direction_2 to_start(Vector_2(w.apex, polygon[w.start.vertex]));
            const Direction_2 to_end(Vector_2(w.apex, polygon[w.end.vertex]));
            CHECK(w.start.direction == (w.start.reversed ? -to_start : to_start));
            CHECK(w.end.direction == (w.end.reversed ? -to_end : to_end));
        }
    }

    // The single-vertex overload agrees with the batch one.
    std::vector<Wedge> single;
    kk.wedges_at_vertex(polygon, 0, std::back_inserter(single));
    CHECK(single.size() == all[0].size());
    for (std::size_t j = 0; j < single.size() && j < all[0].size(); ++j) {
        CHECK(single[j].start.direction == all[0][j].start.direction);
        CHECK(single[j].end.direction == all[0][j].end.direction);
    }
}

// Points are classified into the two closed outer sectors, the interior of A / A~, and the apex.
static void test_side_of_wedge(const std::vector<Point_2> &polygon,
                               const std::vector<Wedge_container> &all) {
    const std::size_t n = polygon.size();

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto &w: all[i]) {
            constexpr K_kernel kk;
            const Vector_2 axis = axis_of(w);

            // The apex is its own case, and belongs to no sector.
            CHECK(kk.side_of_wedge(w, w.apex) == K_kernel::WEDGE_APEX);
            CHECK(!kk.is_in_wedge(w, w.apex));

            // A and A~ occupy the same pair of half-planes, so only `is_in_wedge` separates them.
            const Point_2 inside = w.apex + axis;
            const Point_2 opposite = w.apex - axis;
            CHECK(kk.side_of_wedge(w, inside) == K_kernel::WEDGE_INTERIOR);
            CHECK(kk.side_of_wedge(w, opposite) == K_kernel::WEDGE_INTERIOR);
            CHECK(kk.is_in_wedge(w, inside));
            CHECK(!kk.is_in_wedge(w, opposite));

            // The sectors are closed: a vertex defining a bounding ray belongs to the sector it bounds.
            CHECK(kk.side_of_wedge(w, polygon[w.start.vertex]) ==
                (w.start.reversed ? K_kernel::WEDGE_LEFT : K_kernel::WEDGE_RIGHT));
            CHECK(kk.side_of_wedge(w, polygon[w.end.vertex]) ==
                (w.end.reversed ? K_kernel::WEDGE_RIGHT : K_kernel::WEDGE_LEFT));

            // No vertex ever lands strictly inside a wedge: `create_clipping_list` requires it.
            for (std::size_t v = 0; v < n; ++v)
                CHECK(kk.side_of_wedge(w, polygon[v]) != K_kernel::WEDGE_INTERIOR);
        }
    }
}

// E(A) holds exactly the polygon edges spanning the double wedge, ordered by distance from the apex.
static void test_clipping_list(const std::vector<Point_2> &polygon,
                               const std::vector<Wedge_container> &all) {
    const std::size_t n = polygon.size();

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto &w: all[i]) {
            constexpr K_kernel kk;
            const Vector_2 axis = axis_of(w);
            const Clipping_list list = kk.create_clipping_list(polygon, w);

            // The edge (s-, s) is appended unconditionally, so E(A) is never empty.
            CHECK(!list.empty());

            std::size_t at_apex = 0;
            for (std::size_t e = 0; e < list.size(); ++e) {
                const Clipping_edge &ce = list[e];

                // The stored index refers to a real polygon edge.
                CHECK(ce.edge < n);

                // The denominator is normalized positive; the comparator and
                // `signed_distance_sign()` are only correct because of it.
                CHECK(ce.param_den > 0);
                CHECK(ce.signed_distance_sign() == CGAL::sign(param_of(ce)));

                // `source` is the right-sector endpoint and `target` the left-sector one.
                const K_kernel::Wedge_side s = kk.side_of_wedge(w, ce.source);
                const K_kernel::Wedge_side t = kk.side_of_wedge(w, ce.target);
                CHECK(s == K_kernel::WEDGE_RIGHT || s == K_kernel::WEDGE_APEX);
                CHECK(t == K_kernel::WEDGE_LEFT || t == K_kernel::WEDGE_APEX);

                // Those two endpoints are exactly the endpoints of polygon edge (p_e, p_{e+1}).
                const Point_2 &a = polygon[ce.edge];
                const Point_2 &b = polygon[(ce.edge + 1) % n];
                CHECK((ce.source == a && ce.target == b) || (ce.source == b && ce.target == a));

                // The parameter really locates where the supporting line crosses the axis.
                CHECK(CGAL::collinear(ce.source, ce.target, w.apex + param_of(ce) * axis));

                if (param_of(ce) == 0)
                    ++at_apex;

                if (e > 0)
                    CHECK(param_of(list[e - 1]) <= param_of(ce)); // sorted by signed distance
            }

            // The edges incident to s cross the axis at the apex itself, at parameter zero.
            CHECK(at_apex >= 1);

            // The output-iterator overload yields the same list.
            Clipping_list copied;
            kk.create_clipping_list(polygon, w, std::back_inserter(copied));
            CHECK(copied.size() == list.size());
            for (std::size_t e = 0; e < copied.size() && e < list.size(); ++e)
                CHECK(copied[e].edge == list[e].edge && param_of(copied[e]) == param_of(list[e]));
        }
    }
}

// A point of A is k-clipped exactly when it lies past E(A)_{k+2}, strictly so for even k.
static void test_is_k_clipped(const std::vector<Point_2> &polygon,
                              const std::vector<Wedge_container> &all) {
    const std::size_t n = polygon.size();

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto &w: all[i]) {
            constexpr K_kernel kk;
            const Vector_2 axis = axis_of(w);
            const Clipping_list list = kk.create_clipping_list(polygon, w);

            std::vector<FT> params;
            params.reserve(list.size());
            for (const auto &e: list)
                params.push_back(param_of(e));

            // The distinct crossings on the A side of the axis, apex included as the first one.
            std::vector<FT> knots;
            knots.emplace_back(0);
            for (const auto &param: params)
                if (param > 0 && param != knots.back())
                    knots.push_back(param);

            // Sample between consecutive crossings, on each crossing (where the orientation is
            // COLLINEAR and the parity rule decides), and beyond the last one.
            std::vector<FT> samples;
            for (std::size_t e = 0; e + 1 < knots.size(); ++e)
                samples.push_back((knots[e] + knots[e + 1]) / 2);
            for (std::size_t e = 1; e < knots.size(); ++e)
                samples.push_back(knots[e]);
            samples.push_back(knots.back() + 1);

            for (const auto &t: samples) {
                const Point_2 x = w.apex + t * axis;

                CHECK(kk.is_in_wedge(w, x)); // every sample has t > 0, so it lies in A

                for (std::size_t k = 0; k < list.size() + 2; ++k) {
                    // x sits on the axis at t and E(A)_{k+2} crosses it at params[k+1], so x is
                    // right of that edge exactly when t is the larger parameter -- strictly for
                    // even k, with equality allowed for odd k.
                    const bool expected = (k + 1 < params.size()) &&
                                          (k % 2 == 0 ? t > params[k + 1] : t >= params[k + 1]);
                    CHECK(kk.is_k_clipped(w, list, x, k) == expected);
                }
            }

            // Nothing outside A is k-clipped: not A~, not the apex, not the bounding vertices.
            const Point_2 outside[] = {
                w.apex - axis, w.apex, polygon[w.start.vertex], polygon[w.end.vertex]
            };
            for (const auto &p: outside)
                for (std::size_t k = 0; k < list.size() + 2; ++k)
                    CHECK(!kk.is_k_clipped(w, list, p, k));
        }
    }
}

int main() {
    const std::vector<Point_2> polygon = test_polygon();

    K_kernel kk;
    std::vector<Wedge_container> all;
    kk.wedges(polygon, all);

    std::cout << "Testing wedges / wedges_at_vertex..." << std::endl;
    test_wedges(polygon, all);

    std::cout << "Testing side_of_wedge / is_in_wedge..." << std::endl;
    test_side_of_wedge(polygon, all);

    std::cout << "Testing create_clipping_list..." << std::endl;
    test_clipping_list(polygon, all);

    std::cout << "Testing is_k_clipped..." << std::endl;
    test_is_k_clipped(polygon, all);

    if (failures != 0) {
        std::cerr << failures << " check(s) failed." << std::endl;
        return 1;
    }

    std::cout << "All checks passed." << std::endl;
    return 0;
}
