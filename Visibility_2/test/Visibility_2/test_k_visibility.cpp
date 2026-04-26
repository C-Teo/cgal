// Copyright (c) 2026 Toronto Metropolitan Unversity (Canada).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Yeganeh Bahoo Torudi <bahoo@torontomu.ca>
//             Teodor Cirstoiu <tcirstoiu@torontomu.ca>
//

/*
cmake -S . \
-B build/release/ \
-DCMAKE_BUILD_TYPE=Release \
-DCGAL_DIR=~/projects/cgal
*/

#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K> Point;
typedef K::Segment_2 Segment;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef K::FT FT; // Used for CGAL exact rationals in geometric tests

// Struct to simplify our depth tree
struct LineIsect
{
  const Point& j_start;
  const Point& j_end;
  int depth;
};

void assign_depths(std::vector<LineIsect>& pairs);

std::pair<Polygon, Polygon> build_pa_pb(const std::vector<LineIsect>& upper_pairs,
                                        const std::vector<LineIsect>& lower_pairs,
                                        const std::vector<Point>& intersections,
                                        FT ell_y,
                                        FT eps) {}

bool DEBUGGING = true;

int main() {
  // create a polygon and put some points in it
  Polygon p;

  p.push_back(Point(0, 0));
  p.push_back(Point(0, 5));
  p.push_back(Point(2, 5));
  p.push_back(Point(2, 2));
  p.push_back(Point(8, 2));
  p.push_back(Point(8, 4));
  p.push_back(Point(6, 4));
  p.push_back(Point(6, 3));
  p.push_back(Point(4, 3));
  p.push_back(Point(4, 5));
  p.push_back(Point(10, 5));
  p.push_back(Point(10, 0));

  // create our agent position
  Point q(5, 3.5); // Point q(5, 2);

  // create a line segment to cut through the polygon and q (horizontal sweep line)
  Segment cut_seg(Point(-2, q.y()), Point(12, q.y()));

  // Rotate p around q until no vertex lies on cut_seg
  bool has_vertex_on_cut = false;

  for(auto vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
    if(cut_seg.has_on(*vit)) {
      has_vertex_on_cut = true;
      break;
    }
  }

  if(has_vertex_on_cut) {
    // Apply a small rotation around q: sin=1/1000, cos=1000/1000
    // | cos θ   -sin θ |     | c/hw   -s/hw |
    // | sin θ    cos θ |  =  | s/hw    c/hw |

    Transformation to_origin(CGAL::TRANSLATION, K::Vector_2(-q.x(), -q.y()));
    Transformation rotation(CGAL::ROTATION, 1, 1000, 1000);
    Transformation from_origin(CGAL::TRANSLATION, K::Vector_2(q.x(), q.y()));
    Transformation transform = from_origin * rotation * to_origin;
    p = CGAL::transform(transform, p);
  }

  if(DEBUGGING) {
    for(auto vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
      std::cout << "Vertex: " << *vit << std::endl;
    }
  }

  // loop through all edges in the polygon
  // - extract the geometric segment and return the intersection if it exists
  std::vector<Point> intersections;

  for(auto eit = p.edges_begin(); eit != p.edges_end(); ++eit) {
    const Segment& edge_seg = *eit;

    if(DEBUGGING)
      std::cout << "Edge from " << edge_seg.source() << " to " << edge_seg.target() << std::endl;

    auto result = CGAL::intersection(edge_seg, cut_seg);

    if(result.has_value()) {
      // assume no segment overlap
      if(const Point* pt = std::get_if<Point>(&*result)) {
        if(DEBUGGING)
          std::cout << "Intersection point: " << *pt << std::endl;
        intersections.push_back(*pt);
      }
    }
  }

  // Find the leftmost intersection and move to the front WITHOUT changing order
  if(!intersections.empty()) {
    auto leftmost_it = std::min_element(intersections.begin(), intersections.end(),
                                        [](const Point& a, const Point& b) { return a.x() < b.x(); });
    std::rotate(intersections.begin(), leftmost_it, intersections.end());
  }

  // Compute min distance epsilon from any vertex to the line
  FT eps = CGAL::abs(p.vertices_begin()->y() - q.y()); // initialize to first vertex
  for(auto it = p.vertices_begin(); it != p.vertices_end(); ++it) {
    FT dist = CGAL::abs(it->y() - q.y()); // Compute strictly vertical distance
    if(dist < eps)
      eps = dist;
  }

  if(DEBUGGING)
    std::cout << "Epsilon: " << eps << std::endl;

  // Create edges (results in Polygon_A and Polygon_B)

  // -- Sort by upper and lower pairs

  // Between consecutive pairs (x2i−1,x2i) of the Jordan sequence, for i ∈{1,...,
  // m/2}, the polygon boundary of P lies above . Similarly, between pairs (x2j,x2j+1),
  // for j ∈{1,...,m/2−1}, and between (xm,x0), the boundary of P lies below .
  int m = (int)intersections.size();

  if(DEBUGGING)
    std::cout << "Number of intersections: " << m << std::endl;

  std::vector<LineIsect> upper_pairs;

  for(int index = 1; index <= m / 2; ++index) {
    const Point& x2i_1 = intersections[(2 * index - 1) - 1];
    const Point& x2i = intersections[(2 * index) - 1];
    upper_pairs.push_back({x2i_1, x2i});
  }

  std::vector<LineIsect> lower_pairs;

  for(int index = 1; index <= m / 2 - 1; ++index) {
    const Point& x2j = intersections[(2 * index) - 1];
    const Point& x2j_1 = intersections[(2 * index + 1) - 1];
    lower_pairs.push_back({x2j, x2j_1});
  }

  lower_pairs.push_back({intersections.back(), intersections.front()});

  // -- Create tree of all depths and then sort
  assign_depths(upper_pairs);
  assign_depths(lower_pairs);

  if(DEBUGGING) {
    std::cout << "Upper pairs:" << std::endl;
    for(const auto& pair : upper_pairs) {
      std::cout << "  (" << pair.j_start << ", " << pair.j_end << ")" << ", " << pair.depth << std::endl;
    }

    std::cout << "Lower pairs:" << std::endl;
    for(const auto& pair : lower_pairs) {
      std::cout << "  (" << pair.j_start << ", " << pair.j_end << ")" << ", " << pair.depth << std::endl;
    }
  }

  // -- Create new three-way edges

  // -- Append top/bottom points with its respective cut edges

  // Convert Polygon_A and Polygon_B to Arrangement_2

  // Lemma2 Theradial decomposition of a simple n-vertex polygon P around
  // a query point q can be computed in Θ(n)time.
  // Utilize CGAL Arr_vertical_decomposition or Arr_trapezoid_ric_point_location

  // Start the algorithm...

  // Draw the scene for debugging purposes
  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(p, scene);
  scene.add_segment(cut_seg.source(), cut_seg.target());
  CGAL::draw_graphics_scene(scene, "K-Visibility Test");

  return EXIT_SUCCESS;
}

// Note: partially overlapping intervals can never appear in a valid simple polygon
void assign_depths(std::vector<LineIsect>& arcs) {
  if(arcs.empty())
    return;

  struct Iv
  {
    double lo, hi;
    int index;
  };

  std::vector<Iv> intervals;
  intervals.reserve(arcs.size());

  for(int i = 0; i < (int)arcs.size(); ++i) {
    double xa = CGAL::to_double(arcs[i].j_start.x());
    double xb = CGAL::to_double(arcs[i].j_end.x());
    intervals.push_back({std::min(xa, xb), std::max(xa, xb), i});
  }

  // Sort by lo ascending, and if tie, by hi descending
  std::sort(intervals.begin(), intervals.end(), [](const Iv& a, const Iv& b) {
    if(a.lo != b.lo)
      return a.lo < b.lo;
    return a.hi > b.hi;
  });

  // After sorting, check no two intervals partially overlap
  for(int i = 1; i < (int)intervals.size(); ++i) {
    // valid cases: disjoint (prev.hi <= cur.lo) or nested (prev.hi >= cur.hi)
    assert(intervals[i - 1].hi <= intervals[i].lo ||
           intervals[i - 1].hi >= intervals[i].hi && "Partially overlapping intervals — polygon may not be simple");
  }

  // Stack tracks the hi endpoints of every interval that is still "open"
  std::stack<double> active;

  for(const Iv& iv : intervals) {
    while(!active.empty() && active.top() <= iv.lo) {
      //  By the nested parenthesis property, any future interval must also start after it ended (since intervals are
      //  processed left to right and the future one's lo is even larger).
      active.pop();
    }

    arcs[iv.index].depth = (int)active.size() + 1;
    active.push(iv.hi);
  }
}

std::pair<Polygon, Polygon> build_pa_pb(const std::vector<LineIsect>& upper_pairs,
                                        const std::vector<LineIsect>& lower_pairs,
                                        const std::vector<Point>& intersections,
                                        FT ell_y,
                                        FT eps) {}