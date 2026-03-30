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
  Point q(5, 3.5);
  // Point q(5, 2);

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
  double eps = std::numeric_limits<double>::max();
  for(auto it = p.vertices_begin(); it != p.vertices_end(); ++it) {
    // Compute strictly vertical distance
    double dist = std::abs(CGAL::to_double(it->y()) - CGAL::to_double(q.y()));

    if(dist < eps) {
      eps = dist;
    }
  }

  if(DEBUGGING)
    std::cout << "Epsilon: " << eps << std::endl;

  // Create edges (results in Polygon_A and Polygon_B)

  // Convert Polygon_A and Polygon_B to Arrangement_2

  // Utilize CGAL Arr_vertical_decomposition or Arr_trapezoid_ric_point_location

  // Draw the scene for debugging purposes
  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(p, scene);
  scene.add_segment(cut_seg.source(), cut_seg.target());
  CGAL::draw_graphics_scene(scene, "K-Visibility Test");

  return EXIT_SUCCESS;
}
