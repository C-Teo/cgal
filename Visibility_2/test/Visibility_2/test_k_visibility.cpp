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
#include <CGAL/Arrangement_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <iterator>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Point_2<K>                                    Point;
typedef CGAL::Arr_segment_traits_2<K>                       Traits;
typedef CGAL::Arrangement_2<Traits>                         Arrangement;

typedef Traits::Segment_2                       Segment;

typedef Arrangement::Halfedge_handle                       Halfedge_handle;
typedef Arrangement::Face_handle                           Face_handle;
typedef Arrangement::Vertex_handle                         Vertex_handle;

typedef std::variant<Vertex_handle, Halfedge_handle, Face_handle>
                                                        Zone_result;

int main() {
  // create a polygon and put some points in it
  Polygon_2 p;

  p.push_back(Point(0, 0));
  p.push_back(Point(0, 5));
  p.push_back(Point(2, 5));
  p.push_back(Point(2, 2));
  p.push_back(Point(8, 2));
  p.push_back(Point(8, 5));
  p.push_back(Point(10, 5));
  p.push_back(Point(10, 0));

  Point q(5,2.5);

  Arrangement arr;
  CGAL::insert(arr, p.edges_begin(), p.edges_end());

  Segment cut(Point(-2, q.y()), Point(12, q.y()));

  std::list<Zone_result> zone_elems;

  zone(arr, cut, std::back_inserter(zone_elems));

  std::size_t zone_actual_comp = zone_elems.size();

  std::cout << zone_actual_comp << std::endl;
  
for (const auto& elem : zone_elems)
{
    switch (elem.index())
    {
        case 0: std::cout << "Vertex_handle\n"; break;
        case 1: std::cout << "Halfedge_handle\n"; break;
        case 2: std::cout << "Face_handle\n"; break;
    }
}

for (const auto& elem : zone_elems) {
  if (const Halfedge_handle* he = std::get_if<Halfedge_handle>(&elem)) {

    Segment edge_seg(
      Point((*he)->source()->point().x(), (*he)->source()->point().y()),
      Point((*he)->target()->point().x(), (*he)->target()->point().y()));

    Segment query_seg(
      Point(cut.source().x(), cut.source().y()),
      Point(cut.target().x(), cut.target().y()));

    auto result = CGAL::intersection(query_seg, edge_seg);

    if (result) {
                if (const Point* p = std::get_if<Point>(&*result))
                {
                    std::cout << "  Edge hit! Intersection at: ("
                              << p->x() << ", " << p->y() << ")\n";
                }
    }
  }
}

  CGAL::draw(arr);

  return EXIT_SUCCESS;
}

