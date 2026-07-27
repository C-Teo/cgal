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
#include <CGAL/Arrangement_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_vertical_decomposition_2.h>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include <set>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K> Point;
typedef K::Segment_2 Segment;
typedef K::Line_2 Line;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef K::FT FT; // Used for CGAL exact rationals in geometric tests
typedef CGAL::Arr_segment_traits_2<K> Traits;
typedef CGAL::Arrangement_2<Traits> Arrangement;

// This is what "Vert_decomp_list" actually is:
typedef Arrangement::Vertex_const_handle Vertex_handle;
typedef std::pair<CGAL::Object, CGAL::Object> Feature_pair;
typedef std::pair<Vertex_handle, Feature_pair> Decomp_result;
typedef std::vector<Decomp_result> Decomp_list;

// Struct to simplify our depth tree
struct LineIsect
{
  const Point& j_start;
  const Point& j_end;
  int depth;
};

void assign_depths(std::vector<LineIsect>& pairs);

std::tuple<Polygon, Polygon, std::set<Point>, std::set<Point>> build_pa_pb(const std::vector<LineIsect>& upper_pairs,
                                                                           const std::vector<LineIsect>& lower_pairs,
                                                                           const Polygon& polygon,
                                                                           FT ell_y,
                                                                           FT eps);

void process(const CGAL::Object& obj, const Point& p_prime, Arrangement& radial_arr, const Point& q, FT eps);
void add_radials(Arrangement& radial_arr,
                 const std::vector<Segment>& partition_edges,
                 const std::set<Point>& artificial,
                 const Polygon& p,
                 const Point& q,
                 FT eps);

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

  // Should do this in a loop until has_vertex_on_cut is false
  // But for testing one rotation should be enough
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
        // intersections.push_back(*pt);
        intersections.emplace_back(*pt);
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
  auto [pa, pb, pa_art, pb_art] = build_pa_pb(upper_pairs, lower_pairs, p, q.y(), eps);

  if(DEBUGGING) {
    std::cout << "Pa vertices:" << std::endl;
    for(auto vit = pa.vertices_begin(); vit != pa.vertices_end(); ++vit)
      std::cout << "  " << *vit << std::endl;
    std::cout << "Pb vertices:" << std::endl;
    for(auto vit = pb.vertices_begin(); vit != pb.vertices_end(); ++vit)
      std::cout << "  " << *vit << std::endl;
  }

  // Draw the scene for debugging purposes
  CGAL::Graphics_scene scene1;
  CGAL::Graphics_scene scene2;
  CGAL::Graphics_scene scene3;
  CGAL::Graphics_scene scene4;

  CGAL::add_to_graphics_scene(p, scene1);
  CGAL::add_to_graphics_scene(pa, scene2);
  CGAL::add_to_graphics_scene(pb, scene3);

  scene1.add_segment(cut_seg.source(), cut_seg.target());
  scene2.add_segment(cut_seg.source(), cut_seg.target());
  scene3.add_segment(cut_seg.source(), cut_seg.target());
  scene4.add_segment(cut_seg.source(), cut_seg.target());

  CGAL::draw_graphics_scene(scene1, "K-Visibility Test Scene1");
  CGAL::draw_graphics_scene(scene2, "K-Visibility Test Scene2");
  CGAL::draw_graphics_scene(scene3, "K-Visibility Test Scene3");

  // Lemma2 Theradial decomposition of a simple n-vertex polygon P around
  // a query point q can be computed in Θ(n)time.
  // Utilize CGAL::Arr_vertical_decomposition or CGAL::decompose
  Arrangement poly_arr;
  for(auto it = p.edges_begin(); it != p.edges_end(); ++it) {
    CGAL::insert(poly_arr, *it);
  }

  Arrangement radial_arr;
  add_radials(radial_arr, std::vector<Segment>(pa.edges_begin(), pa.edges_end()), pa_art, p, q, eps);
  add_radials(radial_arr, std::vector<Segment>(pb.edges_begin(), pb.edges_end()), pb_art, p, q, -eps);

  CGAL::add_to_graphics_scene(poly_arr, scene4);

  for(auto it = radial_arr.edges_begin(); it != radial_arr.edges_end(); ++it) {
    scene4.add_segment(it->curve().source(), it->curve().target());
  }

  CGAL::draw_graphics_scene(scene4, "K-Visibility Test Scene4");

  // Start the algorithm
  // . . .

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

std::tuple<Polygon, Polygon, std::set<Point>, std::set<Point>> build_pa_pb(const std::vector<LineIsect>& upper_pairs,
                                                                           const std::vector<LineIsect>& lower_pairs,
                                                                           const Polygon& polygon,
                                                                           FT ell_y,
                                                                           FT eps) {

  std::vector<Point> pa_pts;
  std::vector<Point> pb_pts;
  std::set<Point> pa_art;
  std::set<Point> pb_art;

  // Add any optimizations like a hash_map of intersection points to pairs for O(1) lookup
  // Keyed by x since all intersection points share y = ell_y
  std::map<FT, const LineIsect*> upper_at;
  for(const auto& pair : upper_pairs)
    upper_at[pair.j_start.x()] = &pair;

  std::map<FT, const LineIsect*> lower_at;
  for(const auto& pair : lower_pairs)
    lower_at[pair.j_start.x()] = &pair;

  std::vector<Point> verts(polygon.vertices_begin(), polygon.vertices_end());
  int n = (int)verts.size();

  // Find a point below the line to start the walk for polygon_b
  int start_b = 0;
  for(int i = 0; i < n; ++i) {
    if(verts[i].y() < ell_y) {
      start_b = i;
      break;
    }
  }

  // Find a point above the line to start the walk for polygon_a
  int start_a = 0;
  for(int i = 0; i < n; ++i) {
    if(verts[i].y() > ell_y) {
      start_a = i;
      break;
    }
  }

  // For the point below (polygon_b) follow the order of the polygon and add vertices to pb_pts
  // (as long as it is below the line)
  for(int j = 0; j < n; ++j) {
    int i = (start_b + j) % n;
    int i_next = (start_b + j + 1) % n;
    const Point& cur = verts[i];
    const Point& nxt = verts[i_next];

    if(cur.y() < ell_y)
      pb_pts.push_back(cur);

    // Check if moving to the next point hits an intersection
    bool cur_above = cur.y() > ell_y;
    bool nxt_above = nxt.y() > ell_y;
    if(cur_above == nxt_above)
      continue;

    // Linear interpolation along edge, find when y(t) = ell_y and plug for x
    FT t = (ell_y - cur.y()) / (nxt.y() - cur.y());
    Point isect(cur.x() + t * (nxt.x() - cur.x()), ell_y);

    if(!cur_above) {
      // If we hit an intersection coming from below, use the upper_pair starting at that
      // intersection to replace the boundary above l with the new three-edge path.
      // Append: isect, (isect.x, ell_y + eps/(2*d)), (j_end.x, ell_y + eps/(2*d))
      // Then continue following the polygon above the line until the next intersection.
      auto it = upper_at.find(isect.x());
      if(it != upper_at.end()) {
        FT len = eps / (2 * it->second->depth);
        Point isect_end(it->second->j_end.x(), ell_y);
        Point bridge1(isect.x(), ell_y + len);
        Point bridge2(it->second->j_end.x(), ell_y + len);
        pb_pts.push_back(isect);
        pb_pts.push_back(bridge1);
        pb_pts.push_back(bridge2);
        pb_pts.push_back(isect_end);
        pb_art.insert(isect);
        pb_art.insert(bridge1);
        pb_art.insert(bridge2);
        pb_art.insert(isect_end);
      }
    }
  }

  // For the point above (polygon_a) follow the order of the polygon and add vertices to pa_pts
  // (as long as it is above the line)
  for(int j = 0; j < n; ++j) {
    int i = (start_a + j) % n;
    int i_next = (start_a + j + 1) % n;
    const Point& cur = verts[i];
    const Point& nxt = verts[i_next];

    if(cur.y() > ell_y)
      pa_pts.push_back(cur);

    // Check if moving to the next point hits an intersection
    bool cur_above = cur.y() > ell_y;
    bool nxt_above = nxt.y() > ell_y;
    if(cur_above == nxt_above)
      continue;

    // Linear interpolation along edge, find when y(t) = ell_y and plug for x
    FT t = (ell_y - cur.y()) / (nxt.y() - cur.y());
    Point isect(cur.x() + t * (nxt.x() - cur.x()), ell_y);

    if(cur_above) {
      // If we hit an intersection coming from above, use the lower_pair starting at that
      // intersection to replace the boundary below l with the new three-edge path.
      // Append: isect, (isect.x, ell_y - eps/(2*d)), (j_end.x, ell_y - eps/(2*d))
      // Then continue following the polygon below the line until the next intersection.
      auto it = lower_at.find(isect.x());
      if(it != lower_at.end()) {
        FT len = eps / (2 * it->second->depth);
        Point isect_end(it->second->j_end.x(), ell_y);
        Point bridge1(isect.x(), ell_y - len);
        Point bridge2(it->second->j_end.x(), ell_y - len);
        pa_pts.push_back(isect);
        pa_pts.push_back(bridge1);
        pa_pts.push_back(bridge2);
        pa_pts.push_back(isect_end);
        pa_art.insert(isect);
        pa_art.insert(bridge1);
        pa_art.insert(bridge2);
        pa_art.insert(isect_end);
      }
    }
  }

  return {Polygon(pa_pts.begin(), pa_pts.end()), Polygon(pb_pts.begin(), pb_pts.end()), pa_art, pb_art};
}

// FQ: Projective Transformation Mapping a Point to Infinity
// https://www.youtube.com/watch?v=E-mLLId3uuY

// In homogeneous coordinates we represent a point as [x, y, w]
// To find the 2D Cartesian result we divide by the third component: (x/w, y/w)

// First we translate point + epsilon so q is (0,0) (for the lower polygon all points are below the line)
// Then we apply the projection matrix M which sends (0,1,0) to (0,0,1) (the point at infinity in the y direction)

// M = | 1 0 0 |
//     | 0 0 1 |
//     | 0 1 0 |

// Which is simplified to (x, y) -> (x/y, 1/y) in Cartesian coordinates
Point fq(const Point& p, const Point& q, FT eps) {
  // Translate
  FT dx = p.x() - q.x();
  FT dy = p.y() - q.y() + eps;

  // Project
  return Point(dx / dy, FT(1) / dy);
}

// To inverse we apply the inverse of the projection matrix M^-1 which sends (0,0,1) back to (0,1,0)
// And then we translate back by -epsilon and the location of q
// The inverse matrix M^-1 is the same as M since it is symmetric and involutory (M = M^-1)
Point fq_inv(const Point& p_prime, const Point& q, FT eps) {
  // Project
  FT x = p_prime.x() / p_prime.y();
  FT y = FT(1) / p_prime.y();

  // Translate
  return Point(x + q.x(), y + q.y() - eps);
}

void add_radials(Arrangement& radial_arr,
                 const std::vector<Segment>& partition_edges,
                 const std::set<Point>& artificial,
                 const Polygon& p,
                 const Point& q,
                 FT eps) {
  Arrangement trans_arr;

  for(const auto& e : partition_edges) {
    Point src = fq(e.source(), q, eps);
    Point tgt = fq(e.target(), q, eps);

    CGAL::insert(trans_arr, Segment(src, tgt));
  }

  // 2. Perform vertical decomposition on the unrolled shape
  Decomp_list vd_list;
  CGAL::decompose(trans_arr, std::back_inserter(vd_list));

  if(DEBUGGING) {
    CGAL::Graphics_scene scene_trans;
    CGAL::add_to_graphics_scene(trans_arr, scene_trans);

    // AI helped with debug here...
    for(const auto& entry : vd_list) {
      Point p_prime = entry.first->point();
      Arrangement::Halfedge_const_handle he;
      Arrangement::Vertex_const_handle vh;

      auto draw_ray = [&](const CGAL::Object& obj) {
        if(CGAL::assign(vh, obj)) {
          scene_trans.add_segment(p_prime, vh->point());
        } else if(CGAL::assign(he, obj)) {
          Point s = he->curve().source(), t = he->curve().target();
          FT u = p_prime.x();
          FT y_int = s.y() + (t.y() - s.y()) * (u - s.x()) / (t.x() - s.x());
          scene_trans.add_segment(p_prime, Point(u, y_int));
        }
      };

      draw_ray(entry.second.first);  // below
      draw_ray(entry.second.second); // above
    }

    CGAL::draw_graphics_scene(scene_trans, "Transformed Arrangement + Decomp Rays");
  }

  // 3. Process rays
  for(const auto& entry : vd_list) {
    // Each entry is:
    // - Point p' where the ray from q hits the arrangement (in transformed space)
    // - A pair of CGAL::Objects representing the feature hit above and below p'
    Point p_prime = entry.first->point();
    Point p = fq_inv(p_prime, q, eps);

    if(artificial.count(p)) {
      std::cout << "Skipping artificial point: (" << CGAL::to_double(p.x()) << ", " << CGAL::to_double(p.y()) << ")"
                << std::endl;
      continue; // three-edge construction point
    }

    process(entry.second.second, p_prime, radial_arr, q, eps);
    process(entry.second.first, p_prime, radial_arr, q, eps);
  }
}

void process(const CGAL::Object& obj, const Point& p_prime, Arrangement& radial_arr, const Point& q, FT eps) {
  Point p = fq_inv(p_prime, q, eps);

  Arrangement::Halfedge_const_handle he;
  Arrangement::Vertex_const_handle vh;

  if(CGAL::assign(vh, obj)) { // Ray hits a vertex

    CGAL::insert(radial_arr, Segment(p, fq_inv(vh->point(), q, eps)));

  } else if(CGAL::assign(he, obj)) { // Ray hits a edge
    // u = x-coordinate of the shooting vertex (the vertical ray)
    // s, t = endpoints of the hit edge
    Point s = he->curve().source();
    Point t = he->curve().target();
    FT u = p_prime.x();

    // Interpolate in transformed space
    FT alpha = (u - s.x()) / (t.x() - s.x());
    FT y_int = s.y() + alpha * (t.y() - s.y());

    // Pull back the intersection point to original space
    CGAL::insert(radial_arr, Segment(p, fq_inv(Point(u, y_int), q, eps)));
  } else {
    return; // Ray hit nothing
  }
}

/*

// 1. compute ell_y and eps
FT ell_y = q.y();
FT eps = ...;

// 2. compute jordan sequence (intersection points in boundary order)
std::vector<Point> jordan = ...;

// 3. fill upper_pairs and lower_pairs
std::vector<LineIsect> upper_pairs = ...;
std::vector<LineIsect> lower_pairs = ...;

// 4. assign depths
assign_depths(upper_pairs);
assign_depths(lower_pairs);

// 5. build Pa and Pb
auto [Pa, Pb] = build_pa_pb(P, jordan, upper_pairs, lower_pairs, ell_y, eps);

// 6. build radial decompositions
Arr arr_a = build_radial_decomposition(Pa, q);
Arr arr_b = build_radial_decomposition(Pb, q);

// 7. merge
Arr radial = merge_radial_decompositions(arr_a, arr_b);

*/
