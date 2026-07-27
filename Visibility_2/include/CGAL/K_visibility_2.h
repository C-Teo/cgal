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

#ifndef CGAL_K_VISIBILITY_REGION
#define CGAL_K_VISIBILITY_REGION

#include <CGAL/license/Visibility_2.h>

#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions_2.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <stack>
#include <tuple>
#include <vector>

namespace CGAL {

template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true >
class K_visibility_region_2 {

public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;

  typedef typename K::FT                                FT;
  typedef typename K::Point_2                           Point_2;
  typedef typename K::Segment_2                         Segment_2;
  typedef typename K::Vector_2                          Vector_2;
  typedef CGAL::Polygon_2<K>                            Polygon_2;
  typedef CGAL::Aff_transformation_2<K>                 Aff_transformation_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;

  typedef RegularizationCategory               Regularization_category;
  typedef CGAL::Tag_false                      Supports_general_polygon_category;
  typedef CGAL::Tag_true                       Supports_simple_polygon_category;

  // This is what "Vert_decomp_list" actually is:
  typedef std::pair<CGAL::Object, CGAL::Object>          Feature_pair;
  typedef std::pair<Vertex_const_handle, Feature_pair>   Decomp_result;
  typedef std::vector<Decomp_result>                     Decomp_list;

  typedef std::set<Point_2>                              Point_set;

  // Two consecutive crossings of the sweep line l by the boundary of P, delimiting one excursion of
  // the boundary into a half-plane, together with the nesting depth of that excursion (1 for an
  // outermost one). Used to simplify our depth tree.
  struct Line_isect
  {
    Point_2 j_start;
    Point_2 j_end;
    int depth = 0;
  };

  /*! Builds the radial decomposition of `polygon` around `q` into `radial_arr`.
   *
   * `polygon` has to be simple and in general position with respect to the horizontal line l
   * through `q`, i.e. no vertex may lie on l -- see rotate_to_general_position().
   */
  static void radial_decomposition(const Polygon_2& polygon, const Point_2& q, Arrangement_2& radial_arr) {
    const FT ell_y = q.y();

    // Minimum distance from any vertex to l: inside the band of height eps around l the boundary of
    // P consists of nothing but the edges crossing l, which is what the cuts below rely on.
    const FT eps = min_vertical_distance(polygon, ell_y);

    std::vector<Point_2> crossings = line_crossings(polygon, ell_y);

    std::vector<Line_isect> upper_pairs;
    std::vector<Line_isect> lower_pairs;
    split_pairs(crossings, upper_pairs, lower_pairs);

    assign_depths(upper_pairs);
    assign_depths(lower_pairs);

    Polygon_2 pa, pb;
    Point_set pa_art, pb_art;
    std::tie(pa, pb, pa_art, pb_art) = build_pa_pb(upper_pairs, lower_pairs, polygon, ell_y, eps);

    add_radials(radial_arr, std::vector<Segment_2>(pa.edges_begin(), pa.edges_end()), pa_art, q);
    add_radials(radial_arr, std::vector<Segment_2>(pb.edges_begin(), pb.edges_end()), pb_art, q);
  }

  /*! Rotates `polygon` about `q` by a small rational angle until no vertex lies on the horizontal
   * line through `q`.
   *
   * The transformation is a similarity rather than an exact rotation (the rational sine and cosine
   * scale the polygon by a factor of sqrt(1 + 1/1000^2) per step), which leaves the boundary simple
   * and the visibility combinatorics unchanged.
   */
  static void rotate_to_general_position(Polygon_2& polygon, const Point_2& q) {
    // | cos θ   -sin θ |     | c/hw   -s/hw |
    // | sin θ    cos θ |  =  | s/hw    c/hw |
    Aff_transformation_2 to_origin(CGAL::TRANSLATION, Vector_2(-q.x(), -q.y()));
    Aff_transformation_2 rotation(CGAL::ROTATION, 1, 1000, 1000);
    Aff_transformation_2 from_origin(CGAL::TRANSLATION, Vector_2(q.x(), q.y()));
    Aff_transformation_2 transform = from_origin * rotation * to_origin;

    while(has_vertex_on_line(polygon, q.y()))
      polygon = CGAL::transform(transform, polygon);
  }

  static bool has_vertex_on_line(const Polygon_2& polygon, FT ell_y) {
    for(auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
      if(it->y() == ell_y)
        return true;
    }
    return false;
  }

  /*! Minimum vertical distance from a vertex of `polygon` to the line y = `ell_y`. */
  static FT min_vertical_distance(const Polygon_2& polygon, FT ell_y) {
    CGAL_precondition(!polygon.is_empty());

    FT eps = CGAL::abs(polygon.vertices_begin()->y() - ell_y); // initialize to first vertex
    for(auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
      FT dist = CGAL::abs(it->y() - ell_y); // Compute strictly vertical distance
      if(dist < eps)
        eps = dist;
    }
    return eps;
  }

  /*! The Jordan sequence: the points where the boundary of `polygon` crosses the line y = `ell_y`,
   * in boundary order, starting at a crossing where the boundary leaves the lower half-plane.
   */
  static std::vector<Point_2> line_crossings(const Polygon_2& polygon, FT ell_y) {
    std::vector<Point_2> crossings;
    int first_upward = -1;

    const int n = (int)polygon.size();
    for(int i = 0; i < n; ++i) {
      const Point_2& cur = polygon.vertex(i);
      const Point_2& nxt = polygon.vertex((i + 1) % n);

      CGAL_precondition(cur.y() != ell_y); // general position, see rotate_to_general_position()

      bool cur_above = cur.y() > ell_y;
      bool nxt_above = nxt.y() > ell_y;
      if(cur_above == nxt_above)
        continue;

      if(!cur_above && first_upward < 0)
        first_upward = (int)crossings.size();

      crossings.push_back(edge_point_at_y(cur, nxt, ell_y));
    }

    // split_pairs() needs the boundary to run above l between the first two crossings, so start the
    // sequence where it leaves the lower half-plane. (Starting at the leftmost crossing instead only
    // gives that for a counterclockwise oriented polygon.)
    if(first_upward > 0)
      std::rotate(crossings.begin(), crossings.begin() + first_upward, crossings.end());

    return crossings;
  }

  /*! Splits a Jordan sequence into the excursions above and below l.
   *
   * Between consecutive pairs (x_{2i-1}, x_{2i}) of the Jordan sequence, for i in {1, ..., m/2}, the
   * polygon boundary of P lies above l. Similarly, between pairs (x_{2j}, x_{2j+1}), for
   * j in {1, ..., m/2 - 1}, and between (x_m, x_0), the boundary of P lies below l.
   */
  static void split_pairs(const std::vector<Point_2>& crossings,
                          std::vector<Line_isect>& upper_pairs,
                          std::vector<Line_isect>& lower_pairs) {
    const int m = (int)crossings.size();
    CGAL_precondition(m % 2 == 0); // the boundary is closed, so it crosses l an even number of times

    if(m == 0)
      return;

    for(int index = 1; index <= m / 2; ++index) {
      const Point_2& x2i_1 = crossings[(2 * index - 1) - 1];
      const Point_2& x2i = crossings[(2 * index) - 1];
      upper_pairs.push_back({x2i_1, x2i, 0});
    }

    for(int index = 1; index <= m / 2 - 1; ++index) {
      const Point_2& x2j = crossings[(2 * index) - 1];
      const Point_2& x2j_1 = crossings[(2 * index + 1) - 1];
      lower_pairs.push_back({x2j, x2j_1, 0});
    }

    lower_pairs.push_back({crossings.back(), crossings.front(), 0});
  }

  // Note: partially overlapping intervals can never appear in a valid simple polygon
  static void assign_depths(std::vector<Line_isect>& arcs) {
    if(arcs.empty())
      return;

    struct Iv
    {
      FT lo, hi;
      int index;
    };

    std::vector<Iv> intervals;
    intervals.reserve(arcs.size());

    for(int i = 0; i < (int)arcs.size(); ++i) {
      FT xa = arcs[i].j_start.x();
      FT xb = arcs[i].j_end.x();
      intervals.push_back({(CGAL::min)(xa, xb), (CGAL::max)(xa, xb), i});
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
      CGAL_assertion_msg(intervals[i - 1].hi <= intervals[i].lo || intervals[i - 1].hi >= intervals[i].hi,
                         "Partially overlapping intervals - polygon may not be simple");
    }

    // Stack tracks the hi endpoints of every interval that is still "open"
    std::stack<FT> active;

    for(const Iv& iv : intervals) {
      while(!active.empty() && active.top() <= iv.lo) {
        //  By the nested parenthesis property, any future interval must also start after it ended
        //  (since intervals are processed left to right and the future one's lo is even larger).
        active.pop();
      }

      arcs[iv.index].depth = (int)active.size() + 1;
      active.push(iv.hi);
    }
  }

  /*! Distance from l at which the excursion of the given depth gets cut off.
   *
   * The cut has to stay inside the band of height eps around l (no vertex of P is closer than eps to
   * l, so inside that band the boundary consists only of the edges crossing l), and a nested
   * excursion has to be cut *further* from l than the excursion enclosing it: the crossing edges of
   * the nested pair reach up to their own cut, so a nested cut closer to l would be sliced by the
   * chord of the enclosing pair. The offset therefore grows with the depth instead of shrinking,
   * while staying strictly inside (0, eps).
   */
  static FT cut_offset(int depth, FT eps) {
    return eps - eps / (2 * depth);
  }

  /*! Point where the edge (`cur` -> `nxt`) crosses the horizontal line y = `target_y`.
   *
   * With target_y = ell_y -/+ cut_offset(...) this places a point on the edge that approached l, at
   * exactly that offset from l, without ever crossing to the other side. Trigonometrically it is
   * stepping back offset / sin(theta) along the edge, theta being the angle the edge makes with l;
   * since l is horizontal that step reduces to a difference in y, so solving for the edge parameter
   * stays exact and the distance from the point to l is a straight y check.
   */
  static Point_2 edge_point_at_y(const Point_2& cur, const Point_2& nxt, FT target_y) {
    FT t = (target_y - cur.y()) / (nxt.y() - cur.y());
    return Point_2(cur.x() + t * (nxt.x() - cur.x()), target_y);
  }

  /*! Splits `polygon` into the sub-polygon Pa above l and the sub-polygon Pb below it.
   *
   * Each excursion of the boundary into the opposite half-plane is replaced by a single horizontal
   * chord that stops short of l, at cut_offset(depth) from it. Returns Pa, Pb and, for each of them,
   * the set of chord endpoints, which are artefacts of the construction rather than vertices of P.
   */
  static std::tuple<Polygon_2, Polygon_2, Point_set, Point_set> build_pa_pb(
      const std::vector<Line_isect>& upper_pairs,
      const std::vector<Line_isect>& lower_pairs,
      const Polygon_2& polygon,
      FT ell_y,
      FT eps) {

    std::vector<Point_2> pa_pts;
    std::vector<Point_2> pb_pts;
    Point_set pa_art;
    Point_set pb_art;

    // Add any optimizations like a hash_map of intersection points to pairs for O(1) lookup
    // Keyed by x since all intersection points share y = ell_y. Both endpoints of a pair map to it,
    // so a crossing resolves to its pair whether the walk enters or leaves the excursion there.
    std::map<FT, const Line_isect*> upper_at;
    for(const auto& pair : upper_pairs) {
      upper_at[pair.j_start.x()] = &pair;
      upper_at[pair.j_end.x()] = &pair;
    }

    std::map<FT, const Line_isect*> lower_at;
    for(const auto& pair : lower_pairs) {
      lower_at[pair.j_start.x()] = &pair;
      lower_at[pair.j_end.x()] = &pair;
    }

    std::vector<Point_2> verts(polygon.vertices_begin(), polygon.vertices_end());
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
      const Point_2& cur = verts[i];
      const Point_2& nxt = verts[i_next];

      if(cur.y() < ell_y)
        pb_pts.push_back(cur);

      // Check if moving to the next point hits an intersection
      bool cur_above = cur.y() > ell_y;
      bool nxt_above = nxt.y() > ell_y;
      if(cur_above == nxt_above)
        continue;

      // Linear interpolation along edge, find when y(t) = ell_y and plug for x
      Point_2 isect = edge_point_at_y(cur, nxt, ell_y);

      // The walk is about to leave the lower half-plane here, or to come back into it: the part of
      // the boundary running above l is replaced by a single chord that stays below l. Stop short of
      // the intersection, on the edge that approached it, at cut_offset(d) below l. The pair's other
      // crossing is reached later in the same walk and stops at the same offset, so the two points
      // become consecutive and the chord between them is the horizontal cut.
      auto it = upper_at.find(isect.x());
      if(it == upper_at.end())
        continue;

      Point_2 stop = edge_point_at_y(cur, nxt, ell_y - cut_offset(it->second->depth, eps));
      pb_pts.push_back(stop);
      pb_art.insert(stop);
    }

    // For the point above (polygon_a) follow the order of the polygon and add vertices to pa_pts
    // (as long as it is above the line)
    for(int j = 0; j < n; ++j) {
      int i = (start_a + j) % n;
      int i_next = (start_a + j + 1) % n;
      const Point_2& cur = verts[i];
      const Point_2& nxt = verts[i_next];

      if(cur.y() > ell_y)
        pa_pts.push_back(cur);

      // Check if moving to the next point hits an intersection
      bool cur_above = cur.y() > ell_y;
      bool nxt_above = nxt.y() > ell_y;
      if(cur_above == nxt_above)
        continue;

      // Linear interpolation along edge, find when y(t) = ell_y and plug for x
      Point_2 isect = edge_point_at_y(cur, nxt, ell_y);

      // Mirror of the walk above: the boundary running below l is replaced by a chord that stays
      // above l, both of its endpoints sitting on the edges of the lower pair at cut_offset(d).
      auto it = lower_at.find(isect.x());
      if(it == lower_at.end())
        continue;

      Point_2 stop = edge_point_at_y(cur, nxt, ell_y + cut_offset(it->second->depth, eps));
      pa_pts.push_back(stop);
      pa_art.insert(stop);
    }

    return {Polygon_2(pa_pts.begin(), pa_pts.end()), Polygon_2(pb_pts.begin(), pb_pts.end()), pa_art, pb_art};
  }

  // FQ: Projective Transformation Mapping a Point to Infinity
  // https://www.youtube.com/watch?v=E-mLLId3uuY

  // In homogeneous coordinates we represent a point as [x, y, w]
  // To find the 2D Cartesian result we divide by the third component: (x/w, y/w)

  // First we translate so q is (0,0)
  // Then we apply the projection matrix M which sends (0,1,0) to (0,0,1) (the point at infinity in
  // the y direction)

  // M = | 1 0 0 |
  //     | 0 0 1 |
  //     | 0 1 0 |

  // Which is simplified to (x, y) -> (x/y, 1/y) in Cartesian coordinates

  // No epsilon shift is needed: Pa lies strictly above l and Pb strictly below it (their cuts stop
  // at cut_offset(d) from l), so dy never vanishes for a point of either polygon.
  static Point_2 fq(const Point_2& p, const Point_2& q) {
    // Translate
    FT dx = p.x() - q.x();
    FT dy = p.y() - q.y();

    // Project
    return Point_2(dx / dy, FT(1) / dy);
  }

  // To inverse we apply the inverse of the projection matrix M^-1 which sends (0,0,1) back to
  // (0,1,0) and then we translate back by the location of q
  // The inverse matrix M^-1 is the same as M since it is symmetric and involutory (M = M^-1)
  static Point_2 fq_inv(const Point_2& p_prime, const Point_2& q) {
    // Project
    FT x = p_prime.x() / p_prime.y();
    FT y = FT(1) / p_prime.y();

    // Translate
    return Point_2(x + q.x(), y + q.y());
  }

  /*! Adds the radials of one half of the partition to `radial_arr`.
   *
   * Lemma 2: the radial decomposition of a simple n-vertex polygon P around a query point q can be
   * computed in Theta(n) time. Under fq the rays through q become vertical, so the radials are
   * obtained from a vertical decomposition of the transformed edges.
   */
  static void add_radials(Arrangement_2& radial_arr,
                          const std::vector<Segment_2>& partition_edges,
                          const Point_set& artificial,
                          const Point_2& q) {
    // 1. Transform the partition so that the rays through q run vertically
    Arrangement_2 trans_arr;

    for(const auto& e : partition_edges) {
      Point_2 src = fq(e.source(), q);
      Point_2 tgt = fq(e.target(), q);

      CGAL::insert(trans_arr, Segment_2(src, tgt));
    }

    // 2. Perform vertical decomposition on the unrolled shape
    Decomp_list vd_list;
    CGAL::decompose(trans_arr, std::back_inserter(vd_list));

    // 3. Process rays
    for(const auto& entry : vd_list) {
      // Each entry is:
      // - Point p' where the ray from q hits the arrangement (in transformed space)
      // - A pair of CGAL::Objects representing the feature hit above and below p'
      Point_2 p_prime = entry.first->point();
      Point_2 p = fq_inv(p_prime, q);

      if(artificial.count(p))
        continue; // chord endpoint of the cut construction, not a vertex of P

      add_radial(entry.second.second, p_prime, radial_arr, q);
      add_radial(entry.second.first, p_prime, radial_arr, q);
    }
  }

  /*! Adds the radial from `p_prime` to the feature `obj` it sees, pulled back to original space. */
  static void add_radial(const CGAL::Object& obj, const Point_2& p_prime, Arrangement_2& radial_arr, const Point_2& q) {
    Point_2 p = fq_inv(p_prime, q);

    Halfedge_const_handle he;
    Vertex_const_handle vh;

    if(CGAL::assign(vh, obj)) { // Ray hits a vertex

      CGAL::insert(radial_arr, Segment_2(p, fq_inv(vh->point(), q)));

    } else if(CGAL::assign(he, obj)) { // Ray hits a edge
      // u = x-coordinate of the shooting vertex (the vertical ray)
      // s, t = endpoints of the hit edge
      Point_2 s = he->curve().source();
      Point_2 t = he->curve().target();
      FT u = p_prime.x();

      // Interpolate in transformed space
      FT alpha = (u - s.x()) / (t.x() - s.x());
      FT y_int = s.y() + alpha * (t.y() - s.y());

      // Pull back the intersection point to original space
      CGAL::insert(radial_arr, Segment_2(p, fq_inv(Point_2(u, y_int), q)));
    } else {
      return; // Ray hit nothing
    }
  }
};

} // namespace CGAL

#endif
