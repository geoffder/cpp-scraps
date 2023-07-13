#include "SlowHull.h"
#include <unordered_set>

double DistanceToVector(const glm::vec3 &pt, const glm::vec3 &v) {
  return glm::length(pt - v * glm::dot(pt, v));
}

inline double DistanceToPlane(const glm::vec3 &pt, const glm::vec4 &plane) {
  return glm::dot({plane.x, plane.y, plane.z}, pt) - plane.w;
}

glm::vec4 PlaneOfPts(const glm::vec3 &a, const glm::vec3 &b,
                     const glm::vec3 &c) {
  auto crx = glm::cross(c - a, b - a);
  auto n = glm::length(crx);
  if (n == 0.)
    return glm::vec4();
  return glm::vec4(crx, glm::dot(crx, a)) / n;
}

// Return a trio of non-collinear indices into pts. If none are found, indices
// will be -1 (and hulling cannot be performed).
glm::ivec3 NonCollinearTriple(const std::vector<glm::vec3> &pts,
                              const double precision) {
  const int len = pts.size();
  int furthest = 1;
  float dist = glm::distance(pts[0], pts[1]);
  for (int i = 2; i < len; i++) {
    double d = glm::distance(pts[0], pts[i]);
    if (d > dist) {
      furthest = i;
      dist = d;
    }
  }
  if (dist <= precision)
    return glm::ivec3(-1);
  const auto n = (pts[0] - pts[furthest]) / dist;
  int third = -1;
  double offset = dist * precision;
  for (int i = 1; i < len; i++) {
    auto off = DistanceToVector(pts[i] - pts[0], n);
    if (off > offset) {
      third = i;
      offset = off;
    }
  }
  if (third < 0)
    return glm::ivec3(-1);
  return glm::ivec3(0, furthest, third);
}

std::pair<int, bool> NonCoplanar(const std::vector<glm::vec3> &pts,
                                 const glm::vec4 &plane,
                                 const double precision) {
  for (int i = 0; i < pts.size(); i++) {
    auto dist = DistanceToPlane(pts[i], plane);
    if (std::abs(dist) > precision)
      return {i, dist > precision};
  }
  return {-1, false}; // All points are coplanar, hull is not a 3d shape
}

struct IVec2Hash {
  std::size_t operator()(const glm::ivec2 &v) const {
    return (v.y << 16) ^ v.x;
  }
};

std::vector<glm::ivec3> Hull(std::vector<glm::vec3> pts,
                             const double precision) {
  const int len = pts.size();
  // FIXME: Should there be a hull specific error?
  if (len < 4)
    return std::vector<glm::ivec3>{};
  const auto trip = NonCollinearTriple(pts, precision);
  const auto plane = PlaneOfPts(pts[trip.x], pts[trip.y], pts[trip.z]);
  const auto [d, dIsAboveTrip] = NonCoplanar(pts, plane, precision);
  if (d < 0)
    return std::vector<glm::ivec3>{}; // All points coplanar.
  const int a = trip.x;
  const int b = dIsAboveTrip ? trip.z : trip.y;
  const int c = dIsAboveTrip ? trip.y : trip.z;
  std::vector<glm::ivec3> triangles;
  std::vector<glm::vec4> planes;
  std::unordered_set<glm::ivec2, IVec2Hash> halfEdges{};
  std::vector<int> dropped;
  std::vector<bool> kept;
  const auto AddTri = [&pts, &triangles, &planes, &dropped, &kept](int a, int b,
                                                                   int c) {
    auto tri = glm::ivec3(c, b, a);
    auto plane = PlaneOfPts(pts[a], pts[b], pts[c]);
    if (dropped.size() > 0) {
      const int idx = dropped[dropped.size() - 1];
      triangles[idx] = tri;
      planes[idx] = plane;
      kept[idx] = true;
      dropped.pop_back();
    } else {
      triangles.push_back(tri);
      planes.push_back(plane);
      kept.push_back(true);
    }
  };
  AddTri(a, b, c);
  AddTri(d, b, a);
  AddTri(c, d, a);
  AddTri(b, d, c);
  for (int i = 0; i < len; i++) {
    if (i == a || i == b || i == c || i == d)
      continue; // skip starting points
    // collect half edges of triangles that are in conflict with the points at
    // idx, pruning the conflicting triangles and their planes in the process
    for (int j = 0; j < triangles.size(); j++) {
      if (kept[j] && DistanceToPlane(pts[i], planes[j]) > precision) {
        halfEdges.insert({triangles[j].x, triangles[j].z});
        halfEdges.insert({triangles[j].z, triangles[j].y});
        halfEdges.insert({triangles[j].y, triangles[j].x});
        // mark triangle/plane as dropped (slot to be reused)
        dropped.push_back(j);
        kept[j] = false;
      }
    }
    // form new triangles with the outer perimeter (horizon) of the set of
    // conflicting triangles and the point at idx
    for (auto e : halfEdges) {
      if (halfEdges.erase({e.y, e.x}) == 0)
        AddTri(e.x, e.y, i);
    }
    halfEdges.clear();
  }
  std::vector<glm::ivec3> out;
  out.reserve((triangles.size() - dropped.size()) * 3);
  for (int i = 0; i < triangles.size(); i++) {
    if (!kept[i])
      continue;
    out.push_back(triangles[i]);
  }
  return out;
}

std::vector<glm::ivec3> Hull(const std::vector<glm::vec3> &pts) {
  float scale = 0;
  for (int i = 0; i < pts.size(); i++) {
    auto abs = glm::abs(pts[i]);
    if (abs.x > scale)
      scale = abs.x;
    if (abs.y > scale)
      scale = abs.y;
    if (abs.z > scale)
      scale = abs.z;
  }
  return Hull(pts, 1e-9 * scale);
}
