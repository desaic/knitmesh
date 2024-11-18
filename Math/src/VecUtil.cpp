#include "VecUtil.h"
#include <algorithm>

bool InBound(int ix, int iy, int iz, const Vec3u& size) {
  return (ix >= 0) && (iy >= 0) && (iz >= 0) && ix < int(size[0]) && iy < int(size[1]) &&
         iz < int(size[2]);
}

bool InBoundu(const Vec3u& gridIdx, const Vec3u& gridSize) {
  return gridIdx[0] < gridSize[0] && gridIdx[1] < gridSize[1] && gridIdx[2] < gridSize[2];
}

bool InBound(int ix, int iy, int iz, const Vec2u& size) {
  return (ix >= 0) && (iy >= 0) && (iz >= 0) && ix < int(size[0]) && iy < int(size[1]);
}

bool InBoundu(const Vec2u& gridIdx, const Vec2u& gridSize) {
  return gridIdx[0] < gridSize[0] && gridIdx[1] < gridSize[1];
}

Vec3i ToSigned(const Vec3u& v) { return Vec3i(int(v[0]), int(v[1]), int(v[2])); }
Vec3u ToUnsigned(const Vec3i& v) { return Vec3u(unsigned(v[0]), unsigned(v[1]), unsigned(v[2])); }

Vec3i Trunc(const Vec3f& v) { return Vec3i(int(v[0]), int(v[1]), int(v[2])); }

Vec3f VecMax(const Vec3f& a, const Vec3f& b) {
  Vec3f c;
  c[0] = std::max(a[0], b[0]);
  c[1] = std::max(a[1], b[1]);
  c[2] = std::max(a[2], b[2]);
  return c;
}

Vec3f VecMin(const Vec3f& a, const Vec3f& b) {
  Vec3f c;
  c[0] = std::min(a[0], b[0]);
  c[1] = std::min(a[1], b[1]);
  c[2] = std::min(a[2], b[2]);
  return c;
}

Vec3f clamp(const Vec3f& v, const Vec3f& lb, const Vec3f& ub) {
  return Vec3f(std::clamp(v[0], lb[0], ub[0]), std::clamp(v[1], lb[1], ub[1]),
    std::clamp(v[2], lb[2], ub[2]));
}
