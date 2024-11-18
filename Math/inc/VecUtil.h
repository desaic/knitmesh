#include <sstream>
#include <string>

#include "Vec2.h"
#include "Vec3.h"

template <typename T>
std::string to_string(const Vec3<T>& v) {
  std::stringstream ss;
  ss << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
  return ss.str();
}

template <typename T>
std::string to_string(const Vec2<T>& v) {
  std::stringstream ss;
  ss << "[" << v[0] << ", " << v[1] << "]";
  return ss.str();
}

bool InBound(int ix, int iy, int iz, const Vec3u& size);
bool InBound(int ix, int iy, int iz, const Vec2u& size);

bool InBoundu(const Vec3u& gridIdx, const Vec3u& gridSize);
bool InBoundu(const Vec2u& gridIdx, const Vec2u& gridSize);

Vec3i ToSigned(const Vec3u& v);
Vec3u ToUnsigned(const Vec3i& v);

Vec3i Trunc(const Vec3f& v);

Vec3f VecMax(const Vec3f& a, const Vec3f& b);
Vec3f VecMin(const Vec3f& a, const Vec3f& b);

Vec3f clamp(const Vec3f& v, const Vec3f& lb, const Vec3f& ub);
