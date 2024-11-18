#ifndef VEC_3_H
#define VEC_3_H

#define _USE_MATH_DEFINES
#include <array>
#include <cmath>
#include <ostream>

template <typename T>
class Vec3 {
 public:
  Vec3() {
    v[0] = (T)0;
    v[1] = (T)0;
    v[2] = (T)0;
  }
  Vec3(T val) {
    v[0] = val;
    v[1] = val;
    v[2] = val;
  }
  Vec3(T v1, T v2, T v3) {
    v[0] = v1;
    v[1] = v2;
    v[2] = v3;
  }
  Vec3(std::array<T, 3> arr) {
    v[0] = arr[0];
    v[1] = arr[1];
    v[2] = arr[2];
  }
  const T& operator[](unsigned i) const { return v[i]; }
  T& operator[](unsigned i) { return v[i]; }

  Vec3<T>& operator+=(const Vec3<T>& v1) {
    v[0] += v1[0];
    v[1] += v1[1];
    v[2] += v1[2];
    return *this;
  }

  Vec3<T>& operator-=(const Vec3<T>& v1) {
    v[0] -= v1[0];
    v[1] -= v1[1];
    v[2] -= v1[2];
    return *this;
  }

  Vec3<T>& operator*=(T v1) {
    v[0] *= v1;
    v[1] *= v1;
    v[2] *= v1;
    return *this;
  }

  Vec3<T>& operator/=(T v1) {
    T inv = (T)(1.0 / v1);
    v[0] *= inv;
    v[1] *= inv;
    v[2] *= inv;
    return *this;
  }

  Vec3<T>& operator/=(Vec3<T> v1) {
    v[0] /= v1[0];
    v[1] /= v1[1];
    v[2] /= v1[2];
    return *this;
  }

  Vec3<T> operator+(const Vec3<T>& v1) {
    auto ret = Vec3<T>(v[0] + v1[0], v[1] + v1[1], v[2] + v1[2]);
    return ret;
  }

  Vec3<T> operator-(const Vec3<T>& v1) {
    auto ret = Vec3<T>(v[0] - v1[0], v[1] - v1[1], v[2] - v1[2]);
    return ret;
  }

  Vec3<T> operator*(T v1) {
    auto ret = Vec3<T>(v[0] * v1, v[1] * v1, v[2] * v1);
    return ret;
  }

  Vec3<T> operator*(const Vec3<T>& v1) const {
    auto ret = Vec3<T>(v[0] * v1[0], v[1] * v1[1], v[2] * v1[2]);
    return ret;
  }

  Vec3<T>& operator*=(const Vec3<T>& v1) {
    v[0] *= v1[0];
    v[1] *= v1[1];
    v[2] *= v1[2];
    return *this;
  }

  Vec3<T> operator/(T v1) {
    auto ret = Vec3<T>(v[0] / v1, v[1] / v1, v[2] / v1);
    return ret;
  }

  Vec3<T> operator/(const Vec3<T>& v1) {
    auto ret = Vec3<T>(v[0] / v1[0], v[1] / v1[1], v[2] / v1[2]);
    return ret;
  }

  bool operator==(const Vec3<T>& v1) const { return v[0] == v1[0] && v[1] == v1[1] && v[2] == v1[2]; }
  T dot(const Vec3<T>& b) const {
    T prod = v[0] * b[0] + v[1] * b[1] + v[2] * b[2];
    return prod;
  }
  Vec3<T> operator-() const { return Vec3<T>(-v[0], -v[1], -v[2]); }

  Vec3<T> cross(const Vec3<T>& b) const {
    Vec3<T> c(v[1] * b[2] - v[2] * b[1], -v[0] * b[2] + v[2] * b[0], v[0] * b[1] - v[1] * b[0]);
    return c;
  }

  T norm2() const { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; }

  T norm() const { return (T)(std::sqrt(float(norm2()))); }

  Vec3<T> absval() const{
    return Vec3<T>(abs(v[0]), abs(v[1]), abs(v[2]));
  }

  void normalize() {
    T n = norm();
    if (n > 0) {
      n = 1.0f / n;
      v[0] *= n;
      v[1] *= n;
      v[2] *= n;
    }
  }

  Vec3 normalizedCopy()const {
    T n = norm();
    if (n > 0) {
      n = 1.0f / n;
    }
    return Vec3(v[0] * n, v[1] * n, v[2] * n);
  }

  void setZero() {
    v[0] = 0.f;
    v[1] = 0.f;
    v[2] = 0.f;
  }

  T* data() {
    return &v[0];
  }

  template<typename U>
  Vec3<U> cast() const {
    return {
        static_cast<U>(v[0]),
        static_cast<U>(v[1]),
        static_cast<U>(v[2]),
    };
  }
  std::array<T, 3> data_array() const {
    return {v[0], v[1], v[2]};
  }

 private:
  T v[3];
};

template <typename T>
std::ostream &operator<<(std::ostream &os, Vec3<T> const &m) {
  return os << m[0] << " " << m[1] << " " << m[2];
}

template <typename T>
Vec3<T> operator+(const Vec3<T>& v0, const Vec3<T>& v1) {
  return Vec3<T>(v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]);
}

template <typename T>
Vec3<T> operator-(const Vec3<T>& v0, const Vec3<T>& v1) {
  return Vec3<T>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]);
}

/// scalar multiplication
template <typename T>
inline Vec3<T> operator*(T s, const Vec3<T>& v) {
  return Vec3<T>(v[0] * s, v[1] * s, v[2] * s);
}

template <typename T>
inline Vec3<T> cWiseMin(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(
    std::min(v1[0], v2[0]),
    std::min(v1[1], v2[1]),
    std::min(v1[2], v2[2])
  );
}

template <typename T>
inline Vec3<T> cWiseMax(const Vec3<T>& v1, const Vec3<T>& v2) {
  return Vec3<T>(
    std::max(v1[0], v2[0]),
    std::max(v1[1], v2[1]),
    std::max(v1[2], v2[2])
  );
}

typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<unsigned> Vec3u;
typedef Vec3<unsigned char> Vec3b;
typedef Vec3<unsigned char> Vec3uc;
typedef Vec3<int> Vec3i;
typedef Vec3<unsigned long long> Vec3ul;

#endif
