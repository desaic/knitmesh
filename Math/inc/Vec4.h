#ifndef VEC_4_H
#define VEC_4_H

#define _USE_MATH_DEFINES
#include <cmath>

template <typename T>
class Vec4 {
 public:
  Vec4() {
    v[0] = (T)0;
    v[1] = (T)0;
    v[2] = (T)0;
    v[3] = (T)0;
  }
  Vec4(T val) {
    v[0] = val;
    v[1] = val;
    v[2] = val;
    v[3] = val;
  }
  Vec4(T v1, T v2, T v3, T v4) {
    v[0] = v1;
    v[1] = v2;
    v[2] = v3;
    v[3] = v4;
  }

  const T& operator[](unsigned i) const { return v[i]; }
  
  T& operator[](unsigned i) { return v[i]; }

  Vec4<T>& operator+=(const Vec4<T>& v1) {
    v[0] += v1[0];
    v[1] += v1[1];
    v[2] += v1[2];
    v[3] += v1[3];
    return *this;
  }

  Vec4<T>& operator-=(const Vec4<T>& v1) {
    v[0] -= v1[0];
    v[1] -= v1[1];
    v[2] -= v1[2];
    v[3] -= v1[3];
    return *this;
  }

  Vec4<T>& operator*=(T v1) {
    v[0] *= v1;
    v[1] *= v1;
    v[2] *= v1;
    v[3] *= v1;
    return *this;
  }

  Vec4<T>& operator/=(T v1) {
    T inv = (T)(1.0 / v1);
    v[0] *= inv;
    v[1] *= inv;
    v[2] *= inv;
    v[3] *= inv;
    return *this;
  }

  Vec4<T> operator+(const Vec4<T>& v1) {
    auto ret = Vec4<T>(v[0] + v1[0], v[1] + v1[1], v[2] + v1[2], v[3] + v1[3]);
    return ret;
  }

  Vec4<T> operator-(const Vec4<T>& v1) {
    auto ret = Vec4<T>(v[0] - v1[0], v[1] - v1[1], v[2] - v1[2], v[3] - v1[3]);
    return ret;
  }

  Vec4<T> operator*(T v1) {
    auto ret = Vec4<T>(v[0] * v1, v[1] * v1, v[2] * v1, v[3] * v1);
    return ret;
  }

  Vec4<T> operator/(T v1) {
    auto ret = Vec4<T>(v[0] / v1, v[1] / v1, v[2] / v1, v[3] / v1);
    return ret;
  }

  bool operator==(const Vec4<T>& v1) {
    return v[0] == v1[0] && v[1] == v1[1] && v[2] == v1[2] && v[3] == v1[3];
  }

  T dot(const Vec4<T>& b) const {
    T prod = v[0] * b[0] + v[1] * b[1] + v[2] * b[2] + v[3] * b[3];
    return prod;
  }

  Vec4<T> operator-() { return Vec4<T>(-v[0], -v[1], -v[2], -v[3]); }

  /// doesn't really work for 4D. assume homogeneous coordinates instead.
  Vec4<T> cross(const Vec4<T>& b) const {
    Vec4<T> c(v[1] * b[2] - v[2] * b[1], -v[0] * b[2] + v[2] * b[0],
              v[0] * b[1] - v[1] * b[0], v[3] * b[3]);
    return c;
  }

  T norm2() const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
  }

  T norm() const { return (T)(std::sqrt(float(norm2()))); }

  void normalize() {
    T n = norm();
    if (n > 0) {
      n = 1.0f / n;
      v[0] *= n;
      v[1] *= n;
      v[2] *= n;
      v[3] *= n;
    }
  }

  void setZero() {
    v[0] = 0.f;
    v[1] = 0.f;
    v[2] = 0.f;
    v[3] = 0.f;
  }

 private:
  T v[4];
};

template <typename T>
Vec4<T> operator+(const Vec4<T>& v0, const Vec4<T>& v1) {
  return Vec4<T>(v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2], v0[3] + v1[3]);
}

template <typename T>
Vec4<T> operator-(const Vec4<T>& v0, const Vec4<T>& v1) {
  return Vec4<T>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2], v0[3] - v1[3]);
}

/// scalar multiplication
template <typename T>
inline Vec4<T> operator*(T s, const Vec4<T>& v) {
  return Vec4<T>(v[0] * s, v[1] * s, v[2] * s, v[3] * s);
}

typedef Vec4<float> Vec4f;
typedef Vec4<double> Vec4d;
typedef Vec4<unsigned> Vec4u;
typedef Vec4<unsigned char> Vec4b;
typedef Vec4<int> Vec4i;
typedef Vec4<unsigned long long> Vec4ul;

#endif
