#ifndef MATRIX3F_H
#define MATRIX3F_H

#include "Vec3.h"

// 3x3 Matrix, stored in column major order (OpenGL style)
class Matrix3f {
 public:
  // default to 0.
  Matrix3f();
  Matrix3f(float m00, float m01, float m02, float m10, float m11, float m12,
           float m20, float m21, float m22);

  // setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2]
  // otherwise, sets the rows
  Matrix3f(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2,
           bool setColumns = true);

  Matrix3f(const Matrix3f& rm);             // copy constructor
  Matrix3f& operator=(const Matrix3f& rm);  // assignment operator
  Matrix3f& operator+=(const Matrix3f& rm);
  //negation operator
  Matrix3f operator-() const;  
  // no destructor necessary

  const float& operator()(int i, int j) const;
  float& operator()(int i, int j);

  Vec3f getRow(int i) const;
  void setRow(int i, const Vec3f& v);

  Vec3f getCol(int j) const;
  void setCol(int j, const Vec3f& v);

  float determinant() const;
  Matrix3f inverse(bool* pbIsSingular = NULL,
                   float epsilon = 0.f) const;  // TODO: invert in place as well
  float squaredNorm() const;

  float trace() const;

  void transpose();
  Matrix3f transposed() const;

  Vec3f getLowestEigenvector();

  // ---- Utility ----
  operator float*();  // automatic type conversion for GL

  static float determinant3x3(float m00, float m01, float m02, float m10,
                              float m11, float m12, float m20, float m21,
                              float m22);

  static Matrix3f Diag(float m00,float m11, float m22);
  static Matrix3f ones();
  static Matrix3f Zero();
  static Matrix3f identity();
  static Matrix3f rotateX(float radians);
  static Matrix3f rotateY(float radians);
  static Matrix3f rotateZ(float radians);
  static Matrix3f scaling(float sx, float sy, float sz);
  static Matrix3f uniformScaling(float s);
  static Matrix3f rotation(const Vec3f& axis, float radians);

  // Returns the rotation matrix represented by a unit quaternion
  // if q is not normalized, it it normalized first
  static Matrix3f rotation(float qx, float qy, float qz, float qw);

 private:
  float m_elements[9];
};

// Matrix-Vector multiplication
// 3x3 * 3x1 ==> 3x1
Vec3f operator*(const Matrix3f& m, const Vec3f& v);
Matrix3f operator*(float scalar, const Matrix3f& x);

// Matrix-Matrix multiplication
Matrix3f operator*(const Matrix3f& x, const Matrix3f& y);

Matrix3f operator*(const Matrix3f& x, const float scalar);

Matrix3f operator+(const Matrix3f& x, const Matrix3f& y);

Matrix3f operator-(const Matrix3f& x, const Matrix3f& y);

Matrix3f OuterProd(const Vec3f& v1, const Vec3f& v2);

#endif  // MATRIX3F_H
