#include "Matrix3f.h"

#include <math.h>

#include <cmath>
#include <cstring>

#include "Vec3.h"

#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a) > (b) ? (a) : (b))

Matrix3f::Matrix3f() {
  for (int i = 0; i < 9; ++i) {
    m_elements[i] = 0.0f;
  }
}

Matrix3f::Matrix3f(float m00, float m01, float m02, float m10, float m11,
                   float m12, float m20, float m21, float m22) {
  m_elements[0] = m00;
  m_elements[1] = m10;
  m_elements[2] = m20;

  m_elements[3] = m01;
  m_elements[4] = m11;
  m_elements[5] = m21;

  m_elements[6] = m02;
  m_elements[7] = m12;
  m_elements[8] = m22;
}

Matrix3f::Matrix3f(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2,
                   bool setColumns) {
  if (setColumns) {
    setCol(0, v0);
    setCol(1, v1);
    setCol(2, v2);
  } else {
    setRow(0, v0);
    setRow(1, v1);
    setRow(2, v2);
  }
}

Matrix3f::Matrix3f(const Matrix3f& rm) {
  memcpy(m_elements, rm.m_elements, 9 * sizeof(float));
}

Matrix3f& Matrix3f::operator=(const Matrix3f& rm) {
  if (this != &rm) {
    memcpy(m_elements, rm.m_elements, 9 * sizeof(float));
  }
  return *this;
}

Matrix3f& Matrix3f::operator+=(const Matrix3f& rm) {
  for (size_t i = 0; i < 9; i++) {
    m_elements[i] += rm.m_elements[i];
  }
  return *this;
}

Matrix3f Matrix3f::operator-() const { 
  Matrix3f m;
  for (unsigned i = 0; i < 9; i++) {
    m.m_elements[i] = -m_elements[i];
  }
  return m; 
}

const float& Matrix3f::operator()(int i, int j) const {
  return m_elements[j * 3 + i];
}

float& Matrix3f::operator()(int i, int j) { return m_elements[j * 3 + i]; }

Vec3f Matrix3f::getRow(int i) const {
  return Vec3f(m_elements[i], m_elements[i + 3], m_elements[i + 6]);
}

void Matrix3f::setRow(int i, const Vec3f& v) {
  m_elements[i] = v[0];
  m_elements[i + 3] = v[1];
  m_elements[i + 6] = v[2];
}

Vec3f Matrix3f::getCol(int j) const {
  int colStart = 3 * j;

  return Vec3f(m_elements[colStart], m_elements[colStart + 1],
               m_elements[colStart + 2]);
}

void Matrix3f::setCol(int j, const Vec3f& v) {
  int colStart = 3 * j;

  m_elements[colStart] = v[0];
  m_elements[colStart + 1] = v[1];
  m_elements[colStart + 2] = v[2];
}

float Matrix3f::determinant() const {
  return Matrix3f::determinant3x3(m_elements[0], m_elements[3], m_elements[6],
                                  m_elements[1], m_elements[4], m_elements[7],
                                  m_elements[2], m_elements[5], m_elements[8]);
}

float det2x2(float m00, float m01, float m10, float m11) {
  return m00 * m11 - m01 * m10;
}

Matrix3f Matrix3f::inverse(bool* pbIsSingular, float epsilon) const {
  float m00 = m_elements[0];
  float m10 = m_elements[1];
  float m20 = m_elements[2];

  float m01 = m_elements[3];
  float m11 = m_elements[4];
  float m21 = m_elements[5];

  float m02 = m_elements[6];
  float m12 = m_elements[7];
  float m22 = m_elements[8];

  float cofactor00 = det2x2(m11, m12, m21, m22);
  float cofactor01 = -det2x2(m10, m12, m20, m22);
  float cofactor02 = det2x2(m10, m11, m20, m21);

  float cofactor10 = -det2x2(m01, m02, m21, m22);
  float cofactor11 = det2x2(m00, m02, m20, m22);
  float cofactor12 = -det2x2(m00, m01, m20, m21);

  float cofactor20 = det2x2(m01, m02, m11, m12);
  float cofactor21 = -det2x2(m00, m02, m10, m12);
  float cofactor22 = det2x2(m00, m01, m10, m11);

  float determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02;

  bool isSingular = (fabs(determinant) < epsilon);
  if (isSingular) {
    if (pbIsSingular != NULL) {
      *pbIsSingular = true;
    }
    return Matrix3f();
  } else {
    if (pbIsSingular != NULL) {
      *pbIsSingular = false;
    }

    float reciprocalDeterminant = 1.0f / determinant;

    return Matrix3f(
        cofactor00 * reciprocalDeterminant, cofactor10 * reciprocalDeterminant,
        cofactor20 * reciprocalDeterminant, cofactor01 * reciprocalDeterminant,
        cofactor11 * reciprocalDeterminant, cofactor21 * reciprocalDeterminant,
        cofactor02 * reciprocalDeterminant, cofactor12 * reciprocalDeterminant,
        cofactor22 * reciprocalDeterminant);
  }
}

float Matrix3f::squaredNorm() const {
  float sum = 0;
  for (size_t i = 0; i < 9; i++) {
    sum += m_elements[i] *  m_elements[i];
  }
  return sum;
}

float Matrix3f::trace() const {
  return (*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2);
}

void Matrix3f::transpose() {
  float temp;

  for (int i = 0; i < 2; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      temp = (*this)(i, j);
      (*this)(i, j) = (*this)(j, i);
      (*this)(j, i) = temp;
    }
  }
}

Matrix3f Matrix3f::transposed() const {
  Matrix3f out;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out(j, i) = (*this)(i, j);
    }
  }

  return out;
}

Matrix3f::operator float*() { return m_elements; }

// static
float Matrix3f::determinant3x3(float m00, float m01, float m02, float m10,
                               float m11, float m12, float m20, float m21,
                               float m22) {
  return (m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) +
          m02 * (m10 * m21 - m11 * m20));
}

Matrix3f Matrix3f::Diag(float m00, float m11, float m22) {
  Matrix3f m;
  for (int i = 1; i < 8; ++i) {
    m.m_elements[i] = 0;
  }
  m.m_elements[0] = m00;
  m.m_elements[4] = m11;
  m.m_elements[8] = m22;
  return m;
}

// static
Matrix3f Matrix3f::ones() {
  Matrix3f m;
  for (int i = 0; i < 9; ++i) {
    m.m_elements[i] = 1;
  }

  return m;
}

Matrix3f Matrix3f::Zero() {
  Matrix3f m;
  for (int i = 0; i < 9; ++i) {
    m.m_elements[i] = 0;
  }
  return m;
}

// static
Matrix3f Matrix3f::identity() {
  Matrix3f m;

  m(0, 0) = 1;
  m(1, 1) = 1;
  m(2, 2) = 1;

  return m;
}

// static
Matrix3f Matrix3f::rotateX(float radians) {
  float c = cos(radians);
  float s = sin(radians);

  return Matrix3f(1, 0, 0, 0, c, -s, 0, s, c);
}

// static
Matrix3f Matrix3f::rotateY(float radians) {
  float c = cos(radians);
  float s = sin(radians);

  return Matrix3f(c, 0, s, 0, 1, 0, -s, 0, c);
}

// static
Matrix3f Matrix3f::rotateZ(float radians) {
  float c = cos(radians);
  float s = sin(radians);

  return Matrix3f(c, -s, 0, s, c, 0, 0, 0, 1);
}

// static
Matrix3f Matrix3f::scaling(float sx, float sy, float sz) {
  return Matrix3f(sx, 0, 0, 0, sy, 0, 0, 0, sz);
}

// static
Matrix3f Matrix3f::uniformScaling(float s) {
  return Matrix3f(s, 0, 0, 0, s, 0, 0, 0, s);
}

// static
Matrix3f Matrix3f::rotation(const Vec3f& axis, float radians) {
  Vec3f normalizedAxis = axis;
  normalizedAxis.normalize();

  float cosTheta = cos(radians);
  float sinTheta = sin(radians);

  float x = normalizedAxis[0];
  float y = normalizedAxis[1];
  float z = normalizedAxis[2];

  return Matrix3f(x * x * (1.0f - cosTheta) + cosTheta,
                  y * x * (1.0f - cosTheta) - z * sinTheta,
                  z * x * (1.0f - cosTheta) + y * sinTheta,
                  x * y * (1.0f - cosTheta) + z * sinTheta,
                  y * y * (1.0f - cosTheta) + cosTheta,
                  z * y * (1.0f - cosTheta) - x * sinTheta,
                  x * z * (1.0f - cosTheta) - y * sinTheta,
                  y * z * (1.0f - cosTheta) + x * sinTheta,
                  z * z * (1.0f - cosTheta) + cosTheta);
}

// static
Matrix3f Matrix3f::rotation(float qx, float qy, float qz, float qw) {
  float norm2 = qx * qx + qy * qy + qz * qz + qw * qw;
  float invNorm2 = 1.0f / norm2;

  float xx = invNorm2 * qx * qx;
  float yy = invNorm2 * qy * qy;
  float zz = invNorm2 * qz * qz;

  float xy = invNorm2 * qx * qy;
  float zw = invNorm2 * qz * qw;

  float xz = invNorm2 * qx * qz;
  float yw = invNorm2 * qy * qw;

  float yz = invNorm2 * qy * qz;
  float xw = invNorm2 * qx * qw;

  return Matrix3f(1.0f - 2.0f * (yy + zz), 2.0f * (xy - zw), 2.0f * (xz + yw),
                  2.0f * (xy + zw), 1.0f - 2.0f * (xx + zz), 2.0f * (yz - xw),
                  2.0f * (xz - yw), 2.0f * (yz + xw), 1.0f - 2.0f * (xx + yy));
}

//////////////////////////////////////////////////////////////////////////
// Operators
//////////////////////////////////////////////////////////////////////////

Vec3f operator*(const Matrix3f& m, const Vec3f& v) {
  Vec3f output(0, 0, 0);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      output[i] += m(i, j) * v[j];
    }
  }

  return output;
}

Matrix3f operator*(const Matrix3f& x, const Matrix3f& y) {
  Matrix3f product;  // zeroes

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        product(i, k) += x(i, j) * y(j, k);
      }
    }
  }

  return product;
}

Matrix3f operator*(const Matrix3f& x, float scalar) {
  Matrix3f product;  // zeroes

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      product(i, j) = x(i, j) * scalar;
    }
  }

  return product;
}

Matrix3f operator*(float scalar, const Matrix3f& x) {
  Matrix3f product;  // zeroes

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      product(i, j) += x(i, j) * scalar;
    }
  }
  return product;
}

Matrix3f operator+(const Matrix3f& x, const Matrix3f& y) {
  Matrix3f sum;  // zeroes

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      sum(i, j) = x(i, j) + y(i, j);
    }
  }

  return sum;
}

Matrix3f operator-(const Matrix3f& x, const Matrix3f& y) {
  Matrix3f diff;  // zeroes

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      diff(i, j) = x(i, j) - y(i, j);
    }
  }

  return diff;
}

// Slightly modified version of  Stan Melax's code for 3x3 matrix
// diagonalization (Thanks Stan!) source:
// http://www.melax.com/diag.html?attredirects=0
Vec3f Matrix3f::getLowestEigenvector() {
  // A must be a symmetric matrix.
  // returns Q and D such that
  // Diagonal matrix D = QT * A * Q; and A = Q*D*QT

  float A[3][3];
  float Q[3][3];
  float D[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A[i][j] = m_elements[i + 3 * j];
    }
  }

  const int maxsteps = 24;  // certainly wont need that many.
  int k0, k1, k2;
  float o[3], m[3];
  float q[4] = {0.0, 0.0, 0.0, 1.0};
  float jr[4];
  float sqw, sqx, sqy, sqz;
  float tmp1, tmp2, mq;
  float AQ[3][3];
  float thet, sgn, t, c;
  for (int i = 0; i < maxsteps; ++i) {
    // quat to matrix
    sqx = q[0] * q[0];
    sqy = q[1] * q[1];
    sqz = q[2] * q[2];
    sqw = q[3] * q[3];
    Q[0][0] = (sqx - sqy - sqz + sqw);
    Q[1][1] = (-sqx + sqy - sqz + sqw);
    Q[2][2] = (-sqx - sqy + sqz + sqw);
    tmp1 = q[0] * q[1];
    tmp2 = q[2] * q[3];
    Q[1][0] = 2.0 * (tmp1 + tmp2);
    Q[0][1] = 2.0 * (tmp1 - tmp2);
    tmp1 = q[0] * q[2];
    tmp2 = q[1] * q[3];
    Q[2][0] = 2.0 * (tmp1 - tmp2);
    Q[0][2] = 2.0 * (tmp1 + tmp2);
    tmp1 = q[1] * q[2];
    tmp2 = q[0] * q[3];
    Q[2][1] = 2.0 * (tmp1 + tmp2);
    Q[1][2] = 2.0 * (tmp1 - tmp2);

    // AQ = A * Q
    AQ[0][0] = Q[0][0] * A[0][0] + Q[1][0] * A[0][1] + Q[2][0] * A[0][2];
    AQ[0][1] = Q[0][1] * A[0][0] + Q[1][1] * A[0][1] + Q[2][1] * A[0][2];
    AQ[0][2] = Q[0][2] * A[0][0] + Q[1][2] * A[0][1] + Q[2][2] * A[0][2];
    AQ[1][0] = Q[0][0] * A[0][1] + Q[1][0] * A[1][1] + Q[2][0] * A[1][2];
    AQ[1][1] = Q[0][1] * A[0][1] + Q[1][1] * A[1][1] + Q[2][1] * A[1][2];
    AQ[1][2] = Q[0][2] * A[0][1] + Q[1][2] * A[1][1] + Q[2][2] * A[1][2];
    AQ[2][0] = Q[0][0] * A[0][2] + Q[1][0] * A[1][2] + Q[2][0] * A[2][2];
    AQ[2][1] = Q[0][1] * A[0][2] + Q[1][1] * A[1][2] + Q[2][1] * A[2][2];
    AQ[2][2] = Q[0][2] * A[0][2] + Q[1][2] * A[1][2] + Q[2][2] * A[2][2];
    // D = Qt * AQ
    D[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
    D[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
    D[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
    D[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
    D[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
    D[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
    D[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
    D[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
    D[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];
    o[0] = D[1][2];
    o[1] = D[0][2];
    o[2] = D[0][1];
    m[0] = fabs(o[0]);
    m[1] = fabs(o[1]);
    m[2] = fabs(o[2]);

    k0 = (m[0] > m[1] && m[0] > m[2]) ? 0
         : (m[1] > m[2])              ? 1
                         : 2;  // index of largest element of offdiag
    k1 = (k0 + 1) % 3;
    k2 = (k0 + 2) % 3;
    if (o[k0] == 0.0) {
      break;  // diagonal already
    }
    thet = (D[k2][k2] - D[k1][k1]) / (2.0 * o[k0]);
    sgn = (thet > 0.0) ? 1.0 : -1.0;
    thet *= sgn;  // make it positive
    t = sgn / (thet + ((thet < 1.E6) ? sqrt(thet * thet + 1.0)
                                     : thet));  // sign(T)/(|T|+sqrt(T^2+1))
    c = 1.0 / sqrt(t * t + 1.0);                // c= 1/(t^2+1) , t=s/c
    if (c == 1.0) {
      break;  // no room for improvement - reached machine precision.
    }
    jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
    jr[k0] =
        sgn *
        sqrt((1.0 - c) /
             2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
    jr[k0] *=
        -1.0;  // since our quat-to-matrix convention was for v*M instead of M*v
    jr[3] = sqrt(1.0f - jr[k0] * jr[k0]);
    if (jr[3] == 1.0) {
      break;  // reached limits of floating point precision
    }
    q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
    q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
    q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
    q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
    mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0] /= mq;
    q[1] /= mq;
    q[2] /= mq;
    q[3] /= mq;
  }

  // find correct eigenvector
  int min_ind = 0;
  for (int i = 0; i < 3; i++) {
    if (D[i][i] < D[min_ind][min_ind]) {
      min_ind = i;
    }
  }

  return Vec3f(Q[0][min_ind], Q[1][min_ind], Q[2][min_ind]);
}

Matrix3f OuterProd(const Vec3f& v1, const Vec3f& v2) {
  // clang-format off
  Matrix3f mat(v1[0] * v2[0], v1[0] * v2[1], v1[0] * v2[2], 
               v1[1] * v2[0], v1[1] * v2[1], v1[1] * v2[2], 
               v1[2] * v2[0], v1[2] * v2[1], v1[2] * v2[2]);
  // clang-format on
  return mat;
}