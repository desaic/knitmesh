#pragma once
#include "Vec3.h"
#include "Matrix3f.h"
#include <vector>

void rotateVerts(const std::vector<float>& verts, std::vector<float>& rotated,
  float rx, float ry, float rz);

void transformVerts(const std::vector<float> & verts, std::vector<float> & dst, 
  const Matrix3f mat);

Matrix3f RotationMatrix(float rx, float ry, float rz);