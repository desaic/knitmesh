#pragma once
#include <vector>

#include "Vec3.h"

struct BBox {
  /// vertices at min and max coordinates
  Vec3f vmin, vmax;
  BBox() : vmin(0.0f, 0.0f, 0.0f), vmax(0.0f, 0.0f, 0.0f) {}
  BBox(Vec3f a, Vec3f b) : vmin(a), vmax(b) {}
  void Merge(const BBox& b);
  bool _init = false;
};

/// oriented bounding box
struct OBBox {
  Vec3f origin;
  Vec3f axes[3];
  OBBox() : origin(0) {
    axes[0][0] = 1;
    axes[1][1] = 1;
    axes[2][2] = 1;
  }
};

struct IntBox {
  /// vertices at min and max coordinates
  Vec3i min, max;
  IntBox() : min(0, 0, 0), max(0, 0, 0) {}
};

void UpdateBBox(const std::vector<float>& verts, BBox& box);

void ComputeBBox(const std::vector<float>& verts, BBox& box);
void ComputeBBox(Vec3f v0, Vec3f v1, Vec3f v2, BBox& box);
void ComputeBBox(const std::vector<Vec3f>& verts, BBox& box);