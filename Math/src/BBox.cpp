#include "BBox.h"
#include <vector>
#include <algorithm>

void UpdateBBox(const std::vector<float>& verts, BBox& box)
{
  const unsigned dim = 3;
  for (size_t i = 0; i < verts.size(); i+=3) {
    for (unsigned j = 0; j < dim; j++) {
      box.vmin[j] = std::min(box.vmin[j], verts[i + j]);
      box.vmax[j] = std::max(box.vmax[j], verts[i + j]);
    }
  }
}

void ComputeBBox(const std::vector<float>& verts, BBox& box)
{  
  if (verts.size() == 0) {
    return;
  }
  const unsigned dim = 3;
  for (size_t i = 0; i < 3; i++) {
    box.vmin[i] = verts[i];
    box.vmax[i] = verts[i];
  }
  box._init = true;
  UpdateBBox(verts, box);
}

void ComputeBBox(const std::vector<Vec3f>& verts, BBox& box) {
  if (verts.size() == 0) {
    return;
  }
  box.vmin = verts[0];
  box.vmax = verts[0];
  const unsigned dim = 3;
  for (size_t i = 0; i < verts.size(); i++) {
    for (unsigned j = 0; j < dim; j++) {
      box.vmin[j] = std::min(box.vmin[j], verts[i][j]);
      box.vmax[j] = std::max(box.vmax[j], verts[i][j]);
    }
  }
  box._init = true;
}

void ComputeBBox(Vec3f v0, Vec3f v1, Vec3f v2, BBox& box) {
  const unsigned dim = 3;
  box.vmin = v0;
  box.vmax = v0;
  for (unsigned j = 0; j < dim; j++) {
    box.vmin[j] = std::min(box.vmin[j], v1[j]);
    box.vmax[j] = std::max(box.vmax[j], v1[j]);
    box.vmin[j] = std::min(box.vmin[j], v2[j]);
    box.vmax[j] = std::max(box.vmax[j], v2[j]);
  }
  box._init = true;
}

void BBox::Merge(const BBox& b) {
  if (!b._init) {
    return;
  }
  if (!_init) {
    vmin = b.vmin;
    vmax = b.vmax;
    _init = true;
  } else {
    for (size_t i = 0; i < 3; i++) {
      vmin[i] = std::min(b.vmin[i], vmin[i]);
      vmax[i] = std::max(b.vmax[i], vmax[i]);
    }
  }
}
