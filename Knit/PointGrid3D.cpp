#include "PointGrid3D.h"

/// @brief
/// @param grid
/// @param points
/// @param origin
/// @param h grid cell size.
void BuildGrid(PointGrid3D& grid, const std::vector<Vec3f>& points) {
  for (size_t i = 0; i < points.size(); i++) {
    Vec3u index = grid.GetIndex(points[i]);
    if (grid.IsIndexValid(index)) {
      grid.GetCell(index).push_back(i);
    }
  }
}

/// @brief get list of nearby points within radius relative to center of
/// bounding box for struct i.
std::vector<size_t> GetNeighbors(const std::vector<Vec3f>& points,int i, float radius,
                                 const PointGrid3D& grid) {
  std::vector<size_t> neighbors;
  Vec3f point = points[i];
  Vec3u index = grid.GetIndex(point);
  if (!grid.IsIndexValid(index)) {
    return neighbors;
  }
  const int MAX_SEARCH_RADIUS = 4;
  int searchRadius = int(radius / grid.VoxelSize()) + 1;
  searchRadius = std::min(MAX_SEARCH_RADIUS, searchRadius);
  for (int dz = -searchRadius; dz <= searchRadius; dz++) {
    for (int dy = -searchRadius; dy <= searchRadius; dy++) {
      for (int dx = -searchRadius; dx <= searchRadius; dx++) {
        Vec3u neighborIndex(index[0] + dx, index[1] + dy, index[2] + dz);
        if (!grid.IsIndexValid(neighborIndex)) {
          continue;
        }
        const std::vector<size_t> nbrIndex = grid.GetCell(neighborIndex);
        for (auto j : nbrIndex) {
          if (j == i) {
            continue;
          }
          const Vec3f& neighbor = points[j];
          Vec3f Pij = neighbor - point;
          if (Pij.norm2() < radius * radius) {
            neighbors.push_back(j);
          }
        }
      }
    }
  }
  return neighbors;
}
