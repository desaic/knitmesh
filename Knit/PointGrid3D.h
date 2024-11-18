#pragma
#include <vector>

#include "Array3D.h"
#include "Vec3.h"
#include "BBox.h"

class PointGrid3D {
 public:
  PointGrid3D(const BBox& box, float dx) : box_(box), dx_(dx) { Allocate(); }

  void Allocate() {
    Vec3f boxSize = box_.vmax - box_.vmin;
    Vec3u size(boxSize[0] / dx_, boxSize[1] / dx_, boxSize[2] / dx_);
    size += Vec3u(1, 1, 1);
    grid_.Allocate(size[0], size[1], size[2], {});
  }

  /// @brief Integer index for a point in the grid
  /// @param x
  /// @return
  Vec3u GetIndex(const Vec3f& x) const {
    Vec3u index;
    for (unsigned d = 0; d < 3; d++) {
      index[d] = unsigned((x[d] - box_.vmin[d]) / dx_);
    }
    return index;
  }

  bool IsIndexValid(const Vec3u& index) const {
    const Vec3u size = grid_.GetSize();
    return index[0] < size[0] && index[1] < size[1] && index[2] < size[2];
  }

  std::vector<size_t>& GetCell(const Vec3u& index) {
    return grid_(index[0], index[1], index[2]);
  }

  const std::vector<size_t>& GetCell(const Vec3u& index) const {
    return grid_(index[0], index[1], index[2]);
  }

  const Vec3u& GetSize() const { return grid_.GetSize(); }
  const BBox& GetBox() const { return box_; }
  float VoxelSize() const { return dx_; }
  const Array3D<std::vector<size_t> >& GetGrid() const { return grid_; } 
 private:
  BBox box_;
  float dx_;
  Array3D<std::vector<size_t> > grid_;
  std::vector<size_t> emptyCell_;
};

void BuildGrid(PointGrid3D& grid, const std::vector<Vec3f>& points);
std::vector<size_t> GetNeighbors(const std::vector<Vec3f>& points, int i, float radius,
                                 const PointGrid3D& grid);