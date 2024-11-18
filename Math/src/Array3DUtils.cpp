#include "Array3DUtils.h"
#include <deque> 

float TrilinearInterp(const Array3Df& dist, const Vec3f& x, const Vec3f& voxSize) {
  float MAX_DIST = dist.GetSize()[0] * voxSize[0] * 10.0f;
  return TrilinearInterpWithDefault(dist, x, voxSize, MAX_DIST);
}

float TrilinearInterpWithDefault(const Array3Df& dist, Vec3f local, const Vec3f& voxSize,
                                 float defaultVal) {
  local = local / voxSize;
  unsigned ix = unsigned(local[0]);
  unsigned iy = unsigned(local[1]);
  unsigned iz = unsigned(local[2]);
  Vec3u gridSize = dist.GetSize();
  if (ix >= gridSize[0] - 1 || iy >= gridSize[1] - 1 || iz >= gridSize[2] - 1) {
    return defaultVal;
  }
  float a = (ix + 1) - local[0];
  float b = (iy + 1) - local[1];
  float c = (iz + 1) - local[2];
  // v indexed by y and z
  float v00 = a * dist(ix, iy, iz) + (1 - a) * dist(ix + 1, iy, iz);
  float v10 = a * dist(ix, iy + 1, iz) + (1 - a) * dist(ix + 1, iy + 1, iz);
  float v01 = a * dist(ix, iy, iz + 1) + (1 - a) * dist(ix + 1, iy, iz + 1);
  float v11 = a * dist(ix, iy + 1, iz + 1) + (1 - a) * dist(ix + 1, iy + 1, iz + 1);
  float v0 = b * v00 + (1 - b) * v10;
  float v1 = b * v01 + (1 - b) * v11;
  float v = c * v0 + (1 - c) * v1;
  return v;
}

static bool InBound(int x, int y, int z, const Vec3u& size) {
  return x >= 0 && y >= 0 && z >= 0 && x < int(size[0]) && y < int(size[1]) && z < int(size[2]);
}

size_t LinearIdx(unsigned x, unsigned y, unsigned z, const Vec3u& size) {
  return x + y * size_t(size[0]) + z * size_t(size[0] * size[1]);
}

Vec3u GridIdx(size_t l, const Vec3u& size) {
  unsigned x = l % size[0];
  unsigned y = (l % (size[0] * size[1])) / size[0];
  unsigned z = l / (size[0] * size[1]);
  return Vec3u(x, y, z);
}

void FloodOutsideSeed(Vec3u seed, const Array3D<short>& dist, float distThresh, Array3D8u& label) {
  Vec3u size = dist.GetSize();
  size_t linearSeed = LinearIdx(seed[0], seed[1], seed[2], size);
  std::deque<size_t> q(1, linearSeed);
  const unsigned NUM_NBR = 6;
  const int nbrOffset[6][3] = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1} };
  while (!q.empty()) {
    size_t linearIdx = q.front();
    q.pop_front();
    Vec3u idx = GridIdx(linearIdx, size);
    label(idx[0], idx[1], idx[2]) = 1;
    for (unsigned ni = 0; ni < NUM_NBR; ni++) {
      int nx = int(idx[0] + nbrOffset[ni][0]);
      int ny = int(idx[1] + nbrOffset[ni][1]);
      int nz = int(idx[2] + nbrOffset[ni][2]);
      if (!InBound(nx, ny, nz, size)) {
        continue;
      }
      unsigned ux = unsigned(nx), uy = unsigned(ny), uz = unsigned(nz);
      uint8_t nbrLabel = label(ux, uy, uz);
      size_t nbrLinear = LinearIdx(ux, uy, uz, size);
      if (nbrLabel == 0) {
        uint16_t nbrDist = dist(ux, uy, uz);
        if (nbrDist >= distThresh) {
          label(ux, uy, uz) = 1;
          q.push_back(nbrLinear);
        }
      }
    }
  }
}

/// <param name="dist">distance grid. at least 1 voxel of padding is assumed.</param>
/// <param name="distThresh">stop flooding if distance is less than this</param>
/// <returns>voxel labels. 1 for outside.</returns>
Array3D8u FloodOutside(const Array3D<short>& dist, float distThresh) {
  Array3D8u label;
  Vec3u size = dist.GetSize();
  label.Allocate(size[0], size[1], size[2]);
  label.Fill(0);
  // process easy voxels first.

  // all 6 faces are labeled to be outside regardless of actual distance
  // assuming the voxel grid has padding around the mesh.
  // all threads hanging from the faces are also labeld as outside quickly.
  for (unsigned z = 0; z < size[2]; z++) {
    for (unsigned y = 0; y < size[1]; y++) {
      if (y == 0 || y == size[1] - 1 || z == 0 || z == size[2] - 1) {
        for (unsigned x = 0; x < size[0]; x++) {
          label(x, y, z) = 1;
        }
        continue;
      }
      label(0, y, z) = 1;
      label(size[0] - 1, y, z) = 1;
      for (unsigned x = 1; x < size[0] - 1; x++) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
      for (unsigned x = size[0] - 1; x > 0; x--) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
    }
  }

  for (unsigned x = 0; x < size[0]; x++) {
    for (unsigned y = 0; y < size[1]; y++) {
      for (unsigned z = 1; z < size[2] - 1; z++) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
      for (unsigned z = size[2] - 1; z > 0; z--) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
    }
  }

  for (unsigned z = 0; z < size[2]; z++) {
    for (unsigned x = 0; x < size[0]; x++) {
      for (unsigned y = 1; y < size[1] - 1; y++) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
      for (unsigned y = size[1] - 1; y > 0; y--) {
        short d = dist(x, y, z);
        if (d <= distThresh) {
          break;
        }
        label(x, y, z) = 1;
      }
    }
  }

  const unsigned NUM_NBR = 6;
  const int nbrOffset[6][3] = { {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  for (unsigned z = 1; z < size[2] - 1; z++) {
    for (unsigned y = 1; y < size[1] - 1; y++) {
      for (unsigned x = 1; x < size[0] - 1; x++) {
        for (unsigned ni = 0; ni < NUM_NBR; ni++) {
          unsigned ux = unsigned(x + nbrOffset[ni][0]);
          unsigned uy = unsigned(y + nbrOffset[ni][1]);
          unsigned uz = unsigned(z + nbrOffset[ni][2]);

          if (!label(x, y, z)) {
            continue;
          }
          uint8_t nbrLabel = label(ux, uy, uz);
          if (nbrLabel == 0) {
            uint16_t nbrDist = dist(ux, uy, uz);
            if (nbrDist >= distThresh) {
              FloodOutsideSeed(Vec3u(x, y, z), dist, distThresh, label);
            }
          }
        }
      }
    }
  }
  return label;
}
