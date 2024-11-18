#pragma once
#include "Vec3.h"
#include "Array3D.h"

class TrigMesh;

struct GridCell {
  static const unsigned NUM_PT = 8;
  // voxel corner coordinates
  //   4----5
  // 6----7
  //
  //   0----1
  // 2----3
   Vec3f p[NUM_PT];
   float val[NUM_PT];	//value of the function at this grid corner
   GridCell() {
     for (unsigned i = 0; i < NUM_PT; i++) {
       val[i] = 1e2f;
     }
   }
};

///\param mesh output a mini mesh for cell
void MarchCube(GridCell &cell, float level, TrigMesh*mesh);

template<typename T>
void MarchOneCube(unsigned x, unsigned y, unsigned z, const Array3D<T>& grid,
                  float level, TrigMesh* surf, float h);

template <typename T>
void MarchingCubes(const Array3D<T>& grid, float level, TrigMesh* surf,
                   float h);
