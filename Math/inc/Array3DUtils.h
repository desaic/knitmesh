#include "Array3D.h"
float TrilinearInterp(const Array3Df& dist, const Vec3f& x, const Vec3f& voxSize);
float TrilinearInterpWithDefault(const Array3Df& dist, Vec3f x, const Vec3f& voxSize,
                                 float defaultVal);
Array3D8u FloodOutside(const Array3D<short>& dist, float distThresh);