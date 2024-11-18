#pragma once

#include "Vec2.h"
#include "Vec3.h"

struct SurfacePoint {
  Vec3f v;
  //v= v0 + e1*l1 + e2*l2
  Vec2f l;
  Vec3f n;
  Vec2f uv;
  Vec3f tx, ty;
};
