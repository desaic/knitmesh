#pragma once
#ifndef POLY_MESH_H
#define POLY_MESH_H
#include <string>
#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "FixedArray.h"

struct PolyFace {
  using TIndex = unsigned;
  TIndex n = 3;
  FixedArray<TIndex> v;
  PolyFace(const std::vector<unsigned> indices) { 
    n = indices.size();
    v.Allocate(n);
    for (unsigned i = 0; i < n; i++) {
      v[i] = indices[i];
    }
  }
  TIndex& operator[](unsigned vi) { return v[vi]; }
  TIndex operator[](unsigned vi) const { return v[vi]; }
  static const int MAX_N = 5;
};

class PolyMesh {
 public:
  PolyMesh();
  using TValue = float;
  using TIndex = PolyFace::TIndex;

  // delete and free all data.
  void Clear();
  size_t NumFaces() const { return f.size(); }
  size_t NumVerts() const { return v.size(); }
  void Scale(float s);
  void Translate(float dx, float dy, float dz);
  void Append(const PolyMesh& m);

  int LoadObj(const std::string& filename);
  void SaveObj(const std::string& filename);
  
  std::vector<PolyFace>& F() { return f; }
  const std::vector<PolyFace>& F() const { return f; }
  PolyFace & F(TIndex i) { return f[i]; }
  const PolyFace& F(TIndex i) const { return f[i]; }
  
  std::vector<Vec3<TValue> >& V() { return v; }
  const std::vector<Vec3<TValue> >& V() const { return v; }
  Vec3<TValue> & V(unsigned i) { return v[i]; }
  const Vec3<TValue>& V(unsigned i) const { return v[i]; }

 private:
  std::vector<Vec3<TValue> > v;
  std::vector<PolyFace> f;
};

PolyMesh SubSet(const PolyMesh& m, const std::vector<size_t>& faces);
PolyMesh MakeQuadCube(const Vec3f& mn, const Vec3f& mx);

#endif