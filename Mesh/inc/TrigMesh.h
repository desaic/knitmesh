#pragma once
#ifndef TRIG_MESH_H
#define TRIG_MESH_H
#include <string>
#include <vector>

#include "Vec2.h"
#include "Vec3.h"

struct Triangle {
  Vec3f v[3];
  const Vec3f& operator[](unsigned i) const { return v[i]; }
  Vec3f& operator[](unsigned i) { return v[i]; }
};

class TrigMesh {
 public:
  TrigMesh();
  using TValue = float;
  using TIndex = unsigned;

  std::vector<TValue> v;
  std::vector<TIndex> t;

  // these are all needed just for SDF
  /// triangle edges
  std::vector<TIndex> te;

  /// normal of each triangle
  std::vector<TValue> nt;

  /// averaged edge normal
  std::vector<Vec3f> ne;

  /// angle weighted vertex normal
  std::vector<Vec3f> nv;

  // texture coordinates
  std::vector<Vec2f> uv;

  // delete and free all data.
  void Clear();
  size_t GetNumTrigs() const { return t.size() / 3; }
  size_t GetNumVerts() const { return v.size() / 3; }
  int LoadStl(const std::string& filename);
  int LoadObj(const std::string& filename);
  int LoadStep(const std::string& filename);
  std::vector<TValue>& V() { return v; }
  std::vector<TIndex>& T() { return t; }
  const std::vector<TValue>& V()const { return v; }
  const std::vector<TIndex>& T()const { return t; }
  // can get face, edge or vertex normal depending on bary centric coord.
  Vec3f GetNormal(unsigned tIdx, const Vec3f& bary)const;
  Vec3f GetTrigNormal(size_t tIdx)const;
  Triangle GetTriangleVerts(size_t tIdx) const;
  Vec2f GetTriangleUV(unsigned tIdx, unsigned j) const;
  // compute triangle area.
  float GetArea(unsigned tIdx) const;
  Vec3f GetRandomSample(const unsigned tIdx, const float r1,
                        const float r2) const;
  // also returns trig normal and a tangent vector
  void GetRandomSample(const unsigned tIdx, const float r1, const float r2,
                       Vec3f& pos, Vec3f& normal, Vec3f& tangent) const;

  void ComputeTrigNormals();
  void ComputeVertNormals();
  /// experimental feature. do not use.
  void ComputePseudoNormals();
  void ComputeTrigEdges();

  unsigned GetIndex(unsigned tIdx, unsigned jIdx) const;
  Vec3f GetVertex(unsigned vIdx) const;
  Vec3f GetTriangleVertex(unsigned tIdx, unsigned jIdx) const;
  void AddVert(float x, float y, float z);
  void scale(float s);
  void translate(float dx, float dy, float dz);
  void append(const TrigMesh& m);

  void SaveObj(const std::string& filename);
  // save out a point cloud as obj.
  static void SaveObj(const std::string& filename,
                      const std::vector<Vec3f>& pts);
  int SaveStlTxt(const std::string& filename);
  int SaveStl(const std::string& filename);
};

TrigMesh SubSet(const TrigMesh& m, const std::vector<size_t>& trigs);
TrigMesh MakeCube(const Vec3f& mn, const Vec3f& mx);
TrigMesh MakePlane(const Vec3f& mn, const Vec3f& mx, const Vec3f& n);

#endif