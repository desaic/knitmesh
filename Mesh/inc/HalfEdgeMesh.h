#pragma once

#include <vector>

#include "PolyMesh.h"

/// @brief does not work with more than 4 billion edges or vertices.
/// change to uint64_t instead.
struct HalfEdge {
  //specialized version
  //next = 3 * (index/3) + ((index+1)%3)
  //face = index/3
  const static unsigned INVALID_IDX;
  unsigned twin = INVALID_IDX;
  unsigned vertex = 0;
  unsigned next = INVALID_IDX;
  bool hasTwin() const { return twin != INVALID_IDX; }
};

class HalfEdgeMesh {
 public:
  HalfEdgeMesh(const PolyMesh& pm);
  unsigned Face(unsigned halfEdge) const;
  unsigned Next(unsigned halfEdgeIdx)const;
  unsigned Vert(unsigned halfEdgeIdx) const;

  const std::vector<PolyFace>& F() const { return m_.F(); }
  const PolyFace& F(PolyMesh::TIndex i) const { return m_.F(i); }
  const std::vector<Vec3<PolyMesh::TValue>>& V() const { return m_.V(); }
  const Vec3<PolyMesh::TValue>& V(unsigned i) const { return m_.V(i); }

  // index of an half edge locally in its face.
  unsigned EdgeIndexInFace(unsigned halfEdge, unsigned face) const;
  unsigned HalfEdgeIndex(unsigned face, unsigned j) const;
  void BuildFaceEdgeIndex(const PolyMesh& pm);
  bool IsBadVert(unsigned v) const;
  unsigned NumFaces() const { return m_.NumFaces(); }
  unsigned NumVerts() const { return m_.NumVerts(); }
  /// @brief JUNCTION is a vertex with interleaved faces and empty spaces.
  enum VERTEX_TOPOLOGY { MANIFOLD, BOUNDARY, JUNCTION };

  VERTEX_TOPOLOGY CheckVertexManifold(unsigned v, const std::vector<unsigned>& trigs);
  /// vertex at T junction or boundary
  /// @return number of bad vertices.
  int IdentifyBadVerts();

  /// <summary>
  /// All outgoing edges from vertex v.
  /// Does not work if mesh is not a closed manifold.
  /// </summary>
  std::vector<unsigned> VertEdges(unsigned v) const;
  /// index of the first half edge pointing away from v
  std::vector<unsigned> v_;
  std::vector<HalfEdge> he;
  /// @brief first half edge index of face i.
  /// also total number of half edges before face i.
  /// the first half edge starts from the first vertex of face i.
  std::vector<unsigned> faceEdgeIndex;
  /// @brief face index for each half edge.
  std::vector<unsigned> edgeFaceIndex;
  /// index into badVertTrigs for bad vertices.
  /// good vertices stores INVALID_IDX
  std::vector<unsigned> badVertIdx;
  /// border vertices explicitly store all outgoing half edge indices.
  /// VertEdges() uses this data for bad vertices instead of looping over edges.
  std::vector<std::vector<unsigned> > badVertEdges;

  static const int MAX_LABEL = 4;
  /// <summary>
  ///  -1 for unvisited. -2 for queued.
  /// 0 1 2 3 for each orientation.
  /// 0 means the direction vector is same as the direction of
  /// edge 0 in face i. the direction vector rotates in counter-
  /// clockwise as label increases. also the offset for the starting
  /// vertex for the curve pattern.
  /// </summary>
  std::vector<int> faceLabels;

  /// vertex normals used for laying down curves.
  /// vertices are extruded in +- normal direction so that 
  /// each quad forms a hexahedron.
  std::vector<Vec3f> vn;
  const PolyMesh& m_;
};
