#include "HalfEdgeMesh.h"
#include "Edge.h"
#include <map>
#include <iostream>

const unsigned HalfEdge::INVALID_IDX = 0xffffffff;

void HalfEdgeMesh::BuildFaceEdgeIndex(const PolyMesh& pm) {
  faceEdgeIndex.resize(pm.NumFaces());
  size_t edgeCount = 0;
  for (size_t i = 0; i < pm.NumFaces(); i++) {
    faceEdgeIndex[i] = edgeCount;
    edgeCount += pm.F(i).n;
  }
  edgeFaceIndex.resize(edgeCount);
  for (size_t f = 0; f < pm.NumFaces(); f++) {
    for (unsigned e = 0; e < pm.F(f).n; e++) {
      edgeFaceIndex[faceEdgeIndex[f] + e] = f;
    }
  }
}

unsigned HalfEdgeMesh::EdgeIndexInFace(unsigned halfEdge, unsigned face) const {
  return halfEdge - faceEdgeIndex[face];
}

unsigned HalfEdgeMesh::HalfEdgeIndex(unsigned face, unsigned j) const {
  return faceEdgeIndex[face] + j;
}

HalfEdgeMesh::HalfEdgeMesh(const PolyMesh& pm) : m_(pm) {
  unsigned numFaces = pm.NumFaces();
  unsigned numVerts = pm.NumVerts();
  v_.resize(numVerts);
  std::fill(v_.begin(), v_.end(), HalfEdge::INVALID_IDX);  
  BuildFaceEdgeIndex(pm);
  he.resize(edgeFaceIndex.size());
  //map from an edge to one half edge 
  std::map<Edge, unsigned> edges;
  for (unsigned i = 0; i < numFaces; i++) {
    const auto& f = pm.F(i);
    for (unsigned vj = 0; vj < f.n; vj++) {
      unsigned srcVert = pm.F(i)[vj];
      Edge e(srcVert, pm.F(i)[(vj + 1) % f.n]);
      auto it = edges.find(e);
      unsigned halfEdgeIdx = faceEdgeIndex[i] + vj;
      unsigned twinIdx = HalfEdge::INVALID_IDX;
      if (it == edges.end()) {
        edges[e] = halfEdgeIdx;
      } else {
        unsigned idx = it->second;
        if (idx == halfEdgeIdx) {
          //shouldn't be possible;
          twinIdx = he[idx].twin;
        } else {
          twinIdx = idx;
        }
        if (twinIdx < he.size()) {
          he[twinIdx].twin = halfEdgeIdx;
        }
      }
      HalfEdge& h = he[halfEdgeIdx];
      if (twinIdx != HalfEdge::INVALID_IDX) {
        h.twin = twinIdx;
      }
      h.vertex = pm.F(i)[vj];
      if (v_[srcVert] == HalfEdge::INVALID_IDX) {
        v_[srcVert] = halfEdgeIdx;
      }
    }
  }
  for (unsigned i = 0; i < numFaces; i++) {
    const auto& f = pm.F(i);
    for (unsigned vj = 0; vj < f.n; vj++) {
      unsigned e0 = faceEdgeIndex[i] + vj;
      unsigned e1 = faceEdgeIndex[i] + ((vj + 1) % f.n);
      he[e0].next = e1;
    }
  }
  IdentifyBadVerts();
}

bool HalfEdgeMesh::IsBadVert(unsigned v) const { return badVertIdx[v] != HalfEdge::INVALID_IDX; }

bool AddToSet(unsigned f, std::vector<bool>& hit, const std::vector<unsigned>& faces) {
  for (unsigned i = 0; i < faces.size(); i++) {
    if (faces[i] == f) {
      if (hit[i]) {
        return false;
      }
      hit[i] = true;
      return true;
    }
  }
  return false;
}

HalfEdgeMesh::VERTEX_TOPOLOGY HalfEdgeMesh::CheckVertexManifold(
    unsigned v, const std::vector<unsigned>& faces) {
  //using unordered_map is ~2x slower
  std::vector<bool> hit(faces.size(), false);
  unsigned edge = v_[v];  
  if (!AddToSet(edgeFaceIndex[edge], hit, faces)) {
    return BOUNDARY;
  }
  while(true) {
    unsigned f = edgeFaceIndex[edge];
    unsigned n = m_.F(f).n;
    unsigned ei = EdgeIndexInFace(edge, f);
    unsigned nextEdge = faceEdgeIndex[f] + (ei + n - 1) % n;
    const auto halfEdge = he[nextEdge];
    if (!halfEdge.hasTwin()) {
      return BOUNDARY;
    }
    unsigned twin = halfEdge.twin;
    edge = twin;
    if (edge == v_[v]) {
      //looped back
      break;
    }
    if (!AddToSet(edgeFaceIndex[edge], hit, faces)) {
      return JUNCTION;
    }
  }
  for (auto b:hit) {
    if (!b) {
      return JUNCTION;
    }
  }
  return MANIFOLD;
}

int HalfEdgeMesh::IdentifyBadVerts() {
  int count = 0;
  //list of Faces around each vertex
  std::vector<std::vector<unsigned> > vertFace;
  unsigned numFaces = m_.NumFaces();
  unsigned numVerts = m_.NumVerts();
  vertFace.resize(numVerts);
  for (unsigned f = 0; f < numFaces; f++) {
    for (unsigned j = 0; j < m_.F(f).n; j++) {
      vertFace[m_.F(f)[j]].push_back(f);
    }
  }
  badVertIdx.resize(numVerts, HalfEdge::INVALID_IDX);
  for (unsigned v = 0; v < numVerts; v++) {
    int ret = CheckVertexManifold(v, vertFace[v]);
    if (ret == MANIFOLD) {
      continue;
    }
    badVertIdx[v] = badVertEdges.size();
    std::vector<unsigned> edges(vertFace[v].size());
    for (unsigned i = 0; i < edges.size();i++) {
      unsigned f = vertFace[v][i];
      unsigned j = 0;
      for (; j < m_.F(f).n; j++) {
        if (m_.F(f)[j] == v) {
          break;
        }
      }
      edges[i] = faceEdgeIndex[f] + j;
    }
    badVertEdges.push_back(edges);
    count++;
  }
  return count;
}

unsigned HalfEdgeMesh::Face(unsigned halfEdgeIdx) const { return halfEdgeIdx / 3; }

unsigned HalfEdgeMesh::Next(unsigned halfEdgeIdx) const {
  unsigned t = halfEdgeIdx / 3;
  unsigned vi = (halfEdgeIdx + 1) % 3;
  return 3 * t + vi;
}

unsigned HalfEdgeMesh::Vert(unsigned halfEdgeIdx) const {
  unsigned t = Face(halfEdgeIdx);
  unsigned vi = halfEdgeIdx % 3;
  return m_.F()[t][vi];
}

std::vector<unsigned> HalfEdgeMesh::VertEdges(unsigned vi) const {
  if (IsBadVert(vi)) {
    return badVertEdges[badVertIdx[vi]];
  }
  std::vector<unsigned> edges;
  unsigned e0 = v_[vi];
  const HalfEdge& h = he[e0];
  edges.push_back(e0);
  unsigned e= he[Next(Next(e0))].twin;
  const unsigned MAX_COUNT = 1000000;
  while (e!=HalfEdge::INVALID_IDX && e!=e0) {
    edges.push_back(e);
    e = he[Next(Next(e))].twin;
  }
  return edges;
}
