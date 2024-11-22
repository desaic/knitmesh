#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "Array2D.h"
#include "Array3D.h"
#include "Array3DUtils.h"
#include "BBox.h"
#include "CurvePatch.h"
#include "DisjSet.h"
#include "Edge.h"
#include "HalfEdgeMesh.h"
#include "ImageIO.h"
#include "ImageUtils.h"
#include "Matrix4f.h"
#include "MeshUtil.h"
#include "PointGrid3D.h"
#include "PolyMesh.h"
#include "TrigMesh.h"

Vec3f PerspectiveTrans(const Matrix4f& m, const Vec3f& v) {
  Vec4f h(v[0], v[1], v[2], 1);
  h = m * h;
  if (h[3] != 0) {
    h[0] /= h[3];
    h[1] /= h[3];
    h[2] /= h[3];
  }
  return Vec3f(h[0], h[1], h[2]);
}

/// <summary>
/// Valid return coordinates are between (0,0) to (width-1, height-1).
/// </summary>
/// <param name="projCoord">only valid within [-1,-1] to [1,1]</param>
/// <param name="width"></param>
/// <param name="height"></param>
/// <returns></returns>
Vec2f ImageCoord(const Vec3f& projCoord, unsigned width, unsigned height) {
  Vec2f coord;
  coord[0] = width * 0.5f * (projCoord[0] + 1);
  coord[1] = height * 0.5f * (projCoord[1] + 1);
  return coord;
}

bool InBound(const Vec2i& index, const Vec2u& imageSize) {
  return index[0] >= 0 && index[0] < int(imageSize[0]) && index[1] >= 0 &&
         index[1] < int(imageSize[1]);
}

void DrawHEMesh(const HalfEdgeMesh& hem) {
  Array2D8u image(800, 800);
  Vec2u imageSize = image.GetSize();
  BBox box;
  ComputeBBox(hem.V(), box);
  Vec3f center = 0.5f * (box.vmax + box.vmin);
  Vec3f boxSize = box.vmax - box.vmin;
  float fovRad = 1.57f;
  float aspect = 1.0f;
  float near = 0.1f;
  float far = near + 10 * boxSize.norm();
  Vec3f eye = center + 2.0f * boxSize * Vec3f(0.1, 0.1, 0.5);
  Matrix4f proj =
      Matrix4f::perspectiveProjection(fovRad, aspect, near, far, false);
  Matrix4f view = Matrix4f::lookAt(eye, center, Vec3f(0, 1, 0));
  Matrix4f prod = proj * view;
  for (size_t i = 0; i < hem.he.size(); i++) {
    const HalfEdge& he = hem.he[i];
    unsigned v0 = he.vertex;
    unsigned v1 = hem.he[he.twin].vertex;
    Vec3f x0 = hem.V(v0);
    Vec3f x1 = hem.V(v1);
    x0 = PerspectiveTrans(prod, x0);
    x1 = PerspectiveTrans(prod, x1);
    Vec2f imageIdx0 = ImageCoord(x0, imageSize[0], imageSize[1]);
    Vec2f imageIdx1 = ImageCoord(x1, imageSize[0], imageSize[1]);
    DrawLine(image, imageIdx0, imageIdx1, 255);
  }

  SavePngGrey("hem.png", image);
}

void FloodLabels(HalfEdgeMesh& hem, size_t fi) {
  std::deque<unsigned> q;
  q.push_back(fi);
  while (!q.empty()) {
    size_t faceIndex = q.front();
    q.pop_front();
    const PolyFace& f = hem.m_.F(faceIndex);
    unsigned numEdges = f.n;
    HalfEdge he = hem.he[hem.faceEdgeIndex[faceIndex]];
    // which orientation is most common among neighbors.
    std::vector<unsigned> orientVote(4, 0);
    for (unsigned i = 0; i < numEdges; i++) {
      unsigned twinIndex = he.twin;
      unsigned neighborFace = hem.edgeFaceIndex[twinIndex];
      int neighborLabel = hem.faceLabels[neighborFace];
      if (neighborLabel < 0) {
        // neighbor unlabeled. add to queue. skip voting.
        if (neighborLabel == -1) {
          hem.faceLabels[neighborFace] = -2;
          q.push_back(neighborFace);
        }
        he = hem.he[he.next];
        continue;
      }
      // local index of the shared edge in neighbor face
      unsigned nei = hem.EdgeIndexInFace(twinIndex, neighborFace);

      // label relative to neighbor's shared edge
      neighborLabel = (neighborLabel + (HalfEdgeMesh::MAX_LABEL - nei)) %
                      HalfEdgeMesh::MAX_LABEL;
      // label relative to this face's shared edge.
      // flips because direction of half edge flips.
      int myLabel = (HalfEdgeMesh::MAX_LABEL / 2 + neighborLabel) %
                    HalfEdgeMesh::MAX_LABEL;
      // label relative to edge 0 of this face.
      myLabel = (myLabel + i) % HalfEdgeMesh::MAX_LABEL;
      orientVote[myLabel]++;
      he = hem.he[he.next];
    }

    unsigned maxCount = 0;
    unsigned maxLabel = 0;
    for (size_t i = 0; i < orientVote.size(); i++) {
      if (orientVote[i] > maxCount) {
        maxCount = orientVote[i];
        maxLabel = i;
      }
    }
    hem.faceLabels[faceIndex] = int(maxLabel);
  }
}

void LabelOrientation(HalfEdgeMesh& hem) {
  size_t nF = hem.NumFaces();
  hem.faceLabels.resize(hem.NumFaces(), -1);
  for (size_t i = 0; i < nF; i++) {
    if (hem.faceLabels[i] >= 0) {
      // already labeled.
      continue;
    }
    FloodLabels(hem, i);
  }
}

// save half edge mesh as obj.
// save the labels as arrows made of line segments
void SaveLabeledMeshObj(const std::string& outFile, const HalfEdgeMesh& hem) {
  std::ofstream out(outFile);
  for (size_t i = 0; i < hem.NumVerts(); i++) {
    out << "v " << hem.V(i)[0] << " " << hem.V(i)[1] << " " << hem.V(i)[2]
        << "\n";
  }
  for (size_t i = 0; i < hem.NumFaces(); i++) {
    out << "f " << (hem.F(i)[0] + 1) << " " << (hem.F(i)[1] + 1) << " "
        << (hem.F(i)[2] + 1) << " " << (hem.F(i)[3] + 1) << "\n";
  }
  if (hem.faceLabels.size() != hem.NumFaces()) {
    std::cout << "no labels\n";
    return;
  }
  unsigned vertCount = hem.NumVerts();
  for (size_t i = 0; i < hem.NumFaces(); i++) {
    int label = hem.faceLabels[i];
    if (label < 0) {
      continue;
    }
    if (hem.F(i).n != 4) {
      // only knows quads for now.
      continue;
    }
    Vec3f quadCenter(0);
    const PolyFace& face = hem.F(i);
    for (size_t j = 0; j < face.n; j++) {
      quadCenter += hem.V(face[j]);
    }
    quadCenter = (1.0f / face.n) * quadCenter;
    Vec3f edge = hem.V(face[(label + 1) % face.n]) - hem.V(face[label]);
    edge = 0.4f * edge;
    Vec3f p = quadCenter + edge;
    out << "v " << quadCenter[0] << " " << quadCenter[1] << " " << quadCenter[2]
        << "\n";
    out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    out << "l " << vertCount + 1 << " " << vertCount + 2 << "\n";
    vertCount += 2;
  }

  //draw normal lines
  if (hem.vn.size() == hem.NumVerts()) {
    for (size_t i = 0; i < hem.NumVerts(); i++) {
      Vec3f v = hem.V(i);
      v += hem.vn[i];
      out << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
      out << "l " << (i + 1) << " " << (vertCount + 1) <<"\n";
      vertCount++;
    }
  }
}

CurvePatch ConvertCurvesToNaturalCoordinates(const CurvePatch& curves) {
  CurvePatch out;
  if (curves.size() == 0 || curves[0].size() == 0) {
    return out;
  }
  BBox box;
  box.vmin = curves[0][0];
  box.vmax = curves[0][0];
  for (size_t c = 0; c < curves.size(); c++) {
    for (size_t i = 0; i < curves[c].size(); i++) {
      const auto& v = curves[c][i];
      for (unsigned d = 0; d < 3; d++) {
        box.vmin[d] = std::min(box.vmin[d], v[d]);
        box.vmax[d] = std::max(box.vmax[d], v[d]);
      }
    }
  }

  out.resize(curves.size());
  Vec3f center = 0.5f * (box.vmax + box.vmin);
  Vec3f scale;
  scale[0] = 2.0f / (box.vmax[0] - box.vmin[0]);
  scale[1] = 2.0f / (box.vmax[1] - box.vmin[1]);
  scale[2] = 2.0f * scale[0];
  for (size_t c = 0; c < curves.size(); c++) {
    out[c].resize(curves[c].size());
    for (size_t i = 0; i < curves[c].size(); i++) {
      auto& v = out[c][i];
      const auto& vin = curves[c][i];
      for (unsigned d = 0; d < 3; d++) {
        v[d] = scale[d] * (vin[d] - center[d]);
      }
    }
  }

  box.vmin = out[0][0];
  box.vmax = out[0][0];
  for (size_t c = 0; c < curves.size(); c++) {
    for (size_t i = 0; i < curves[c].size(); i++) {
      const auto& v = out[c][i];
      for (unsigned d = 0; d < 3; d++) {
        box.vmin[d] = std::min(box.vmin[d], v[d]);
        box.vmax[d] = std::max(box.vmax[d], v[d]);
      }
    }
  }
  //  std::cout << box.vmin[0] << " " << box.vmax[0] << "\n";
  // std::cout << box.vmin[1] << " " << box.vmax[1] << "\n";
  // std::cout << box.vmin[2] << " " << box.vmax[2] << "\n";

  return out;
}

//quad with vertex normals
class Quad {
 public:
  //4 vertices in xy plane. z is centered.
  std::array<Vec3f, 4> v;
  //positive z direction for each vertex.
  std::array<Vec3f, 4> n;
  unsigned size() const { return v.size(); }
};

class Hex {
 public:
  std::array<Vec3f, 8> v;
  Vec3f& operator[](unsigned i) { return v[i]; }
  const Vec3f& operator[](unsigned i) const { return v[i]; }
  unsigned size() const { return v.size(); }
};

Hex QuadToHex(const Quad& q) {
  Hex h;
  for (unsigned i = 0; i < q.size(); i++) {
    h[i] = q.v[i] - q.n[i];
    h[i + q.size()] = q.v[i] + q.n[i];
  }
  return h;
}

float HexShapeFun(const Vec3f& X, unsigned vi) {
  const float N[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                         {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  return 0.125f * (1 + N[vi][0] * X[0]) * (1 + N[vi][1] * X[1]) *
         (1 + N[vi][2] * X[2]);
}

Vec3f TrilinearInterp(const Vec3f& x, const Hex& hex) {
  Vec3f vout(0);
  for (unsigned i = 0; i < hex.size(); i++) {
    vout += hex.v[i] * HexShapeFun(x, i);
  }
  return vout;
}

CurvePatch TrilinearInterp(const CurvePatch& curves, const Hex& hex) {
  CurvePatch patch(curves.size());
  for (size_t c = 0; c < patch.size(); c++) {
    patch[c].resize(curves[c].size());
    for (size_t i = 0; i < curves[c].size(); i++) {
      patch[c][i] = TrilinearInterp(curves[c][i], hex);
    }
  }
  return patch;
}

std::vector<CurvePatch> PutCurvesOnLabeledMesh(const CurvePatch& patch,
  const HalfEdgeMesh& hem, std::vector<PatchModifier>& mods) {
  std::vector<CurvePatch> faceCurves(hem.NumFaces());

  std::ofstream dout("F:/meshes/stitch/debugPoints.obj");

  for (size_t i = 0; i < hem.NumFaces(); i++) {
    Quad quad;
    CurvePatch curves = patch;
    const PolyFace& face = hem.F(i);
    if (face.n != quad.size()) {
      //only knows how to deal with quads
      continue;
    }
    int lab = hem.faceLabels[i];
    if (lab < 0) {
      //should never be
      lab = 0;
    }
    //go though 4 neighbors and modify points on each edge
    auto he = hem.he[hem.faceEdgeIndex[i]];
    for (unsigned n = 0; n < face.n; n++) {
      if (!he.hasTwin()) {
        //@TODO handle boundary edges.
        he = hem.he[he.next];
        continue;
      }
      int myEdge = (n + 4 - lab) % 4;
      unsigned neighbor = hem.edgeFaceIndex[he.twin];
      int neighborLab = hem.faceLabels[neighbor];
      if (neighborLab < 0) {
        //impossible
        he = hem.he[he.next];
        continue;
      }
      int neighborEdgeIndex = hem.EdgeIndexInFace(he.twin, neighbor);
      int neighborEdge = (neighborEdgeIndex + 4 - neighborLab) % 4;
      for (size_t mi = 0; mi < mods.size(); mi++) {
        if (mods[mi].myEdge == myEdge && mods[mi].neighborEdge == neighborEdge) {
          mods[mi].Apply(curves);
          //only one modifier can be applied to each edge.
          break;
        }
      }
      he = hem.he[he.next];
    }

    for (unsigned j = 0; j < quad.size(); j++) {
      unsigned vIdx = face[(j + lab) % quad.size()];
      quad.v[j] = hem.V(vIdx);
      quad.n[j] = hem.vn[vIdx];
    }
    Hex hex = QuadToHex(quad);
    faceCurves[i] = TrilinearInterp(curves, hex);
    for (unsigned c = 0; c < faceCurves[i].size(); c++) {
      for (unsigned j = 0; j < faceCurves[i][c].size(); j++) {
        Vec3f p = faceCurves[i][c][j];
        dout << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
      }
    }
  }
  return faceCurves;
}

void ComputeVertNormals(HalfEdgeMesh& hem) {
  hem.vn.resize(hem.NumVerts(), Vec3f(0));
  for (size_t fi = 0; fi < hem.NumFaces(); fi++) {
    const auto& face = hem.F(fi);
    Vec3f n(0);
    for (unsigned vi = 0; vi < face.n - 2; vi++) {
      // not numerically stable. whatever.
      Vec3f e1 = hem.V(face[vi + 1]) - hem.V(face[vi]);
      Vec3f e2 = hem.V(face[vi + 2]) - hem.V(face[vi]);
      Vec3f cross = e1.cross(e2);      
      n += cross;
    }
    n.normalize();
    for (unsigned vi = 0; vi < face.n; vi++) {
      hem.vn[face[vi]] += n;
    }
  }
  for (unsigned v = 0; v < hem.NumVerts(); v++) {
    hem.vn[v].normalize();
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: stitcher meshFile.obj\n";
    return 0;
  }

  std::string meshFile = argv[1];
  PolyMesh pm;
  pm.LoadObj(argv[1]);
  std::cout << pm.NumVerts() << " " << pm.NumFaces() << "\n";

  HalfEdgeMesh hem(pm);
  // DrawHEMesh(hem);
  LabelOrientation(hem);
  ComputeVertNormals(hem);
  //SaveLabeledMeshObj("F:/meshes/stitch/labeled.obj", hem);
  auto curves =
      LoadCurvePatch("F:/github/knitmesh/Data/scaled_patch.txt");
  std::vector<PatchModifier> mods = GeneratePatchModifiers(curves);

  std::vector<CurvePatch> meshCurves = PutCurvesOnLabeledMesh(curves, hem, mods);

  return 0;
}
