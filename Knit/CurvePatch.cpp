#include "CurvePatch.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <map>

#include "Array2D.h"
#include "Array3D.h"
#include "Array3DUtils.h"
#include "DisjSet.h"
#include "Edge.h"
#include "ImageIO.h"
#include "ImageUtils.h"
#include "PointGrid3D.h"
#include "TrigMesh.h"
#include "Spline.h"

void SaveCurvePatch(const std::string& file, const CurvePatch& patch) {
  std::ofstream out(file);
  out << "num_curves " << patch.size() << "\n";
  for (size_t i = 0; i < patch.size(); i++) {
    out << patch[i].size() << "\n";
    for (size_t j = 0; j < patch[i].size(); j++) {
      Vec3f v = patch[i][j];
      out << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
  }
}

void SaveCurvePatchObj(const std::string& objfile, const CurvePatch& patch) {
  std::ofstream out(objfile);
  for (size_t i = 0; i < patch.size(); i++) {    
    for (size_t j = 0; j < patch[i].size(); j++) {
      Vec3f v = patch[i][j];
      out <<"v "<< v[0] << " " << v[1] << " " << v[2] << "\n";
    }
  }
}

void BuildUnitPatch() {
  std::string meshFile = "F:/meshes/stitch/singleStitch.obj";
  TrigMesh patch;
  patch.LoadObj(meshFile);
  DisjSet curves;
  BBox box;
  ComputeBBox(patch.v, box);

  std::vector<Vec3f> points(patch.GetNumVerts());
  for (size_t i = 0; i < points.size(); i++) {
    points[i] = patch.GetVertex(i);
  }

  PointGrid3D pointGrid(box, 0.1);
  BuildGrid(pointGrid, points);

  Vec3u size = pointGrid.GetSize();
  // Array3D8u grid;
  // float voxRes[3] = {0.1, 0.1, 0.1};
  // grid.Allocate(size[0], size[1],size[2], 0);
  // const auto& ptGrid = pointGrid.GetGrid();
  // for (size_t i = 0; i < grid.GetData().size(); i++) {
  //   grid.GetData()[i] = ptGrid.GetData()[i].size() > 0;
  // }
  curves.Create(points.size());
  float distThresh = 0.2f;
  std::set<Edge> edges;
  for (size_t i = 0; i < points.size(); i++) {
    std::vector<size_t> nbrs =
        GetNeighbors(points, int(i), distThresh, pointGrid);
    if (nbrs.size() == 0) {
      continue;
    }
    float minDist = (points[i] - points[nbrs[0]]).norm2();
    size_t closestIndex = nbrs[0];
    for (size_t j = 1; j < nbrs.size(); j++) {
      float dist = (points[i] - points[nbrs[j]]).norm2();
      if (dist < minDist) {
        minDist = dist;
        closestIndex = nbrs[j];
      }
    }
    curves.Union(closestIndex, i);
    edges.insert(Edge(i, closestIndex));
    // find second closest if there is one
    if (nbrs.size() == 1) {
      continue;
    }
    bool found = false;
    size_t secondClose = 0;
    float secondMin = 1.5 * minDist;
    for (size_t j = 0; j < nbrs.size(); j++) {
      if (nbrs[j] == closestIndex) {
        continue;
      }
      float dist = (points[i] - points[nbrs[j]]).norm2();
      if (dist < secondMin) {
        found = true;
        secondClose = nbrs[j];
        secondMin = dist;
      }
    }
    if (found) {
      curves.Union(i, secondClose);
    }
    edges.insert(Edge(i, secondClose));
  }

  std::vector<unsigned> roots = curves.GetRoots();
  std::vector < std::vector<unsigned> > adjList(points.size());
  std::vector<std::vector<unsigned> > paths(roots.size());

  for (const auto& edge : edges) {
    adjList[edge.v[0]].push_back(edge.v[1]);
    adjList[edge.v[1]].push_back(edge.v[0]);
  }

  for (size_t c = 0; c < roots.size(); c++) {
    std::deque<unsigned> l;
    // move along root and adjList[c][0] till the end
    unsigned v = roots[c];
    unsigned prev = v;
    if (adjList[v].size() > 1) {
      prev = adjList[v][1];
    }
    while (1) {
      l.push_back(v);
      unsigned next = adjList[v][0];
      if (next == prev) {
        if (adjList[v].size() < 2) {
          break;
        }
        next = adjList[v][1];
      }
      prev = v;
      v = next;
    }
    v = roots[c];
    if (adjList[v].size() == 1) {
      continue;
    }
    // walk the other way
    prev = v;
    v = adjList[v][1];
    while (1) {
      l.push_front(v);
      unsigned next = adjList[v][0];
      if (next == prev) {
        if (adjList[v].size() < 2) {
          break;
        }
        next = adjList[v][1];
      }
      prev = v;
      v = next;
    }
    paths[c] = std::vector<unsigned>(l.begin(), l.end());
  }
  //std::ofstream out("F:/meshes/stitch/lines.obj");
  //for (size_t i = 0; i < points.size(); i++) {
  //  out << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2]
  //      << "\n";
  //}

  //for (size_t c = 0; c < roots.size(); c++) {
  //  for (size_t i = 0; i < paths[c].size() - 1; i++) {
  //    out << "l " << int(1 + paths[c][i]) << " " << int(1 + paths[c][i + 1])
  //        << "\n";
  //  }
  //}
  //out.close();

  std::vector<std::vector<Vec3f> > cpatch(paths.size());
  for (size_t i = 0; i < paths.size(); i++) {
    cpatch[i].resize(paths[i].size());
    for (size_t j = 0; j < paths[i].size(); j++) {
      cpatch[i][j] = points[paths[i][j]];
    }
  }
  SaveCurvePatch("F:/meshes/stitch/curve_patch.txt", cpatch);
}

CurvePatch LoadCurvePatch(const std::string& filename) {
  std::ifstream in(filename);
  CurvePatch curves;
  size_t numCurves;
  std::string token;
  in >> token >> numCurves;
  curves.resize(numCurves);
  for (size_t i = 0; i < numCurves; i++) {
    size_t numPoints;
    in >> numPoints;
    curves[i].resize(numPoints);
    for (size_t j = 0; j < numPoints; j++) {
      Vec3f p;
      in >> p[0] >> p[1] >> p[2];
      curves[i][j] = p;
    }
  }
  return curves;
}

void DrawCurvePatch(const CurvePatch& curvePatch) {
  Array2D8u image(800,800);
  float scale = 100;
  for (size_t c = 0; c < curvePatch.size();c++) {
    for (size_t i = 0; i < curvePatch[c].size() - 1; i++) {
      unsigned x0 = (4 + curvePatch[c][i][0]) * scale;
      unsigned y0 = (4 + curvePatch[c][i][1]) * scale;
      uint8_t color = 20 * (curvePatch[c][i][2] + 4);
      unsigned x1 = (4 + curvePatch[c][i + 1][0]) * scale;
      unsigned y1 = (4 + curvePatch[c][i + 1][1]) * scale;
      DrawLine(image, Vec2f(x0, y0), Vec2f(x1, y1), color);
    }
  }
  SavePngGrey("F:/meshes/stitch/curve_patch.png",image);
}

struct DoubleLinkList {
  struct Node {
    Node() {} 
    Node(int k) : key(k) {}
    int key = -1;
    int prev = -1;
    int next = -1;
  };
  
  int head = -1, tail = -1;
  std::map<int, std::shared_ptr<Node> > nodes;
  
  void Connect(std::shared_ptr<Node> a, std::shared_ptr<Node> b) { 
    a->next = b->key;
    b->prev = a->key;
  }

  void InsertAfter(std::shared_ptr<Node> pos, std::shared_ptr<Node> newNode) {
    int c = pos->next;
    Connect(pos, newNode);
    if (c >= 0) {
      Connect(newNode, nodes[c]);
    } else {
      tail = newNode->key;
    }
    nodes[newNode->key] = newNode;
  }
  void InsertAfter(int key, std::shared_ptr<Node> n) {
    auto prev = LookupByKey(key);
    if (prev) {
      InsertAfter(prev, n);
    }
  }
  void InsertStart(std::shared_ptr<Node> newNode) {
    int a = head;
    if (a<0) {
      head = newNode->key;
      tail = newNode->key;
      nodes[newNode->key] = newNode;
      return;
    }    
    
    Connect(newNode, nodes[head]);
    head = newNode->key;
    nodes[newNode->key] = newNode;
  }

  std::shared_ptr<Node> LookupByKey(int key) {
    if (!nodes.contains(key)) {
      return nullptr;
    }
    return nodes[key];
  }
  int Next(int key) { 
    if (!nodes.contains(key)) {
      return -1;
    }
    return nodes[key]->next; 
  }
};

void TestLinkList() {
  using Node = DoubleLinkList::Node;
  DoubleLinkList ll;
  std::shared_ptr<Node> n = std::make_shared<Node>(123);
  ll.InsertStart(n);
  n = std::make_shared<Node>(456);
  ll.InsertStart(n);
  n = std::make_shared<Node>(789);
  ll.InsertAfter(456, n);
  n = std::make_shared<Node>(101112);
  ll.InsertAfter(123, n);
  int key = ll.head;
  while (key >= 0) {
    std::cout << key << " ";
    key = ll.Next(key);
  }
  std::cout << key << "\n";
}


enum class BoundarySide { NEG_X = 0, POS_X, NEG_Y, POS_Y };

int FindEndpointBoundary(const Vec3f& pt, const BBox& tileBox) {
  std::vector<float> boundaryDists(4, 0);
  boundaryDists[0] = std::abs(pt[0] - tileBox.vmin[0]);
  boundaryDists[1] = std::abs(pt[0] - tileBox.vmax[0]);
  boundaryDists[2] = std::abs(pt[1] - tileBox.vmin[1]);
  boundaryDists[3] = std::abs(pt[1] - tileBox.vmax[1]);
  float minDist = boundaryDists[0];
  int minIndex = 0;
  for (size_t i = 1; i < boundaryDists.size(); i++) {
    if (boundaryDists[i] < minDist) {
      minDist = boundaryDists[i];
      minIndex = int(i);
    }
  }
  return minIndex;
}

bool PointXYInBox(const Vec3f& pt, const BBox& box) {
  return pt[0] >= box.vmin[0] && pt[0] <= box.vmax[0] && pt[1] >= box.vmin[1] &&
         pt[1] <= box.vmax[1];
}

Vec3f SnapPointToBox(const Vec3f& pt, const BBox& box) {
  int boundary = FindEndpointBoundary(pt, box);
  Vec3f out = pt;
  switch (boundary) {
    case int(BoundarySide::NEG_X):
      out = Vec3f(box.vmin[0], pt[1], pt[2]);
      break;
    case int(BoundarySide::POS_X):
      out = Vec3f(box.vmax[0], pt[1], pt[2]);
      break;
    case int(BoundarySide::NEG_Y):
      out = Vec3f(pt[0], box.vmin[1], pt[2]);
      break;
    case int(BoundarySide::POS_Y):
      out = Vec3f(pt[0], box.vmax[1], pt[2]);
      break;
  }
  return out;
}

std::vector<Vec3f> SnapCurveInBox(const std::vector<Vec3f>& curve,
                                  const BBox& box) {
  std::vector<Vec3f> out;

  // i j first and last point on a curve inside box.
  int i = 0;
  for (; i < int(curve.size()); i++) {
    if (PointXYInBox(curve[i], box)) {
      break;
    }
  }
  int j = curve.size() - 1;
  for (; j >= 0; j--) {
    if (PointXYInBox(curve[j], box)) {
      break;
    }
  }
  // add a new start and a new end on the boundary
  out.resize(j - i + 3);
  out[0] = SnapPointToBox(curve[i], box);
  for (size_t k = i; k <= j; k++) {
    out[k - i + 1] = curve[k];
  }
  out.back() = SnapPointToBox(curve[j], box);
  return out;
}

void MoveCurveEndPoint(std::vector<Vec3f>& curve, int endPointIndex,
                       const Vec3f& targetPoint, const BBox& box) {
  Vec3f p0 = curve[endPointIndex];
  Vec3f disp = targetPoint - p0;
  Vec3f boxSize = box.vmax - box.vmin;
  if (std::abs(disp[0]) > 0.2f * boxSize[0]) {
    disp[0] = 0;
  } else if (std::abs(disp[1]) > 0.2f * boxSize[1]) {
    disp[1] = 0;
  }

  const int PointsToMove = 10;
  if (endPointIndex == 0) {
    for (int i = 0; i < PointsToMove; i++) {
      if (i >= curve.size()) {
        break;
      }
      float weight = 1.0f - 0.1f * i;
      curve[i] += weight * disp;
    }
  } else {
    for (int i = endPointIndex; i > endPointIndex - 10; i--) {
      if (i < 0) {
        break;
      }
      float weight = 1.0f - 0.1f * (endPointIndex - i);
      curve[i] += weight * disp;
    }
  }
}

CurvePatch MakeTilable(const CurvePatch& patch, const BBox& tileBox) {
  CurvePatch out(patch.size());
  // find matching pairs of curve beginnings and ends
  struct EndPoint {
    int curveId = 0;
    int pointId = 0;
    int boundary = 0;
    Vec3f pos;
  };

  std::vector<EndPoint> ends(2 * patch.size());
  for (size_t i = 0; i < patch.size(); i++) {
    std::vector<Vec3f> newCurve;
    Vec3f start = patch[i][0];
    Vec3f end = patch[i].back();
    ends[2 * i].boundary = FindEndpointBoundary(start, tileBox);
    ends[2 * i + 1].boundary = FindEndpointBoundary(end, tileBox);
    out[i] = SnapCurveInBox(patch[i], tileBox);

    ends[2 * i].curveId = int(i);
    ends[2 * i].pointId = 0;
    ends[2 * i].pos = out[i][0];
    ends[2 * i + 1].curveId = int(i);
    ends[2 * i + 1].pointId = int(out[i].size() - 1);
    ends[2 * i + 1].pos = out[i].back();
  }

  // match pairs in x edges and y edges
  std::vector<std::vector<unsigned> > boundaryPoints(4);
  for (size_t i = 0; i < ends.size(); i++) {
    boundaryPoints[ends[i].boundary].push_back(i);
  }

  std::vector<Edge> pairs;
  // pairs in x
  for (size_t i = 0; i < boundaryPoints[0].size(); i++) {
    int closest = 0;
    EndPoint p0 = ends[boundaryPoints[0][i]];
    float minDist = -1;
    for (size_t j = 0; j < boundaryPoints[1].size(); j++) {
      EndPoint p1 = ends[boundaryPoints[1][j]];
      Vec3f pos1 = p1.pos;
      pos1[0] = p0.pos[0];
      float dist = (pos1 - p0.pos).norm();
      if (minDist < 0 || dist < minDist) {
        minDist = dist;
        closest = int(boundaryPoints[1][j]);
      }
    }
    pairs.push_back(Edge(boundaryPoints[0][i], closest));
  }
  // pairs in y
  for (size_t i = 0; i < boundaryPoints[2].size(); i++) {
    int closest = 0;
    EndPoint p0 = ends[boundaryPoints[2][i]];
    float minDist = -1;
    for (size_t j = 0; j < boundaryPoints[3].size(); j++) {
      EndPoint p1 = ends[boundaryPoints[3][j]];
      Vec3f pos1 = p1.pos;
      pos1[1] = p0.pos[1];
      float dist = (pos1 - p0.pos).norm();
      if (minDist < 0 || dist < minDist) {
        minDist = dist;
        closest = int(boundaryPoints[3][j]);
      }
    }
    pairs.push_back(Edge(boundaryPoints[2][i], closest));
  }

  // move end points to middle points of each edge
  for (size_t ei = 0; ei < pairs.size(); ei++) {
    EndPoint e0 = ends[pairs[ei].v[0]];
    EndPoint e1 = ends[pairs[ei].v[1]];
    Vec3f mid = 0.5f * (e0.pos + e1.pos);
    MoveCurveEndPoint(out[e0.curveId], e0.pointId, mid, tileBox);
    MoveCurveEndPoint(out[e1.curveId], e1.pointId, mid, tileBox);
  }

    
  //std::ofstream outf("F:/meshes/stitch/tile.obj");
  //for (size_t i = 0; i < out.size(); i++) {
  //  for (size_t j = 0; j < out[i].size(); j++) {
  //    outf << "v " << out[i][j][0] << " " << out[i][j][1] << " " << out[i][j][2]
  //        << "\n";
  //  }
  //}
  //outf.close();

  return out;
}

void MakeTilablePatch() {
  std::vector<std::vector<Vec3f> > patch =
      LoadCurvePatch("F:/meshes/stitch/curve_patch.txt");
  BBox tileBox;
  tileBox.vmin = Vec3f(-5.9761, -5.1287, 0);
  tileBox.vmax = Vec3f(5.9688, 5.0846, 0);
  tileBox.vmin = 0.5f * tileBox.vmin;
  tileBox.vmax = 0.5f * tileBox.vmax;

  std::vector<std::vector<Vec3f> > tile = MakeTilable(patch, tileBox);
  SaveCurvePatch("F:/meshes/stitch/stitchTile.txt", tile);
  DrawCurvePatch(tile);
}

std::vector<Vec3f> ComputeFiller(const std::vector<Vec3f> & curve, unsigned i0, float delta) {
  std::vector<Vec3f> filler;
  SplineCR1 spline;
  spline.p[1] = curve[i0];
  spline.p[2] = curve[i0 + 1];
  if (i0 > 0) {
    spline.p[0] = curve[i0 - 1];
  } else {
    spline.p[0] = spline.p[1] + delta * (spline.p[1] - spline.p[2]);
  }
  if (i0 < curve.size() - 2) {
    spline.p[3] = curve[i0 + 2];
  } else {
    spline.p[3] = spline.p[2] + delta * (spline.p[2] - spline.p[1]);
  }
  float dist12 = (spline.p[2] - spline.p[1]).norm();
  unsigned numSamples = unsigned(dist12 / delta) + 1;
  float dt = 1.0f / numSamples;
  for (unsigned s = 0; s < numSamples; s++) {
    Vec3f sample = spline.Eval(s * dt);
    filler.push_back(sample);
  }
  return filler;
}

void FillInGaps(CurvePatch& patch, float maxDist) {
  float delta = 0.9f * maxDist;
  for (size_t c = 0; c < patch.size(); c++) {
    std::vector<Vec3f> newCurve;
    for (size_t i = 0; i < patch[c].size() - 1; i++) {
      Vec3f edge = patch[c][i + 1] - patch[c][i];
      if (edge.norm2() > maxDist * maxDist) {
        std::vector<Vec3f> filler = ComputeFiller(patch[c], i, delta);
        newCurve.insert(newCurve.end(), filler.begin(), filler.end());
      } else {
        newCurve.push_back(patch[c][i]);
      }      
    }
    newCurve.push_back(patch[c].back());
    patch[c] = newCurve;
  }
}