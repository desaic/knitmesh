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

  CurvePatch cpatch(paths.size());
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
  if (!in.good()) {
    return CurvePatch();
  }
  size_t numCurves;
  std::string token;
  in >> token >> numCurves;
  CurvePatch curves(numCurves);
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

//find which side each endpoint is on.
std::vector<EndPoint> FindEndPointSides(const CurvePatch&patch, const BBox & box) {
  std::vector<EndPoint> ends(2 * patch.size());
  for (size_t i = 0; i < patch.size(); i++) {
    std::vector<Vec3f> newCurve;
    Vec3f start = patch[i][0];
    Vec3f end = patch[i].back();
    ends[2 * i].boundary = FindEndpointBoundary(start, box);
    ends[2 * i + 1].boundary = FindEndpointBoundary(end, box);
    ends[2 * i].curveId = int(i);
    ends[2 * i].pointId = 0;
    ends[2 * i].pos = start;
    ends[2 * i + 1].curveId = int(i);
    ends[2 * i + 1].pointId = int(patch[i].size() - 1);
    ends[2 * i + 1].pos = end;
  }
  return ends;
}

CurvePatch MakeTilable(const CurvePatch& patch, const BBox& tileBox) {
  CurvePatch out(patch.size());
  // find matching pairs of curve beginnings and ends
  for (size_t i = 0; i < patch.size(); i++){
    out[i] = SnapCurveInBox(patch[i], tileBox);
  }
  std::vector<EndPoint> ends = FindEndPointSides(out, tileBox);

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
  CurvePatch patch =
      LoadCurvePatch("F:/meshes/stitch/curve_patch.txt");
  BBox tileBox;
  tileBox.vmin = Vec3f(-5.9761, -5.1287, 0);
  tileBox.vmax = Vec3f(5.9688, 5.0846, 0);
  tileBox.vmin = 0.5f * tileBox.vmin;
  tileBox.vmax = 0.5f * tileBox.vmax;

  CurvePatch tile = MakeTilable(patch, tileBox);
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

CurvePatch RotateCW(const CurvePatch& patch, int quarterTurns) {
  CurvePatch out(patch.size());
  out.c = patch.c;
  for (int turn = 0; turn < quarterTurns; turn++) {
    for (size_t c = 0; c < out.size(); c++) {
      auto & curve = out[c];
      for (size_t i = 0; i < curve.size(); i++) {
        float x = curve[i][0];
        float y = curve[i][1];
        curve[i][0] = y;
        curve[i][1] = -x;
      }
    }
  }
  out.ends = patch.ends;
  for (size_t i = 0; i < out.ends.size(); i++) {
    int c = out.ends[i].curveId;
    int p = out.ends[i].pointId;
    out.ends[i].pos = out[c][p];
    if (out.ends[i].boundary >= 0) {
      out.ends[i].boundary = (out.ends[i].boundary - quarterTurns + 4) % 4;
    }
  }

  return out;
}

void PatchModifier::Apply(CurvePatch& curves) {
  for (const auto & mod:curveMods) {
    auto& curve = curves.c[mod.curveId];
    if (mod.op == CurveMod::OP::REPLACE) {
      if (mod.endPoint % 2 == 0) {
        //starting point
        for (size_t i = 0; i < mod.p.size(); i++) {
          curve[i] = mod.p[i];
        }
      }
      else {
        for (size_t i = 0; i < mod.p.size(); i++) {
          curve[curve.size() - i - 1] = mod.p[i];
        }
      }
    }
    else if(mod.op == CurveMod::OP::CONNECT){
      std::vector<Vec3f > newCurve;
      if (mod.endPoint % 2 == 0) {
        for (size_t i = 0; i < mod.p.size(); i++) {
          newCurve.push_back(mod.p[mod.p.size() - i - 1]);
        }
        newCurve.insert(newCurve.end(), curve.begin(), curve.end());
      }
      else {
        newCurve.insert(newCurve.end(), curve.begin(), curve.end());
        newCurve.insert(newCurve.end(), mod.p.begin(), mod.p.end());
      }
      curve = newCurve;
    }
  }
}

struct CurveMatching {
  //pairs between me and neighbor.
  std::vector<std::pair<int, int>>cross;
  //pair up unmatched end points on my side
  std::vector<std::pair<int, int>>myPair;
  //pair up unpatched neighbor points.
  std::vector<std::pair<int, int>>nbrPair;
};

//match end points between edge 0 of me and edge 2 of nbr in xz plane.
CurveMatching MatchEndPoints(const CurvePatch&me, const CurvePatch & nbr) {
  std::vector<int> myEnds, nbrEnds;
  for (size_t i = 0; i < me.ends.size(); i++) {
    if (me.ends[i].boundary == 0) {
      myEnds.push_back(i);
    }
  }
  for (size_t j = 0; j < nbr.ends.size(); j++) {
    if (nbr.ends[j].boundary == 2) {
      nbrEnds.push_back(j);
    }
  }
  
  CurveMatching matches;
  // is this neighbor endpoint already matched.
  std::vector<uint8_t> nbrMatched(nbrEnds.size(), 0);
  std::vector<uint8_t> myMatched(myEnds.size(), 0);
  
  if (myEnds.size() >= nbrEnds.size()) {
    for (size_t j = 0; j < nbrEnds.size(); j++) {
      float minDist = -1;
      int match = 0;
      for (size_t i = 0; i < myEnds.size(); i++) {
        Vec3f v = me.ends[myEnds[i]].pos - nbr.ends[nbrEnds[j]].pos;
        v[1] = 0;
        float dist = v.norm();
        if (i == 0 || dist < minDist) {
          minDist = dist;
          match = int(i);
        }
      }
      matches.cross.push_back(std::make_pair(myEnds[match], nbrEnds[j]));
      myMatched[match] = 1;
      nbrMatched[j] = 1;
    }
    for (size_t i = 0; i < myEnds.size(); i++) {
      float minDist = -1;
      int match = 0;
      if (myMatched[i]) {
        continue;
      }
      for (size_t j = i+1; j < myEnds.size(); j++) {
        if (myMatched[j]) {
          continue;
        }
        Vec3f v = me.ends[myEnds[i]].pos - me.ends[myEnds[j]].pos;
        v[1] = 0;
        float dist = v.norm();
        if (j == 0 || dist < minDist) {
          minDist = dist;
          match = int(j);
        }
      }
      if (match > 0) {
        matches.myPair.push_back(std::make_pair(myEnds[i], myEnds[match]));
      }
    }
  }
  else {
    for (size_t i = 0; i < myEnds.size(); i++) {
      float minDist = -1;
      int match = 0;
      for (size_t j = 0; j < nbrEnds.size(); j++) {
        Vec3f v = me.ends[myEnds[i]].pos - nbr.ends[nbrEnds[j]].pos;
        v[1] = 0;
        float dist = v.norm();
        if (match == 0 || dist < minDist) {
          minDist = dist;
          match = int(j);
        }
      }
      matches.cross.push_back(std::make_pair(myEnds[i], nbrEnds[match]));
      nbrMatched[match] = 1;
      myMatched[i] = 1;
    }

    for (size_t i = 0; i < nbrEnds.size(); i++) {
      float minDist = -1;
      int match = 0;
      if (nbrMatched[i]) {
        continue;
      }
      for (size_t j = i + 1; j < nbrEnds.size(); j++) {
        if (nbrMatched[j]) {
          continue;
        }
        Vec3f v = nbr.ends[nbrEnds[i]].pos - nbr.ends[nbrEnds[j]].pos;
        v[1] = 0;
        float dist = v.norm();
        if (match == 0 || dist < minDist) {
          minDist = dist;
          match = int(j);
        }
      }
      if (match > 0) {
        matches.nbrPair.push_back(std::make_pair(nbrEnds[i], nbrEnds[match]));
      }
    }
  }
  return matches;
}

std::vector<CurveMod> BlendEndPoints(const CurvePatch & my, int myEndIndex, const CurvePatch&nbr,
  int nbrEndIndex) {
  const auto& myEnd = my.ends[myEndIndex];
  const auto& nbrEnd = nbr.ends[nbrEndIndex];
  std::vector<CurveMod> mods(2);
  Vec3f middle = 0.5f * (myEnd.pos + nbrEnd.pos);
  mods[0].curveId = myEnd.curveId;
  mods[0].endPoint = myEndIndex;
  mods[0].op = CurveMod::OP::REPLACE;
  
  //number of points being moved
  const size_t NUM_MOD_POINTS = 10;
  size_t numPoints = std::min(NUM_MOD_POINTS, my.c[myEnd.curveId].size());
  mods[0].p.resize(numPoints);
  int pointId = myEnd.pointId;
  int step = 1;
  if (myEnd.pointId > 0) {
    step = -1;
  }

  const auto& myc = my.c[myEnd.curveId];
  for (size_t i = 0; i < numPoints; i++) {
    float weight = (numPoints - i) / float(numPoints);    
    mods[0].p[i] = weight * middle + (1 - weight) * myc[pointId];
    mods[0].p[i][1] = myc[pointId][1];
    pointId += step;
  }

  //nbr side
  mods[1].curveId = nbrEnd.curveId;
  mods[1].endPoint = nbrEndIndex;
  mods[1].op = CurveMod::OP::REPLACE;

  numPoints = std::min(NUM_MOD_POINTS, nbr.c[nbrEnd.curveId].size());
  mods[1].p.resize(numPoints);
  pointId = nbrEnd.pointId;
  step = 1;
  if (nbrEnd.pointId > 0) {
    step = -1;
  }

  const auto& nc = nbr.c[nbrEnd.curveId];
  for (size_t i = 0; i < numPoints; i++) {
    float weight = (numPoints - i) / float(numPoints);
    mods[1].p[i] = weight * middle + (1 - weight) * nc[pointId];
    mods[1].p[i][1] = nc[pointId][1];
    pointId += step;
  }

  return mods;
}

CurveMod ConnectEndPoints(const CurvePatch& curves, int endIndex1, int endIndex2) {
  const auto& end1 = curves.ends[endIndex1];
  const auto& end2 = curves.ends[endIndex2];
  CurveMod mod;
  Vec3f middle = 0.5f * (end1.pos + end2.pos);
  mod.curveId = end1.curveId;
  mod.endPoint = endIndex1;
  mod.endPointConn = endIndex2;
  mod.op = CurveMod::OP::CONNECT;

  //number of points being added
  const size_t NUM_MOD_POINTS = 20;  
  mod.p.resize(NUM_MOD_POINTS);
  int pointId = end1.pointId;

  SplineCR1 spline;
    
  Vec3f move;
  float moveDist= 0.05;
  int side = end1.boundary;
  switch (side) {
  case 0:
    move = Vec3f(0, -1, 0);
    break;
  case 1:
    move = Vec3f(1, 0, 0);
    break;
  case 2:
    move = Vec3f(0, 1, 0);
    break;
  case 3:
    move = Vec3f(-1, 0, 0);
    break;
  default:
    break;
  }

  spline.p.resize(5);
  spline.p[0] = end1.pos;
  spline.p[1] = end1.pos + moveDist *move;
  spline.p[2] = middle + 0.2f * move;
  spline.p[3] = end2.pos + moveDist * move;
  spline.p[4] = end2.pos;

  for (size_t i = 0; i < NUM_MOD_POINTS; i++) {
    Vec3f pt = spline.Eval(2 * i / float(NUM_MOD_POINTS));
    mod.p.push_back(pt);
  }
  return mod;
}

std::vector<PatchModifier> GeneratePatchModifiers(const CurvePatch& curves) {
  std::vector<PatchModifier> mods;
  BBox box;
  box._init = true;
  box.vmin = Vec3f(-1, -1, -1);
  box.vmax = Vec3f(1, 1, 1);
  CurvePatch localCopy = curves;
  localCopy.ends = FindEndPointSides(curves, box);
  //     3            2
  //  0     1  =>  3     1
  //     2            0
  // relabel from axis order to edge order.
  const int RELABEL[4] = { 3,1,0,2 };
  for (size_t i = 0; i < localCopy.ends.size(); i++) {
    if (localCopy.ends[i].boundary >= 0) {
      localCopy.ends[i].boundary = RELABEL[localCopy.ends[i].boundary];
    }
  }

  for (int myEdge = 0; myEdge < 3; myEdge++) {
    //rotate this curve patch so that edge i is facing down.
    int myRot = myEdge;
    CurvePatch myCopy = RotateCW(localCopy, myRot);
    
    for (int neighborEdge = myEdge + 1; neighborEdge < 4; neighborEdge++) {
      //number of 90deg clockwise turns so that neighbor's edge is facing up (2).
      int neighborRot = ((neighborEdge - 2) + 4) % 4;

      PatchModifier mod;
      mod.myEdge = myEdge;
      mod.neighborEdge = neighborEdge;
      CurvePatch neighborCopy = RotateCW(localCopy, neighborRot);
      PatchModifier patchMod, nbrMod;
      CurveMatching matches = MatchEndPoints(myCopy, neighborCopy);
      for (auto match : matches.cross) {
        std::vector<CurveMod> curveMods = BlendEndPoints(myCopy, match.first, neighborCopy, match.second);
        patchMod.curveMods.push_back(curveMods[0]);
        nbrMod.curveMods.push_back(curveMods[1]);
      }
      for (auto match : matches.myPair) {
        CurveMod curveMod = ConnectEndPoints(myCopy, match.first,match.second);
        patchMod.curveMods.push_back(curveMod);        
      }
      for (auto match : matches.nbrPair) {
        CurveMod curveMod = ConnectEndPoints(neighborCopy, match.first, match.second);
        nbrMod.curveMods.push_back(curveMod);
      }
      patchMod.Apply(myCopy);
      nbrMod.Apply(neighborCopy);
      SaveCurvePatchObj("F:/meshes/stitch/myCopy" + std::to_string(myEdge) + "_" + std::to_string(neighborEdge) + ".obj", myCopy);
      SaveCurvePatchObj("F:/meshes/stitch/neighborCopy" + std::to_string(myEdge) + "_" + std::to_string(neighborEdge) + ".obj", neighborCopy);

      mods.push_back(mod);

    }
  }
  return mods;
}
