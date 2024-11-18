#include "Spline.h"
#include "BBox.h"

int SplineCR1::NumSegments() const {
  if (isLoop) {
    return p.size();
  } else {
    return int(p.size()) - 3;
  }
}

Vec3f SplineCR1::Eval(float s) const {
  if (isLoop) {
    return EvalLoop(s);
  }
  int numSegs = NumSegments();
  if (numSegs <= 0) {
    return Vec3f(0);
  }
  s = std::max(0.0f, s);
  s = std::min(s, float(numSegs));
  //Figure 3 of the paper 
  Vec3f C12;
  // curve is parameterized s.t.
  //si - s_i-1 = 0
  unsigned i0 = unsigned(s);
  // s-s1
  float s1 = s - i0;
  // s2-s
  float s2 = 1 - s1;
  float s0 = s1 + 1;
  float s3 = s2 + 1;
  Vec3f L01 = -s1 * p[i0] + s0 * p[i0 + 1];
  Vec3f L12 = s2 * p[i0 + 1] + s1 * p[i0 + 2];
  Vec3f L23 = s3 * p[i0 + 2] - s2 * p[i0 + 3];
  
  Vec3f L012 = 0.5f * (s2 * L01 + s0 * L12);
  Vec3f L123 = 0.5f * (s3 * L12 + s1 * L23);
  C12 = s2 * L012 + s1 * L123;
  return C12;
}

Vec3f SplineCR1::EvalLoop(float s) const {
  int numSegs = NumSegments();
  if (numSegs <= 0) {
    return Vec3f(0);
  }
  s = std::max(0.0f, s);
  s = std::min(s, float(numSegs));
  // Figure 3 of the paper
  Vec3f C12;
  // curve is parameterized s.t.
  // si - s_i-1 = 0
  unsigned i1 = std::min(unsigned(numSegs - 1), unsigned(s));
  unsigned i0 = 0;
  if (i1 > 0) {
    i0 = i1 - 1;
  } else {
    i0 = p.size() - 1;
  }
  unsigned i2 = (i1 + 1) % p.size();
  unsigned i3 = (i1 + 2) % p.size();
  // s-s1
  float s1 = s - i1;
  // s2-s
  float s2 = 1 - s1;
  float s0 = s1 + 1;
  float s3 = s2 + 1;
  Vec3f L01 = -s1 * p[i0] + s0 * p[i1];
  Vec3f L12 = s2 * p[i1] + s1 * p[i2];
  Vec3f L23 = s3 * p[i2] - s2 * p[i3];
  Vec3f L012 = 0.5f * (s2 * L01 + s0 * L12);
  Vec3f L123 = 0.5f * (s3 * L12 + s1 * L23);
  C12 = s2 * L012 + s1 * L123;
  return C12;
}

Vec3u SplineGrid::GridIndex(Vec3f x) const {
  x -= _origin;
  Vec3u index;
  Vec3u size = _grid.GetSize();
  for (unsigned d = 0; d < 3; d++) {
    x[d] = std::max(x[d], 0.0f);
    index[d] = unsigned(x[d] / _voxSize[d]);
    index[d] = std::min(index[d], size[d] - 1);
  }
  return index;
}

void SplineGrid::AddPoint(const CurvePoint& pt) { 
  Vec3u index = GridIndex(pt.x);
  size_t listIndex = _grid(index[0], index[1], index[2]);
  if (listIndex == 0) {
    //create a new list
    listIndex = _points.size();
    _points.push_back(std::vector<CurvePoint>(1, pt));
    _grid(index[0], index[1],index[2]) = listIndex;
  } else {
    _points[listIndex].push_back(pt);
  }
}

void SplineGrid::AddPoints(const SplineMesh& splines, unsigned samplePerSeg) {
  std::vector<CurvePoint> samples = SampleCurve(splines, samplePerSeg);
  _points.resize(1);
  if (samples.size() == 0) {
    return;
  }
  BBox box(samples[0].x, samples[0].x);
  for (size_t i = 1; i < samples.size(); i++) {
    for (unsigned d = 0; d < 3; d++) {
      box.vmax[d] = std::max(box.vmax[d], samples[i].x[d]);
      box.vmin[d] = std::min(box.vmin[d], samples[i].x[d]);
    }
  }
  Vec3f margin = 0.5f * _voxSize;
  _origin = box.vmin - margin;
  //bounding box size plus margin
  Vec3f boxSize = box.vmax - box.vmin + 2.0f * margin;
  Vec3u gridSize(boxSize[0] / _voxSize[0], boxSize[1] / _voxSize[1], boxSize[2] / _voxSize[2]);
  gridSize += Vec3u(1);
  _grid.Allocate(gridSize, 0);
  for (unsigned i = 0; i < samples.size();i++) {
    AddPoint(samples[i]);
  }
}

std::vector<CurvePoint> SampleCurve(const SplineMesh& splines,
  unsigned ptPerSeg) {
  size_t numCurves = splines.Size();
  std::vector<CurvePoint> samples;
  ptPerSeg = std::max(1u, ptPerSeg);
  for (size_t i = 0; i < numCurves; i++) {
    const auto& spline = splines[i];
    for (size_t seg = 0; seg < spline.NumSegments(); seg++) {
      for (size_t j = 0; j < ptPerSeg; j++) {
        float s = seg + float(j) / ptPerSeg;
        CurvePoint pt(i,s,spline.Eval(s));
        samples.push_back(pt);
      }
    }
    //add last point of last segment
    float s = spline.NumSegments();
    samples.push_back(CurvePoint(i, s, spline.Eval(s)));
  }
  return samples;
}

/// @brief If nearby list already contains another point from 
/// the same curve in the same segment, update the min distance
/// to that segment. otherwise insert pt into nearby.
/// @param newDist 
/// @param pt 
/// @param nearby 
void UpdateDistToSegment(const Vec3f coord, float newDist, const CurvePoint& pt, 
  std::vector<CurvePoint> & nearby) {
  for (size_t i = 0; i < nearby.size(); i++) {
    if (nearby[i].curveId == pt.curveId && int(nearby[i].s) == int(pt.s)) {
      float oldDist = (coord - nearby[i].x).norm();
      if (newDist < oldDist) {
        nearby[i] = pt;
      }
      return;
    }
  }
  nearby.push_back(pt);
}

std::vector<CurvePoint> SplineGrid::FindNearbyPoints(Vec3f coord, float minDist,
                                                     float maxDist) const {
  maxDist = std::max(maxDist, 0.0f);
  Vec3i searchRadius(std::ceil(maxDist / _voxSize[0]), std::ceil(maxDist / _voxSize[1]),
                     std::ceil(maxDist / _voxSize[2]));
  std::vector<CurvePoint> nearby;
  Vec3u index = GridIndex(coord);
  int x0 = int(index[0]);
  int y0 = int(index[1]);
  int z0 = int(index[2]);
  Vec3u size = _grid.GetSize();
  
  for (int z = z0 - searchRadius[2]; z <= z0 + searchRadius[2]; z++) {
    if (z < 0 || z >= int(size[2])) {
      continue;
    }
    for (int y = y0 - searchRadius[1]; y <= y0 + searchRadius[2]; y++) {
      if (y < 0 || y >= int(size[1])) {
        continue;
      }
      for (int x = x0 - searchRadius[0]; x <= x0 + searchRadius[0]; x++) {
        if (x < 0 || x >= int(size[0])) {
          continue;
        }
        size_t listIndex = _grid(x, y, z);
        if (listIndex == 0) {
          //empty list
          continue;
        }
        const auto& pointList = _points[listIndex];
        for (size_t i = 0; i < pointList.size(); i++) {
          float dist = (coord - pointList[i].x).norm();
          if (dist < minDist) {
            nearby.push_back(pointList[i]);
          //  //early return when a sample point is very close.
            return nearby;
          }
          if (dist > maxDist) {
            continue;
          }
          UpdateDistToSegment(coord, dist, pointList[i], nearby);
        }
      }
    }
  }
  return nearby;
}

PointSplineDistance BetweenCurves(const Vec3f& coord, const std::vector<CurvePoint>& nearby,
                                  const SplineMesh& splines, float waxGap, float yarnRad) {
  PointSplineDistance pDist;
  if (nearby.empty()) {
    return pDist;
  }
  std::vector<float> dists(nearby.size());
  for (size_t i = 0; i < dists.size(); i++) {
    dists[i] = (nearby[i].x - coord).norm();
  }
  if (nearby.size() < 2) {
    pDist.dist = dists[0];
    pDist.betweenCurves = false;
    return pDist;
  }
  float minDist = dists[0];
  int closest = 0;
  for (size_t i = 1; i < nearby.size(); i++) {
    if (dists[i] < minDist) {
      minDist = dists[i];
      closest = i;
    }
  }
  pDist.dist = minDist;
  float distBound = minDist + waxGap;
  for (size_t i = 0; i < nearby.size(); i++) {
    if (i == closest) {
      continue;
    }
    // this point is on the boundary between two curves
    // or between two far away segments within the same curve.
    if (((nearby[closest].x - nearby[i].x).norm() > yarnRad) && dists[i] < distBound) {
      pDist.betweenCurves = true;
      break;
    }
  }
  return pDist;
}
