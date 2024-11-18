#pragma once
#include "Array3D.h"
#include "Vec3.h"
#include <vector>

/// @brief C1 continuous Catmull-Rom spline
/// with cubic polynomials.
class SplineCR1 {
 public:
  SplineCR1() : p(4) {}
  void Resize(unsigned s) { p.resize(s); }
   /// @brief http://www.cemyuksel.com/research/catmullrom_param/
   /// @param s from 0 to n-1 where n is number of segments.
  Vec3f Eval(float s) const;
  
  /// @brief interpolate at s treating the curve as a loop.
  /// @param s 
  /// @return 
  Vec3f EvalLoop(float s) const;
  /// @brief 
  /// @return negative count if not enough control points. 
  int NumSegments() const;
  size_t NumControlPoints() const { return p.size(); }
  /// control points. must have at least 4 for open curves.
  /// number of segments = p-3. p0 p1 p2 p3 for segment 0.
  std::vector<Vec3f> p;
  bool isLoop = false;
};

class SplineMesh {
 public:
  void Allocate(size_t num) { splines.resize(num); }
  const SplineCR1& operator[](unsigned i) const { return splines[i]; }
  SplineCR1& operator[](unsigned i) { return splines[i]; }
  size_t Size() const { return splines.size(); }
  std::vector<SplineCR1> splines;
};

//point on a curve 
struct CurvePoint {
  CurvePoint(int curveIndex, float curveParam, Vec3f pos)
      : curveId(curveIndex), s(curveParam), x(pos) {}
  
  int curveId = 0;
  //segment id = floor(s).
  /// curve parameterization
  float s;
  /// position
  Vec3f x;
};

struct CurveSet {
  std::vector<std::vector<CurvePoint> > curves;

};

/// @brief acceleration for finding nearby spline points.
class SplineGrid {
 public:
  SplineGrid() : _points(1), _origin(0), _voxSize(1) {}

  void AddPoint(const CurvePoint& pt);

  // sample points from the spline and add them to the grid
  void AddPoints(const SplineMesh& splines, unsigned samplePerSeg);

  /// finds one closest point for each distince segment.
  /// if closest point is closer than minDist, search will stop early.
  std::vector<CurvePoint> FindNearbyPoints(Vec3f x, float minDist, float maxDist) const;

  Vec3u GridIndex(Vec3f x) const;

  // index into list of list of curve points.
  // 0 is reserved for empty list.
  Array3D<size_t> _grid;
  std::vector<std::vector<CurvePoint> > _points;
  // origin of the grid.
  Vec3f _origin;
  /// cell size of grid.
  Vec3f _voxSize;
};

std::vector<CurvePoint> SampleCurve(const SplineMesh& splines, unsigned ptPerSeg);

struct PointSplineDistance {
  /// @brief distance to closest point on splines
  float dist = 1e10f;
  bool betweenCurves = false;
};

/// @brief is this point on the boundary between two curves
/// or between two far away segments within the same curve.
/// @param coord
/// @param nearby
/// @param splines original spline mesh where the nearby points are from.
/// @param waxGap
/// @param yarnRad
/// @return 1 if coord is between curves or segments. 0 otherwise.
PointSplineDistance BetweenCurves(const Vec3f& coord, const std::vector<CurvePoint>& nearby,
                                  const SplineMesh& splines, float waxGap, float yarnRad);