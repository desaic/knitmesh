#pragma once

#include <string>
#include <vector>
#include "BBox.h"
#include "Vec3.h"

struct EndPoint {
  int curveId = 0;
  int pointId = 0;
  int boundary = 0;
  Vec3f pos;
};

class CurvePatch {
public:
  
  CurvePatch() {}
  CurvePatch(size_t s) :c(s) {}
  //curves
  std::vector<std::vector<Vec3f> > c;
  //2 endpoint labels for each curve.
  std::vector<EndPoint> ends;
  size_t size()const {
    return c.size();
  }
  void resize(size_t s) {
    c.resize(s);
  }

  std::vector<Vec3f>& operator[](unsigned index) {
    return c[index];
  }
  const std::vector<Vec3f>& operator[](unsigned index) const {
    return c[index];
  }
};

void BuildUnitPatch();
CurvePatch LoadCurvePatch(const std::string& filename);
void DrawCurvePatch(const CurvePatch& curvePatch);
void SaveCurvePatch(const std::string& file, const CurvePatch& patch);
void SaveCurvePatchObj(const std::string& objfile, const CurvePatch& patch);
CurvePatch MakeTilable(const CurvePatch& patch, const BBox& tileBox);
/// add more sample points between points further away than maxDist.
void FillInGaps(CurvePatch& patch, float maxDist);
void MakeTilablePatch();
std::vector<EndPoint> FindEndPointSides(const CurvePatch& patch, const BBox& box);

// rotate curve patch clockwise by 90 degs and also relabel the endpoints if 
// they are initialized.
CurvePatch RotateCW(const CurvePatch& patch, int quarterTurns);

class PatchModifier{
public:
  enum class OP {REPLACE = 0,  // replaces the points at an endpoint. 
    CONNECT //connect two end points by adding intermediate points
  };
  // end points are indexed by curves.
  // end points 0 1 for start and end of curve 0 etc.
  struct ModData{
    OP op = OP::REPLACE;
    int endPoint = 0;
    // for new connection within the same patch
    int endPointConn = -1;
    std::vector<Vec3f> p;
  };
  // apply modifier only when both myEdge and neighborEdge are matched.
  // undefined behavior when two modifiers with the same "myEdge" applied
  // to a patch.
  int myEdge = 0;
  int neighborEdge = -1;
  void Apply(CurvePatch& patch);
  std::vector<ModData> mods;
};

/// @brief 
/// @param curves in reference cube [-1,1].
/// @return 
std::vector<PatchModifier> GeneratePatchModifiers(const CurvePatch& curves);