#pragma once

#include <string>
#include <vector>
#include "BBox.h"
#include "Vec3.h"

using CurvePatch = std::vector<std::vector<Vec3f> > ;

void BuildUnitPatch();
CurvePatch LoadCurvePatch(const std::string& filename);
void DrawCurvePatch(const CurvePatch& curvePatch);
void SaveCurvePatch(const std::string& file, const CurvePatch& patch);
void SaveCurvePatchObj(const std::string& objfile, const CurvePatch& patch);
CurvePatch MakeTilable(const CurvePatch& patch, const BBox& tileBox);
/// add more sample points between points further away than maxDist.
void FillInGaps(CurvePatch& patch, float maxDist);
void MakeTilablePatch();