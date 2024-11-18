#pragma once

#include "Array2D.h"
#include "Vec4.h"
#include <cmath>

bool Inbound(const Vec2u& size, int x, int y);

void DrawLine(Array2D4b& img, Vec2f p0, Vec2f p1, Vec4b color);

void ConvertChan3To4(const Array2D8u& rgb, Array2D8u& rgba);

void RGBToGrey(uint32_t& c);

void MakeCheckerPatternRGBA(Array2D8u& out);