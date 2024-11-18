#pragma once

#include <string>
#include <vector>

#include "Vec2.h"
#include "Vec3.h"

/// <summary>
/// 
/// </summary>
/// <param name="filename"></param>
/// <param name="t">triangle indices</param>
/// <param name="v">vertex coordinates</param>
/// <param name="uv">uv coordinates</param>
/// <returns>0 on success. -1 file not found.</returns>
int LoadOBJFile(std::string filename, std::vector<uint32_t>& t,
                std::vector<float>& v, std::vector<Vec2f>& uv);