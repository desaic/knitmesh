#pragma once
#include "Spline.h"
#include "Result.h"
#include <filesystem>
/// @brief http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html
/// @param path 
/// @return 
Result<SplineMesh> LoadBCC(const std::filesystem::path& path);
int SaveSplineTxt(const SplineMesh& splines, const std::filesystem::path& outPath);
Result<SplineMesh> LoadSplineTxt(const std::filesystem::path& inPath);