#include "SplineIO.h"

#include <fstream>
#include <filesystem>

struct BCCHeader {
  char sign[3];
  unsigned char byteCount;
  char curveType[2];
  char dimensions;
  char upDimension;
  uint64_t curveCount;
  uint64_t totalControlPointCount;
  char fileInfo[40];
};

int CheckHeader(const BCCHeader& header) {
  if (header.sign[0] != 'B') return -1;       // Invalid file signature
  if (header.sign[1] != 'C') return -1;       // Invalid file signature
  if (header.sign[2] != 'C') return -1;       // Invalid file signature
  if (header.byteCount != 0x44) return -2;    // Only supporting 4-byte integers and floats
  if (header.curveType[0] != 'C') return -3;  // Not a Catmull-Rom curve
  if (header.curveType[1] != '0') return -4;  // Not uniform parameterization
  if (header.dimensions != 3) return -5;      // Only curves in 3D
  return 0;
}

Result<SplineMesh> LoadBCC(const std::filesystem::path& path) {
  BCCHeader header;
  std::ifstream in(path, std::ifstream::binary);
  in.read((char*)(&header), sizeof(header));
  int errCode = CheckHeader(header);
  if (errCode < 0) {
    return Result<SplineMesh>(Error("invalid file header"));
  }
  SplineMesh splines;
  splines.Allocate(header.curveCount);
  for (uint64_t i = 0; i < header.curveCount; i++) {
    int curveControlPointCount;
    in.read((char*)(&curveControlPointCount), sizeof(int));
    splines[i].isLoop = curveControlPointCount < 0;
    if (curveControlPointCount < 0) curveControlPointCount = -curveControlPointCount;
    splines[i].Resize(curveControlPointCount);
    in.read((char*)splines[i].p.data(), 3 * sizeof(float) * curveControlPointCount);
  }
  return splines;
}

int SaveSplineTxt(const SplineMesh& splines, const std::filesystem::path& outPath) {
  std::ofstream out(outPath);
  if (!out.good()) {
    return -1;
  }
  out << splines.Size() << "\n";
  for (size_t i = 0; i < splines.Size(); i++) {
    out << int(splines[i].isLoop) << " " <<splines[i].NumControlPoints()<<"\n";
    for (size_t j = 0; j < splines[i].p.size(); j++) {
      const Vec3f& p = splines[i].p[j];
      out << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
  }
  return 0;
}

Result<SplineMesh> LoadSplineTxt(const std::filesystem::path& inPath) {
  std::ifstream in(inPath);
  if (!in.good()) {
    return Result<SplineMesh>(Errc::no_such_file_or_directory);
  }
  size_t numSplines = 0;
  in >> numSplines;
  const unsigned MAX_NUM_SPLINES = 100000000;
  const unsigned MAX_NUM_POINTS = 100000000;
  if (numSplines == 0 || numSplines > MAX_NUM_SPLINES) {
    return Result<SplineMesh>(Errc::file_too_large);
  }
  SplineMesh splines;
  splines.Allocate(numSplines);
  for (size_t i = 0; i < splines.Size(); i++) {
    int isLoop = 0;
    in >> isLoop;
    splines[i].isLoop = (isLoop > 0);
    int numControlPoints = 0;
    in>>numControlPoints;
    if (numControlPoints <= 0 || numControlPoints > MAX_NUM_POINTS) {
      return Result<SplineMesh>(Errc::file_too_large);
    }
    splines[i].Resize(numControlPoints);
    for (size_t j = 0; j < splines[i].p.size(); j++) {
      Vec3f& p = splines[i].p[j];
      in >> p[0] >> p[1] >> p[2];
    }
  }
  return splines;
}