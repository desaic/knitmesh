#include "PolyMesh.h"

#include <fstream>
#include <sstream>

PolyMesh::PolyMesh() {}

void PolyMesh::Clear() {
  v.clear();
  f.clear();
}

void PolyMesh::Scale(float s) {
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = s * v[i];
  }
}

void PolyMesh::Translate(float dx, float dy, float dz) {
  for (size_t i = 0; i < v.size(); i++) {
    v[i][0] += dx;
    v[i][1] += dy;
    v[i][2] += dz;
  }
}

void PolyMesh::Append(const PolyMesh& m) {
  size_t offset = v.size();
  size_t ot = f.size();
  v.insert(v.end(), m.v.begin(), m.v.end());
  f.insert(f.end(), m.f.begin(), m.f.end());
  for (size_t ii = ot; ii < f.size(); ii++) {
    for (size_t jj = 0; jj < f[ii].n; jj++) {
      f[ii].v[jj] += (TIndex)(offset);
    }
  }
}

int PolyMesh::LoadObj(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.good()) {
    return -1;
  }
  Clear();
  // not very efficient without buffering data in memory.
  std::string line;
  while (std::getline(in, line)) {
    std::string dummy;
    if (line[0] == 'v') {
      Vec3f vin;
      std::stringstream ss(line);
      ss >> dummy >> vin[0] >> vin[1] >> vin[2];
      v.push_back(vin);
    } else if (line[0] == 'f') {
      std::stringstream ss(line);
      std::vector<unsigned> indices;
      ss >> dummy;
      for (unsigned i = 0; i < PolyFace::MAX_N; i++) {
        unsigned index;
        if (!(ss >> index)) {
          break;
        }
        // obj index is 1-based.
        indices.push_back(index - 1);
      }
      if (indices.size() >= 2 && indices.size() <= PolyFace::MAX_N) {
        f.push_back(PolyFace(indices));
      }
    }
  }
  return 0;
}

void PolyMesh::SaveObj(const std::string& filename) {}

PolyMesh SubSet(const PolyMesh& m, const std::vector<size_t>& faces) {
  PolyMesh sub;
  return sub;
}

PolyMesh MakeQuadCube(const Vec3f& mn, const Vec3f& mx) {
  PolyMesh cube;
  return cube;

}