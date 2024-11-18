#include "ObjLoader.h"

#include <fstream>
#include <iostream>
#include <sstream>
int LoadOBJFile(std::string filename, std::vector<uint32_t>& t,
                std::vector<float>& v, std::vector<Vec2f>& uv) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    return -1;
  }
  t.clear();
  v.clear();
  uv.clear();
  // vertex list in uv plane.
  std::vector<Vec2f> uv_verts;
  // uv vertex index for each triangle
  std::vector<int> uv_t;
  std::string curline;

  while (std::getline(file, curline)) {
    auto iss = std::istringstream(curline);
    std::string first_token;
    iss >> first_token;
    if (first_token == "v") {
      for (int i = 0; i < 3; i++) {
        float vf;
        iss >> vf;
        v.push_back(vf);
      }
    } else if (first_token == "f") {
      for (int i = 0; i < 3; i++) {
        std::string idxs;
        iss >> idxs;
        std::istringstream idxss(idxs);

        std::string idx;
        std::getline(idxss, idx, '/');
        t.push_back(std::stoi(idx) - 1);

        // uv indices
        idx.clear();
        std::getline(idxss, idx, '/');
        if (!idx.empty()) {
          uv_t.push_back(std::stoi(idx) - 1);
        }
      }
    } else if (first_token == "vt") {
      Vec2f vt;
      iss >> vt[0] >> vt[1];
      uv_verts.push_back(vt);
    }
  }

  // order uvs properly
  for (int uv_idx : uv_t) {
    uv.push_back(uv_verts[uv_idx]);
  }
  return 0;
}