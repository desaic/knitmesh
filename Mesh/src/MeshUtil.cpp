#include "MeshUtil.h"

#include "Array3D.h"
#include "TrigMesh.h"
#include "MarchingCubes.h"

#include <iostream>
#include <fstream>
#include <thread>
#include <unordered_map>
struct QuadMesh {
    std::vector<std::array<unsigned,4> > q;
    std::vector<Vec3f> v;
};
inline uint64_t gridToLinearIdx(int ix, int iy, int iz, const int* gridSize) {
  return ix * (uint64_t)gridSize[1] * gridSize[2] + iy * (uint64_t)gridSize[2] +
         iz;
}
void saveCalibObj(const float * table, const int * tsize,
                  std::string filename, int subsample)
{
    std::ofstream out (filename);
    for(int ix = 0; ix<tsize[0]; ix+=subsample){
        for(int iy = 0; iy<tsize[1]; iy+=subsample){
            for(int iz = 0; iz<tsize[2]; iz+=subsample){
                out<<"v";
                int tableIdx = gridToLinearIdx(ix, iy, iz, tsize);
                for(int d = 0; d<3; d++){
                    out<<" "<<table[tableIdx*3 + d];
                }
                out<<"\n";
            }
        }
    }
    out.close();
}

void savePointsObj(const std::string & fileName, const std::vector<std::array<float, 3 > > & points)
{
	std::ofstream out(fileName);
	if (!out.good()) {
		return;
	}
	for (size_t i = 0; i < points.size(); i++) {
		out << "v " << points[i][0] <<" "<< points[i][1] << " " << points[i][2] << "\n";
	}
	out.close();
}

void addQuad(QuadMesh& m, const std::vector<Vec3f>& v) {
    std::array<unsigned,4> q;
    for (int i = 0; i < 4; i++) {
        q[i] = i + m.v.size();
    }
    m.q.push_back(q);
    m.v.insert(m.v.end(), v.begin(), v.end());
}

void growRect(const std::vector <uint8_t> & plane, int x0, int y0, int & x1, int & y1, int nx, int ny) {
    std::vector<uint8_t> colFlag(ny, 1);
    int y = y0;
    for (; y < ny; y++) {
        if (plane[x0*ny + y] == 0) {
            break;
        }
    }
    y1 = y - 1;

    bool fullRow = true;
    for (int x = x0 + 1; x < nx; x++) {
        for (int y = y0; y <= y1; y++) {
            if (plane[y + x * ny] == 0) {
                fullRow = false;
                break;
            }
        }
        if (!fullRow) {
            x1 = x - 1;
            break;
        }
    }
    if (fullRow) {
        x1 = nx - 1;
    }
}

void addLayer(QuadMesh * mpt, int layer0, int layer1, int faceIdx, const Array3D8u * volpt, int mat, const Vec3f * dxpt) {
    QuadMesh & m = (*mpt);
    const Array3D8u&vol = (*volpt);
    const Vec3f & dx = (*dxpt);
    Vec3u gridSize = vol.GetSize();
    //-1 for negative side and 1 for positive side face
    int faceOffset = faceIdx % 2;
    int sign = 2 * (faceOffset) - 1;
    int zAxis = faceIdx / 2;
    int xAxis = ( zAxis + 1 ) % 3;
    int yAxis = ( zAxis + 2 ) % 3;
    int nx = gridSize[xAxis];
    int ny = gridSize[yAxis];
    int nz = gridSize[zAxis];
    int volIdx[3];
    int nbrIdx[3];

    std::vector <uint8_t> plane(nx*ny, 0);
    //copy visible voxels in this layer to 2D plane
    for (int layer = layer0; layer < layer1; layer++) {
        volIdx[zAxis] = layer;
        nbrIdx[zAxis] = layer + sign;

        for (int x = 0; x < nx; x++) {
            volIdx[xAxis] = x;
            nbrIdx[xAxis] = x;
            for (int y = 0; y < ny; y++) {
                volIdx[yAxis] = y;
                nbrIdx[yAxis] = y;
                uint8_t val = vol(volIdx[0], volIdx[1], volIdx[2]);
                if (val != mat) {
                    continue;
                }
                if (nbrIdx[zAxis] >= 0 && nbrIdx[zAxis] < nz) {
                    int nbrVal = vol(nbrIdx[0], nbrIdx[1], nbrIdx[2]);
                    if (nbrVal == mat) {
                        continue;
                    }
                }
                plane[x * ny + y] = mat;
            }
        }
        for (int px = 0; px < nx; px++) {
            for (int py = 0; py < ny; py++) {
                if (plane[px*ny + py] != mat) {
                    continue;
                }
                int x0 = px;
                int y0 = py;
                int x1 = px;
                int y1 = py;
                growRect(plane, x0, y0, x1, y1, nx, ny);
                for (int rx = x0; rx <= x1; rx++) {
                    for (int ry = y0; ry <= y1; ry++) {
                        plane[rx*ny + ry] = 0;
                    }
                }
                std::vector<Vec3f> q(4);
                for (int qi = 0; qi < 4; qi++) {
                    q[qi][zAxis] = layer + faceOffset;
                }
                q[0][xAxis] = x0;
                q[0][yAxis] = y0;

                q[2][xAxis] = x1 + 1;
                q[2][yAxis] = y1 + 1;
                if (sign == 1) {
                    q[1][xAxis] = x1 + 1;
                    q[1][yAxis] = y0;

                    q[3][xAxis] = x0;
                    q[3][yAxis] = y1 + 1;
                }
                else {
                    q[3][xAxis] = x1 + 1;
                    q[3][yAxis] = y0;

                    q[1][xAxis] = x0;
                    q[1][yAxis] = y1 + 1;
                }
                for (int qi = 0; qi < 4; qi++) {
                    for (int dim = 0; dim < 3; dim++) {
                        q[qi][dim] *= dx[dim];
                    }
                }
                addQuad(m, q);
            }
        }
    }
}

void appendQuads(QuadMesh & m, const QuadMesh & src) {
    size_t nV = m.v.size();
    size_t nq = m.q.size();
    m.v.insert(m.v.end(), src.v.begin(), src.v.end());
    m.q.insert(m.q.end(), src.q.begin(), src.q.end());
    for (size_t i = 0; i < src.q.size(); i++) {
        for (int j = 0; j < 4; j++) {
            m.q[i + nq][j] += nV;
        }
    }
}

void addFaces(QuadMesh& m, int faceIdx, const Array3D8u& vol, int mat,
              const Vec3f& dx) {
    int axis = faceIdx / 2;
    int nLayer = vol.GetSize()[axis];
    int nThread = 10;
    std::vector<QuadMesh> localM(nThread);
    std::vector<std::thread> threads(nThread);
    int layerPerThread = nLayer / nThread + ((nLayer % nThread) != 0);
    for (int t = 0; t < nThread; t++) {
        int layer0 = t * layerPerThread;
        int layer1 = std::min((t+1)*layerPerThread, nLayer);
        threads[t] = std::thread(addLayer, &localM[t], layer0, layer1, faceIdx, &vol, mat, &dx);
    }
    for (int t = 0; t < nThread; t++) {
        if (threads[t].joinable()) {
            threads[t].join();
        }
        appendQuads(m, localM[t]);
    }
}

void save_obj(const QuadMesh & m, const std::string & filename) {
    std::ofstream out(filename);
    if (!out.good()) {
        std::cout << "cannot open output file" << filename << "\n";
        return;
    }

    std::string vTok("v");
    std::string fTok("f");
    char bslash = '/';
    std::string tok;
    const std::vector<Vec3f>* vert = &m.v;
    for (size_t ii = 0; ii<vert->size(); ii++) {
        out << vTok << " " << (*vert)[ii][0] << " " << (*vert)[ii][1] << " " << (*vert)[ii][2] << "\n";
    }
    for (size_t ii = 0; ii<m.q.size(); ii++) {
        out << fTok << " " << m.q[ii][0] + 1 << " " <<
            m.q[ii][1] + 1 << " " << m.q[ii][2] + 1 << " "<<m.q[ii][3] + 1 << "\n";
    }
    out << "#end\n";
    out.close();
}

void SaveVolAsObjMesh(std::string outfile, const Array3D8u& vol, float* voxRes, int mat) {
  Vec3f origin(0);
  SaveVolAsObjMesh(outfile, vol, voxRes,
    (float*)(&origin), mat);
}

void SaveVolAsObjMesh(std::string outfile, const Array3D8u& vol, float* voxRes,
                      float* origin, int mat) {
  QuadMesh mout;
  // convert zres to mm.
  Vec3f dx = {voxRes[0], voxRes[1], voxRes[2]};
  const int N_FACE = 6;
  for (int f = 0; f < N_FACE; f++) {
    addFaces(mout, f, vol, mat, dx);
  }
  for (size_t i = 0; i < mout.v.size(); i ++) {
    mout.v[i] += Vec3f(origin[0],origin[1],origin[2]);
  }
  save_obj(mout, outfile.c_str());
}

uint8_t LinearIndex(int x, int y, int z, std::vector<uint8_t>& vol, std::array<unsigned, 3>& gsize) {
  return vol[gsize[0] * (gsize[1] * z + y) + x];
}

void MarchOneCube(int x, int y, int z, TrigMesh& m, std::vector<uint8_t>& vol, std::array<unsigned, 3>& gsize) {
    GridCell cell;
    cell.p[0] = Vec3f(0, 0, 0);
    float h = 0.6f;
    cell.p[0][0] += x * h;
    cell.p[0][1] += y * h;
    cell.p[0][2] += z * h;

    for (unsigned i = 1; i < GridCell::NUM_PT; i++) {
      cell.p[i] = cell.p[0];
    }
    cell.p[1][1] += h;

    cell.p[2][0] += h;
    cell.p[2][1] += h;

    cell.p[3][0] += h;

    cell.p[4][2] += h;

    cell.p[5][1] += h;
    cell.p[5][2] += h;

    cell.p[6][0] += h;
    cell.p[6][1] += h;
    cell.p[6][2] += h;

    cell.p[7][0] += h;
    cell.p[7][2] += h;

    cell.val[0] = LinearIndex(x,y,z,vol,gsize);
    cell.val[1] = LinearIndex(x,y+1,z,vol,gsize);
    cell.val[2] = LinearIndex(x+1,y+1,z,vol,gsize);
    cell.val[3] = LinearIndex(x+1,y,z,vol,gsize);
    cell.val[4] = LinearIndex(x,y,z+1,vol,gsize);
    cell.val[5] = LinearIndex(x,y+1,z+1,vol,gsize);
    cell.val[6] = LinearIndex(x+1,y+1,z+1,vol,gsize);
    cell.val[7] = LinearIndex(x+1,y,z+1,vol,gsize);
    MarchCube(cell,0.5, &m);
}

TrigMesh MarchCubes(std::vector<uint8_t>& vol, std::array<unsigned, 3>& gsize) {
    TrigMesh m;
    for (unsigned z = 0; z < gsize[2]/2; z++) {
      for (unsigned y = 0; y < gsize[1]/2; y++) {
        for (unsigned x = 0; x < gsize[0]/2; x++) {
          MarchOneCube(x, y, z, m, vol, gsize);
        }
      }
    }
    return m;
}

void ComputeUVTangents(Vec3f& tx, Vec3f& ty, const Vec3f& e1, const Vec3f& e2, const Vec2f& uve1,
                       const Vec2f& uve2) {
  float l1, l2;
  float det = uve1[0] * uve2[1] - uve1[1] * uve2[0];
  if (det == 0) {
    // uv triangle is degenerate.
    tx = e1;
    ty = e2;
    return;
  }
  l1 = uve2[1] / det;
  l2 = -uve1[1] / det;
  tx = l1 * e1 + l2 * e2;

  l1 = -uve2[0] / det;
  l2 = uve1[0] / det;
  ty = l1 * e1 + l2 * e2;
}

void SamplePointsOneTrig(unsigned tIdx, const TrigMesh& m, std::vector<SurfacePoint>& points,
                         float spacing) {
  Vec3f v0 = m.GetTriangleVertex(tIdx, 0);
  Vec3f v1 = m.GetTriangleVertex(tIdx, 1);
  Vec3f v2 = m.GetTriangleVertex(tIdx, 2);
  Vec3f e1 = v1 - v0;
  Vec3f e2 = v2 - v0;
  Vec3f cross = e1.cross(e2);
  float norm = cross.norm();
  float area = 0.5f * norm;
  Vec3f n;

  if (norm > 0) {
    n = (1.0f / norm) * cross;
  } else {
    // degenerate triangle.
    return;
  }

  Vec2f uv0 = m.GetTriangleUV(tIdx, 0);
  Vec2f uv1 = m.GetTriangleUV(tIdx, 1);
  Vec2f uv2 = m.GetTriangleUV(tIdx, 2);
  Vec2f uve1 = uv1 - uv0;
  Vec2f uve2 = uv2 - uv0;
  Vec3f tx, ty;
  // same within a triangle.
  ComputeUVTangents(tx, ty, e1, e2, uve1, uve2);
  const float oneThird = 1.0f / 3.0f;
  SurfacePoint sp;
  sp.v = oneThird * (v0 + v1 + v2);
  sp.n = n;
  sp.uv = oneThird * (uv0 + uv1 + uv2);
  sp.tx = tx;
  sp.ty = ty;
  points.push_back(sp);
  if (area < spacing * spacing) {
    // small triangle. one center point is enough.
    return;
  }

  float e1Len = e1.norm();
  Vec3f x = (1.0f / e1Len) * e1;
  float e2Len = e2.norm();
  Vec3f y = (1.0f / e2Len) * e2;

  Vec3f e3 = v2 - v1;

  Vec3f e3mid = 0.5f * (v1 + v2);
  Vec3f centerLine = v0 - e3mid;
  float e3len = e3.norm();
  Vec3f e3Unit = (1.0f / e3len) * e3;
  Vec3f projection = centerLine.dot(e3Unit) * e3Unit;
  Vec3f hVec = centerLine - projection;
  float h = hVec.norm();
  int numHSteps = h / spacing + 1;
  for (int hi = 0; hi < numHSteps; hi++) {
    // between 0 to 1
    float y = (hi + 0.5f) / numHSteps;
    float xlen = y * e3len;
    int numXSteps = int(xlen / spacing);
    if (numXSteps == 0) {
      // very narrow area, nothing to sample.
      continue;
    }
    for (int xi = 0; xi < numXSteps; xi++) {
      float x = (xi + 0.5f) / numXSteps;
      Vec3f point = v0 + (1 - y) * (x * e1 + (1 - x) * e2);
      SurfacePoint sp;
      sp.v = point;
      sp.uv = (1 - y) * x * uv1 + (1 - y) * (1 - x) * uv2 + y * uv0;
      sp.n = n;
      sp.tx = tx;
      sp.ty = ty;
      points.push_back(sp);
    }
  }
}

void SamplePoints(const TrigMesh& m, std::vector<SurfacePoint>& points,
                  float spacing) {
    size_t numTrigs = m.t.size() / 3;
    for (size_t ti = 0; ti < numTrigs; ti++) {
      SamplePointsOneTrig(ti, m, points, spacing);
    }
}

float SurfaceArea(const TrigMesh& m) {
  size_t numTrigs = m.t.size() / 3;
  double area = 0;
  for (size_t ti = 0; ti < numTrigs; ti++) {
    Vec3f v0 = m.GetTriangleVertex(ti, 0);
    Vec3f v1 = m.GetTriangleVertex(ti, 1);
    Vec3f v2 = m.GetTriangleVertex(ti, 2);
    v1 = v1 - v0;
    v2 = v2 - v0;
    v1 = v1.cross(v2);
    float trigArea = v1.norm();
    area += trigArea;
  }
  return 0.5f * float(area);
}

float UVArea(const TrigMesh& m) {
  size_t numTrigs = m.t.size() / 3;
  double area = 0;
  for (size_t ti = 0; ti < numTrigs; ti++) {
    Vec2f uv0 = m.GetTriangleUV(ti, 0);
    Vec2f uv1 = m.GetTriangleUV(ti, 1);
    Vec2f uv2 = m.GetTriangleUV(ti, 2);
    uv1 = uv1 - uv0;
    uv2 = uv2 - uv0;

    float trigArea = uv1[0] * uv2[1] - uv1[1] * uv2[0];
    area += std::abs(trigArea);
  }
  return 0.5f * float(area);
}

static float MergeEps = 1e-3f;
class Vec3fKey {
 public:
  Vec3fKey() {}
  Vec3fKey(float x, float y, float z) : _x(x), _y(y), _z(z) {}
  float _x = 0;
  float _y = 0;
  float _z = 0;
  bool operator==(const Vec3fKey& other) const {
    return std::abs(_x - other._x) <= MergeEps && std::abs(_y - other._y) <= MergeEps &&
           std::abs(_z - other._z) <= MergeEps;
  }
};
struct Vec3fHash {
  size_t operator()(const Vec3fKey& k) const {
    return size_t((92821 * 92831) * int64_t(k._x / MergeEps) + 92831 * int64_t(k._y / MergeEps) +
                  int64_t(k._z / MergeEps));
  }
};

//remove all squished triangles
bool DegenerateTrig(uint32_t* t) { return t[0] == t[1] || t[0] == t[2] || t[1] == t[2]; }

void MergeCloseVertices(TrigMesh& m) {
  if (m.v.size() == 0 || m.t.size() == 0) {
    return;
  }
  /// maps from a vertex to its new index.
  std::vector<float> v_out;
  std::vector<unsigned> t_out;
  std::unordered_map<Vec3fKey, size_t, Vec3fHash> vertMap;
  size_t nV = m.v.size() / 3;
  std::vector<uint32_t> newVertIdx(nV, 0);
  size_t vCount = 0;
  for (size_t i = 0; i < m.v.size(); i += 3) {
    Vec3fKey vi(m.v[i], m.v[i + 1], m.v[i + 2]);
    auto it = vertMap.find(vi);
    size_t newIdx = 0;
    if (it == vertMap.end()) {
      vertMap[vi] = vCount;
      newIdx = vCount;
      v_out.push_back(m.v[i]);
      v_out.push_back(m.v[i + 1]);
      v_out.push_back(m.v[i + 2]);
      vCount++;
    } else {
      newIdx = it->second;
    }
    newVertIdx[i / 3] = newIdx;
  }
  for (size_t i = 0; i < m.t.size(); i+=3) {
    unsigned vertIdx[3];
    vertIdx[0] = newVertIdx[m.t[i]];
    vertIdx[1] = newVertIdx[m.t[i + 1]];
    vertIdx[2] = newVertIdx[m.t[i + 2]];
    if(DegenerateTrig(vertIdx)){
      continue;
    }
    t_out.push_back(vertIdx[0]);
    t_out.push_back(vertIdx[1]);
    t_out.push_back(vertIdx[2]);
  }
  if (v_out.size() != m.v.size()) {
    m.v = v_out;
    m.t = t_out;
    m.ComputeTrigNormals();
  }
}
