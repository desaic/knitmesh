#include <string>

#include "Array2D.h"
#include "Vec3.h"

Array2Df dlmread(const std::string& file);

void dlmwrite(const std::string& outfile, const Array2Df& mat);

void Add(std::vector<Vec3f>& dst, const std::vector<Vec3f>& src);

void Add(std::vector<float>& dst, const std::vector<float>& src);
void AddTimes(Array2Df& dst, const Array2Df& src, float c);
void AddTimes(std::vector<Array2Df>& dst, const std::vector<Array2Df>& src, float c);
void AddTimes(std::vector<Vec3f>& dst, const std::vector<double>& src, float c);

void Scale(Array2Df& dst, float s);
void Scale(std::vector<Array2Df>& dst, float s);

std::vector<Vec3f> operator*(float c, const std::vector<Vec3f>& v);

void Fix(std::vector<Vec3f>& dx, const std::vector<bool> fixedDOF);

Vec3f MaxAbs(const std::vector<Vec3f>& v);

/// <summary>
/// Matrix vector product
/// </summary>
void MVProd(const Array2Df& M, const float* v, unsigned vsize, float* out,
  unsigned outSize);

/// <summary>
/// M^T*v
/// </summary>
void MTVProd(const Array2Df& M, const float* v, unsigned vsize, float* out,
  unsigned outSize);

void ptwiseProd(const std::vector<double>& a, const std::vector<double>& b,
                std::vector<double>& c);

void sub(std::vector<double>& src, const std::vector<double>& b);

/// <summary>
/// subtract a constant from src
/// </summary>
void AddConstant(std::vector<float>& src, float b);

float Sum(const std::vector<float>& v);

std::vector<double> ToDouble(const std::vector<float>& a);

double Linf(const std::vector<double>& R);

double L2(const std::vector<double>& R);