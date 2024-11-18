#ifndef ARRAY_2D_H
#define ARRAY_2D_H

#include <string>
#include <vector>

#include "Vec2.h"

//for color images
#include "Vec4.h"

/// a minimal class that interprets a 1D vector as
/// a 2D array
template <typename T>
class Array2D {
 public:
  Array2D() {
    size[0] = 0;
    size[1] = 0;
  }

  Array2D(unsigned width, unsigned height) { Allocate(width, height); }

  void Allocate(unsigned width, unsigned height) {
    size[0] = width;
    size[1] = height;
    data.resize(width * size_t(height));
  }

  void Allocate(Vec2u size) {
    Allocate(size[0], size[1]);
  }

  void Allocate(unsigned width, unsigned height, const T& initVal) {
    size[0] = width;
    size[1] = height;
    data.resize(width * size_t(height), initVal);
  }

  void Fill(const T& val) { std::fill(data.begin(), data.end(), val); }
  ///\param x column
  ///\param y row
  T& operator()(unsigned int x, unsigned int y) {
    return data[x + y * size[0]];
  }

  const T& operator()(unsigned int x, unsigned int y) const {
    return data[x + y * size[0]];
  }

  const Vec2u& GetSize() const { return size; }
  void SetSize(unsigned w, unsigned h) {
    size[0] = w;
    size[1] = h;
  }
  std::vector<T>& GetData() { return data; }
  const std::vector<T>& GetData() const { return data; }
  T* DataPtr() { return data.data(); }
  const T* DataPtr() const { return data.data(); }

  bool Empty()const { return data.size() == 0; }

  Array2D<T>& operator+=(const Array2D<T>& a1) {
    if (a1.GetSize() != GetSize()) {
      return *this;
    }
    for (size_t i = 0; i < data.size(); i++) {
      data[i] += a1.GetData()[i];
    }
    return *this;
  }

  Array2D<T>& operator*=(T rhs) {
    for (size_t i = 0; i < data.size(); i++) {
      data[i] *= rhs;
    }
    return *this;
  }

  unsigned Rows()const {
    return size[1];
  }
  unsigned Cols() const {
    return size[0];
  }

 private:
  Vec2u size;
  std::vector<T> data;
};

typedef Array2D<float> Array2Df;
typedef Array2D<unsigned char> Array2D8u;

typedef Array2D<unsigned char> Slice;
typedef Array2D<Vec4b>Array2D4b;
void flipY(Slice& slice);

/// scales the image using nearest neighbor
void scale(float scalex, float scaley, Slice& slice);

void addMargin(Array2D<uint8_t>* arr, const Vec2i& margin);


/// scalar multiplication
template <typename T>
Array2D<T> operator*(T s, const Array2D<T>& v) {
  Array2D<T> a(v.GetSize()[0],v.GetSize()[1]);
  for (size_t i = 0; i < a.GetData().size(); i++) {
    a.GetData()[i] = s * v.GetData()[i];
  }
  return a;
}

template <typename T>
Array2D<T> operator+(const Array2D<T>& a0, const Array2D<T>& a1) {
  Array2D<T> sum;
  if (a0.GetSize() != a1.GetSize()) {
    return sum;
  }
  for (size_t i = 0; i < a0.GetData().size(); i++) {
    sum.GetData()[i] = a0.GetData()[i] + a1.GetData()[i];
  }
  return sum;
}

#endif
