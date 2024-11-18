#include "Array2D.h"
#include <algorithm>

void flipY(Slice& slice) {
  Vec2u size = slice.GetSize();
  unsigned halfYSize = size[1] / 2;
  std::vector<unsigned char> tmpRow(size[0]);
  for (unsigned y = 0; y < halfYSize; y++) {
    unsigned dsty = size[1] - y - 1;
    for (unsigned x = 0; x < size[0]; x++) {
      tmpRow[x] = slice(x, y);
      slice(x, y) = slice(x, dsty);
    }
    for (unsigned x = 0; x < size[0]; x++) {
      slice(x, dsty) = tmpRow[x];
    }
  }
}

/// scales the image stored inside this object.
void scale(float scalex, float scaley, Slice& slice) {
  float eps = 0.01f;
  float maxScale = 100.0f;
  scalex = std::clamp(scalex, eps, maxScale);
  scaley = std::clamp(scaley, eps, maxScale);
  Vec2u newSize;
  const Vec2u& size = slice.GetSize();
  newSize[0] = unsigned(size[0] * scalex);
  newSize[1] = unsigned(size[1] * scaley);
  if (newSize[0] == 0 || newSize[1] == 0) {
    return;
  }
  Slice newData;
  newData.Allocate(newSize[0], newSize[1]);
  for (unsigned dsty = 0; dsty < newSize[1]; dsty++) {
    unsigned srcy = unsigned(((dsty + 0.5) / newSize[1]) * size[1]);
    if (srcy >= size[1]) {
      continue;
    }
    for (unsigned dstx = 0; dstx < newSize[0]; dstx++) {
      unsigned srcx = unsigned(((dstx + 0.5) / newSize[0]) * size[0]);
      if (srcx >= size[0]) {
        continue;
      }
      newData(dstx, dsty) = slice(srcx, srcy);
    }
  }
  slice = std::move(newData);
}

void addMargin(Array2D<unsigned char>* arr, const Vec2i& margin) {
  Vec2u oldSize = arr->GetSize();
  Vec2u newSize = oldSize;
  newSize[0] += 2 * margin[0];
  newSize[1] += 2 * margin[1];
  Slice newSlice(newSize[0], newSize[1]);
  newSlice.Fill(0);
  for (unsigned row = 0; row < oldSize[1]; row++) {
    unsigned newRow = row + unsigned(margin[1]);
    size_t startIdx = size_t(newRow) * newSize[0] + unsigned(margin[0]);
    uint8_t* newRowPtr = &newSlice.DataPtr()[startIdx];
    uint8_t* oldRowPtr = &arr->DataPtr()[size_t(row) * oldSize[0]];
    std::memcpy(newRowPtr, oldRowPtr, oldSize[0]);
  }
  (*arr) = std::move(newSlice);
}
