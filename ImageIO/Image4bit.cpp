#include "Image4bit.h"

#include <algorithm>
#include <utility>

void Image4bit::ExtendLeft(unsigned dist) {
  Vec2i oldSize = _size;
  ImageData oldData = std::move(_data);
  _size = Vec2i({oldSize[0] + int(dist), oldSize[1]});
  _size[0] = _size[0] + ((8 - _size[0] % 8) % 8);
  _data.Alloc(_size[0], _size[1]);
  DualPixel dp;
  std::fill(_data.data.begin(), _data.data.end(), dp);
  for (size_t row = 0; row < oldSize[1]; row++) {
    size_t colOffset = _data.scanLine - oldData.scanLine;
    memcpy(&(_data.data[row * _data.scanLine + colOffset]), &(oldData.data[row * oldData.scanLine]),
           oldData.scanLine);
  }
}
