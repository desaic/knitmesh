#include "ImageUtils.h"
#include "Vec4.h"

bool Inbound(const Vec2u& size, int x, int y) {
  return x >= 0 && x < size[0] && y >= 0 && y < size[1];
}

void DrawLine(Array2D4b& img, Vec2f p0, Vec2f p1, Vec4b color)
{
  int x0 = int(p0[0]);
  int x1 = int(p1[0]);
  int dx = std::abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int y0 = int(p0[1]);
  int y1 = int(p1[1]);
  int dy = -std::abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;
  Vec2u size = img.GetSize();
  while (1) {
    if (Inbound(size, x0, y0)) {
      img(x0, y0) = color;
    }
    if (x0 == x1 && y0 == y1) {
      break;
    }
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1) {
        break;
      }
      error += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      if (y0 == y1) {
        break;
      }
      error += dx;
      y0 += sy;
    }
  }
}

void DrawLine(Array2D8u& img, Vec2f p0, Vec2f p1, uint8_t color) {
  int x0 = int(p0[0]);
  int x1 = int(p1[0]);
  int dx = std::abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int y0 = int(p0[1]);
  int y1 = int(p1[1]);
  int dy = -std::abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;
  Vec2u size = img.GetSize();
  while (1) {
    if (Inbound(size, x0, y0)) {
      img(x0, y0) = color;
    }
    if (x0 == x1 && y0 == y1) {
      break;
    }
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1) {
        break;
      }
      error += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      if (y0 == y1) {
        break;
      }
      error += dx;
      y0 += sy;
    }
  }
}

void ConvertChan3To4(const Array2D8u& rgb, Array2D8u& rgba) {
  Vec2u inSize = rgb.GetSize();
  unsigned realWidth = inSize[0] / 3;
  rgba.Allocate(4 * realWidth, inSize[1]);
  size_t numPix = realWidth * inSize[1];
  const auto& inData = rgb.GetData();
  auto& outData = rgba.GetData();
  for (size_t i = 0; i < numPix; i++) {
    outData[4 * i] = inData[3 * i];
    outData[4 * i + 1] = inData[3 * i + 1];
    outData[4 * i + 2] = inData[3 * i + 2];
    outData[4 * i + 3] = 255;
  }
}

void RGBToGrey(uint32_t& c) {
  uint8_t red = (c >> 16) & 0xff;
  uint8_t green = (c >> 8) & 0xff;
  uint8_t blue = (c) & 0xff;
  float grey = 0.299 * red + 0.587 * green + 0.114 * blue;
  if (grey > 255) {
    grey = 255;
  }
  c=uint32_t(grey);
}


void MakeCheckerPatternRGBA(Array2D8u& out) {
  unsigned numCells = 4;
  unsigned cellSize = 4;
  unsigned width = numCells * cellSize;
  out.Allocate(width * 4, width);
  out.Fill(200);
  for (size_t i = 0; i < width; i++) {
    for (size_t j = 0; j < width; j++) {
      if ((i / cellSize + j / cellSize) % 2 == 0) {
        continue;
      }
      for (unsigned chan = 0; chan < 3; chan++) {
        out(4 * i + chan, j) = 50;
      }
    }
  }
}
