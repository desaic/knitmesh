#pragma once
#ifndef IMAGE_4BIT_H
#define IMAGE_4BIT_H

#include <array>
#include <cinttypes>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

/// Image is automatically padded to multiple of 8 bytes at allocation.
class Image4bit {
 public:
  using Vec2i = std::array<int, 2>;
  struct DualPixel {
    DualPixel() : val(0) {}

    DualPixel(uint8_t p1, uint8_t p2) : val(0) {
      val |= p1 << 4;
      val |= (p2 & 0xf);
    }
    uint8_t pixel1() const { return val >> 4; }
    uint8_t pixel2() const { return val & 0xf; }
    void setPixel1(uint8_t v) { val = (val & 0xf) | (v << 4); }
    void setPixel2(uint8_t v) { val = (val & 0xf0) | (v & 0xf); }
    uint8_t val;
  };

  struct ImageData {
    ImageData() : scanLine(0) {}

    ImageData(int width, int height) { Alloc(width, height); }

    void Alloc(int width, int height) {
      const int bitsPerPixel = 4;
      scanLine = static_cast<size_t>(std::ceil(bitsPerPixel * width / 32.0)) * 4;
      data.resize(scanLine * height);
    }

    size_t GetSize() const { return data.size() * sizeof(data[0]); }

    const char* GetData() const { return reinterpret_cast<const char*>(&data[0]); }
    /// number of bytes in each row including padding
    size_t scanLine;
    std::vector<DualPixel> data;
  };

  Image4bit(const Vec2i& size)
      :  // Make sure the width is a multple of 8 so that the scan line size
         // matches exactly the width of the image, i.e., we prepad the image.
        _size({size[0] + ((8 - size[0] % 8) % 8), size[1]}),
        _data(_size[0], _size[1]) {}

  const Vec2i& GetSize() const { return _size; }

  const ImageData* GetData() const { return &_data; }

  ImageData* GetData() { return &_data; }

  void SetPixel(int x, int y, uint8_t value) {
    size_t index = y * _data.scanLine + x / 2;
    if (x % 2 == 0) {
      _data.data[index].setPixel1(value);
    } else {
      _data.data[index].setPixel2(value);
    }
  }

  uint8_t GetPixel(int x, int y) const {
    size_t index = y * _data.scanLine + x / 2;
    if (x % 2 == 0) {
      return _data.data[index].pixel1();
    } else {
      return _data.data[index].pixel2();
    }
  }

  void ExtendLeft(unsigned x);

 private:
  size_t _GetIndex(int x, int y) const { return (y * _size[0] + x) * sizeof(uint8_t); }

  Vec2i _size;
  ImageData _data;
};

#endif  // SF_IMAGE_H
