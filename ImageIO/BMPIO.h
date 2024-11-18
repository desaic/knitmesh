#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

class SfBMPIO {
 public:
#pragma pack(push, 1)
  struct Header {
    Header() : type{'B', 'M'}, fileSize(0), reserved1{0, 0}, reserved2{0, 0}, imageOffset(0) {}

    char type[2];
    uint32_t fileSize;
    char reserved1[2];
    char reserved2[2];
    uint32_t imageOffset;
  };
#pragma pack(pop)

#pragma pack(push, 1)
  struct DIBHeader {
    DIBHeader() { std::memset(this, 0, sizeof(*this)); }

    DIBHeader(int _width, int _height)
        : headerSize(sizeof(DIBHeader)),
          width(_width),
          height(_height),
          colorPlanes(1),  // always one bit color plane
          bitsPerPixel(4),
          compression(0),  // none
          imageSize(0),
          xResolution(31000),  // not used
          yResolution(16000),
          numColors(0),          // defaults to 2^bits
          numImportantColors(0)  // every color is important
    {}

    uint32_t headerSize;
    int32_t width;
    int32_t height;
    uint16_t colorPlanes;
    uint16_t bitsPerPixel;
    uint32_t compression;
    uint32_t imageSize;
    int32_t xResolution;  // pixel per meter
    int32_t yResolution;  // pixel per meter
    uint32_t numColors;
    uint32_t numImportantColors;
  };
#pragma pack(pop)
};
