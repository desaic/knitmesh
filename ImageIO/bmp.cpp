#include "bmp.h"

#include "Palette.h"

int encode_bmp_8bit(const unsigned char* imageData, unsigned int width, unsigned int height,
                    std::vector<unsigned char>& out, const Palette_t& palette) {
  unsigned char bmp_file_header[14] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 4, 0, 0};
  unsigned char bmp_info_header[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8, 0};
  unsigned y, padding;

  unsigned char colorPalette[1024];  // @todo(mike) confirm sizes here

  for (int i = 0; i < 256; i++) {
    // @todo(mike) figure out why these are backwards
    colorPalette[4 * i] = palette[i].blue;
    colorPalette[4 * i + 1] = palette[i].green;
    colorPalette[4 * i + 2] = palette[i].red;
    colorPalette[4 * i + 3] = palette[i].alpha;
  }

  // add padding width
  padding = (4 - (width) % 4) % 4;
  unsigned outWidth = width + padding;
  unsigned headerSize = 14 + 40 + 1024;
  const unsigned int size = headerSize + outWidth * height;

  out.resize(size);
  bmp_file_header[2] = (unsigned char)(size);
  bmp_file_header[3] = (unsigned char)(size >> 8);
  bmp_file_header[4] = (unsigned char)(size >> 16);
  bmp_file_header[5] = (unsigned char)(size >> 24);

  bmp_info_header[4] = (unsigned char)(width);
  bmp_info_header[5] = (unsigned char)(width >> 8);
  bmp_info_header[6] = (unsigned char)(width >> 16);
  bmp_info_header[7] = (unsigned char)(width >> 24);
  bmp_info_header[8] = (unsigned char)(height);
  bmp_info_header[9] = (unsigned char)(height >> 8);
  bmp_info_header[10] = (unsigned char)(height >> 16);
  bmp_info_header[11] = (unsigned char)(height >> 24);
  // 256 color table
  bmp_info_header[33] = 1;

  memcpy(out.data(), bmp_file_header, 14);
  memcpy(out.data() + 14, bmp_info_header, 40);
  memcpy(out.data() + 54, colorPalette, 1024);

  for (y = 0; y < height; y++) {
    unsigned srcIdx = (height - y - 1) * width;
    unsigned dstIdx = y * outWidth + headerSize;
    memcpy(out.data() + dstIdx, imageData + srcIdx, width);
  }
  return 0;
}