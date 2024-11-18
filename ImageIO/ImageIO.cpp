#include "ImageIO.h"

#include <filesystem>

#include "lodepng.h"

int LoadPngColor(const std::string& filename, Array2D8u& slice) {
  std::filesystem::path p(filename);
  if (!(std::filesystem::exists(p) && std::filesystem::is_regular_file(p))) {
    return -1;
  }
  std::vector<unsigned char> buf;
  lodepng::State state;
  // input color type
  state.info_raw.colortype = LCT_RGB;
  state.info_raw.bitdepth = 8;
  // output color type
  state.info_png.color.colortype = LCT_RGB;
  state.info_png.color.bitdepth = 8;
  unsigned error = lodepng::load_file(buf, filename.c_str());
  if (error) {
    return -1;
  }
  unsigned w, h;
  slice.GetData().clear();
  error = lodepng::decode(slice.GetData(), w, h, state, buf);
  if (error) {
    return -2;
  }
  slice.SetSize(w * 3, h);
  return 0;
}

int SavePngColor(const std::string& filename, const Array2D8u& arr) {
  std::vector<unsigned char> buf;
  lodepng::State state;
  state.info_raw.colortype = LCT_RGB;
  state.info_raw.bitdepth = 8;
  state.info_png.color.colortype = LCT_RGB;
  state.info_png.color.bitdepth = 8;
  unsigned error =
      lodepng::encode(buf, arr.GetData(), arr.GetSize()[0] / 3, arr.GetSize()[1], state);
  if (error) {
    return -2;
  }
  error = lodepng::save_file(buf, filename.c_str());
  if (error) {
    return -1;
  }
  return 0;
}

int SavePngColor(const std::string& filename, const Array2D4b& arr) {
  std::vector<unsigned char> buf;
  lodepng::State state;
  state.info_raw.colortype = LCT_RGBA;
  state.info_raw.bitdepth = 8;
  state.info_png.color.colortype = LCT_RGB;
  state.info_png.color.bitdepth = 8;
  unsigned error = lodepng::encode(buf, (unsigned char*)arr.DataPtr(), arr.GetSize()[0],
                                   arr.GetSize()[1], state);
  if (error) {
    return -2;
  }
  error = lodepng::save_file(buf, filename.c_str());
  if (error) {
    return -1;
  }
  return 0;
}

int LoadPngGrey(const std::string& filename, Array2D8u& arr) {
  std::vector<unsigned char> buf;
  lodepng::State state;
  state.info_raw.colortype = LCT_GREY;
  state.info_raw.bitdepth = 8;
  state.info_png.color.colortype = LCT_GREY;
  state.info_png.color.bitdepth = 8;
  unsigned error = lodepng::load_file(buf, filename.c_str());
  if (error) {
    return -1;
  }
  unsigned w, h;
  arr.GetData().clear();
  error = lodepng::decode(arr.GetData(), w, h, state, buf);
  if (error) {
    return -2;
  }
  arr.SetSize(w, h);
  return 0;
}

int SavePngGrey(const std::string& filename, const Array2D8u& arr) {
  std::vector<unsigned char> buf;
  lodepng::State state;
  state.info_raw.colortype = LCT_GREY;
  state.info_raw.bitdepth = 8;
  state.info_png.color.colortype = LCT_GREY;
  state.info_png.color.bitdepth = 8;
  unsigned error = lodepng::encode(buf, arr.GetData(), arr.GetSize()[0], arr.GetSize()[1], state);
  if (error) {
    return -2;
  }
  error = lodepng::save_file(buf, filename.c_str());
  if (error) {
    return -1;
  }
  return 0;
}

void EncodeBmpGray(const unsigned char* img, unsigned int width, unsigned int height,
                   std::vector<unsigned char>& out) {
  unsigned char bmp_file_header[14] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 4, 0, 0};
  unsigned char bmp_info_header[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8, 0};
  unsigned y, padding;

  unsigned char colorPalette[1024];

  for (unsigned i = 0; i < 256; i++) {
    colorPalette[i * 4] = i;
    colorPalette[i * 4 + 1] = i;
    colorPalette[i * 4 + 2] = i;
    colorPalette[i * 4 + 3] = 0;
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
    memcpy(out.data() + dstIdx, img + srcIdx, width);
  }
}
