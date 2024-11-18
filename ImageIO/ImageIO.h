#pragma once
#include <string>

#include "Array2D.h"

/// <summary>
/// Load a color png into a byte array.
/// </summary>
/// <param name="filename"></param>
/// <param name="slice">output. Width is 3*image width because rgb.</param>
/// <returns>0 on success. -1 file not found. -2 failure to decode
/// file.</returns>
int LoadPngColor(const std::string& filename, Array2D8u& slice);

/// <summary>
/// save a color image as png.
/// Every 3 bytes treated as RGB for a single pixel.
/// Image width is slice.size[0]/3.
/// </summary>
/// <param name="filename"></param>
/// <param name="arr">input. Width is 3*image width because rgb.</param>
/// <returns>0 on success. -1 could not open output file.
/// -2 failure to encode file.</returns>
int SavePngColor(const std::string& filename, const Array2D8u& arr);
int SavePngColor(const std::string& filename, const Array2D4b& arr);
/// Load a gray scale png to a 2D byte array.
/// <returns>0 on success. -1 file not found. -2 failure to decode
/// file.</returns>
int LoadPngGrey(const std::string& filename, Array2D8u& arr);

/// Save a byte array to gray-scale png.
/// <returns>0 on success. -1 fail to write file. -1 failure to encode
/// file.</returns>
int SavePngGrey(const std::string& filename, const Array2D8u& arr);

/// Encode a byte array into bytes for a bmp file in @param out.
void EncodeBmpGray(const unsigned char* img, unsigned int width, unsigned int height,
                   std::vector<unsigned char>& out);