#pragma once

#include <vector>

#include "Palette.h"

int encode_bmp_8bit(const unsigned char* imageData, unsigned int width, unsigned int height,
                    std::vector<unsigned char>& out, const Palette_t& palette);