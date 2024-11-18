#pragma once

#include <array>

struct Color_t {
  uint8_t red = 0u;
  uint8_t green = 0u;
  uint8_t blue = 0u;
  uint8_t alpha = 0u;  // unused
};
using Palette_t = std::array<Color_t, 256>;
