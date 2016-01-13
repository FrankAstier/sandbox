/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Frank Astier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef UTILS_FAST_HASH_HPP
#define UTILS_FAST_HASH_HPP

#include <stddef.h>
#include <inttypes.h>
#include <string>

namespace utils {

inline uint64_t fasthash64(const void *buf, size_t len, uint64_t seed =0) {

  auto mix = [](uint64_t h) {
    h ^= h >> 23;
    h *= 0x2127599bf4325c37ULL;
    h ^= h >> 47;
    return h;
  };

  const uint64_t m = 0x880355f21e6d1965ULL;
  const uint64_t *pos = (const uint64_t *) buf;
  const uint64_t *end = pos + (len / 8);
  const unsigned char *pos2;
  uint64_t h = seed ^(len * m);
  uint64_t v;

  while (pos != end) {
    v = *pos++;
    h ^= mix(v);
    h *= m;
  }

  pos2 = (const unsigned char *) pos;
  v = 0;

  switch (len & 7) {
    case 7:
      v ^= (uint64_t) pos2[6] << 48;
    case 6:
      v ^= (uint64_t) pos2[5] << 40;
    case 5:
      v ^= (uint64_t) pos2[4] << 32;
    case 4:
      v ^= (uint64_t) pos2[3] << 24;
    case 3:
      v ^= (uint64_t) pos2[2] << 16;
    case 2:
      v ^= (uint64_t) pos2[1] << 8;
    case 1:
      v ^= (uint64_t) pos2[0];
    default:
      h ^= mix(v);
      h *= m;
  }

  return mix(h);
}

inline uint64_t fasthash64(const std::string& str, uint64_t seed =0) {
  return fasthash64(&str[0], str.size(), seed);
}

} // end namespace utils

#endif //UTILS_FAST_HASH_HPP
