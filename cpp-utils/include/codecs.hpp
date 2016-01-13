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

#ifndef CPP_UTILS_CODECS_H
#define CPP_UTILS_CODECS_H

#include <assert.h>
#include <stdint.h>
#include <stddef.h>

#include <type_traits>

namespace utils {

/**
 * Encode unsigned integer val on a variable number of bytes at position it.
 */
template <typename T, typename It>
inline It vbe_encode(It it, T val) {
  static_assert(std::is_unsigned<T>::value, "vbe requires unsigned integers");
  while (val > 0x7F) {
    *it++ = (val & 0x7F) | uint8_t(0x80);
    val >>= 7;
  }
  *it++ = val;
  return it;
}

/**
 * Predict the encoded size of unsigned integer val.
 */
template <typename T>
inline size_t vbe_encoded_size(T val) {
  static_assert(std::is_unsigned<T>::value, "vbe requires unsigned integers");
  size_t n = 1;
  while (val > 0x7F) {
    val >>= 7; ++n;
  }
  return n;
}

/**
 * Decode unsigned integer val from position it, and return 1 after end of last
 * encoded byte of val.
 */
template <typename T, typename It>
inline It vbe_decode(It it, T& val) {
  static_assert(std::is_unsigned<T>::value, "vbe requires unsigned integers");
  val = *it & 0x7f;
  for (uint_fast8_t i = 1; *it & 0x80; ++i)
    val |= (*(++it) & 0x7F) << (i * 7);
  return ++it;
}

} // end namespace utils

#endif //CPP_UTILS_CODECS_H
