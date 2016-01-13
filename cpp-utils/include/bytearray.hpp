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

#ifndef CPP_UTILS_BYTEARRAY_HPP
#define CPP_UTILS_BYTEARRAY_HPP

#include <assert.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <string.h>

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

/**
 * A class that stores a pointer to an array of bytes, and its size.
 * This is essentially a pair<size_t, char*>, with convenient methods.
 */
template <typename S = size_t, typename C = char>
class bytearray {
public:
  typedef S size_type;
  typedef C value_type;
  static_assert(sizeof(value_type) == 1, "bytearray requires size 1 value_type");

  typedef const value_type* pointer_type;
  typedef const value_type* const_iterator;

  bytearray()
  : _size(0), _bytes(nullptr) { }

  bytearray(pointer_type data)
  : _size(strlen(data)), _bytes(new value_type[_size]) {
    std::copy(data, data + _size, _bytes);
    assert(_size > 0);
    assert(_bytes != nullptr);
  }

  bytearray(const std::string& data)
  : _size(data.size()), _bytes(new value_type[_size]) {
    std::copy(data.begin(), data.end(), _bytes);
    assert(_size > 0);
    assert(_bytes != nullptr);
  }

  bytearray(const bytearray& other)
  : _size(other._size), _bytes(new value_type[other._size]) {
    std::copy(other._bytes, other._bytes + _size, _bytes);
    assert(_size > 0);
    assert(_bytes != nullptr);
  }

  bytearray(bytearray&& other)
  : _size(other._size), _bytes(other._bytes) {
    other._size = 0;
    other._bytes = nullptr;
    assert(_size > 0);
    assert(_bytes != nullptr);
  }

  bytearray& operator=(const bytearray& other) {
    if (&other != this) {
      assert(other._size > 0);
      assert(other._bytes != nullptr);
      _size = other._size;
      _bytes = new value_type[_size];
      std::copy(other._bytes, other._bytes + _size, _bytes);
    }
    return *this;
  }

  bytearray& operator=(bytearray&& other) {
    if (&other != this) {
      assert(other._size > 0);
      assert(other._bytes != nullptr);
      _size = other._size;
      _bytes = other._bytes;
      other._size = 0;
      other._bytes = nullptr;
    }
    return *this;
  }

  ~bytearray() {
    delete[] _bytes;
    _size = 0;
    _bytes = nullptr;
  }

  operator std::string() const {
    return std::string(begin(), end());
  }

  const_iterator begin() const { return bytes(); }
  const_iterator end() const { return bytes() + size(); }

  size_type size() const { return _size; }
  pointer_type bytes() const { return _bytes; }

  std::ostream& prettyPrint(std::ostream& out) const {
    out << "\"";
    for (size_t i = 0; i < size(); ++i)
      out << _bytes[i]; // use +x.bytes[i] to see integer value of bytes
    out << "\"";
    return out;
  }

private:
  size_type  _size;
  value_type* _bytes;

  template <typename U, typename V>
  friend std::istream& operator>>(std::istream&, bytearray<U,V>&);
};

template <typename S, typename C>
inline bool operator==(const bytearray<S,C>& a, const bytearray<S,C>& b) {
  return a.size() == b.size() && std::equal(a.bytes(), a.bytes() + a.size(), b.bytes());
}

template <typename S, typename C>
inline bool operator==(const bytearray<S,C>& a, typename bytearray<S,C>::pointer_type b) {
  return a.size() == strlen(b) && std::equal(a.bytes(), a.bytes() + a.size(), b);
}

template <typename S, typename C>
inline bool operator==(typename bytearray<S,C>::pointer_type b, const bytearray<S,C>& a) {
  return a == b;
}

template <typename S, typename C>
inline std::ostream& operator<<(std::ostream& out, const bytearray<S,C>& x) {
  auto it = vbe_encode(std::ostreambuf_iterator<char>(out), x.size());
  std::copy(x.bytes(), x.bytes() + x.size(), it);
  return out;
}

template <typename S, typename C>
inline std::istream& operator>>(std::istream& in, bytearray<S,C>& x) {
  size_t size = 0;
  auto it = vbe_decode(std::istreambuf_iterator<char>(in), size);
  if (x._size < size) {
    delete[] x._bytes;
    x._bytes = new char[size];
  }
  x._size = size;
  for (size_t i = 0; i < size; ++i)
    x._bytes[i] = *it++;
  return in;
}

} // end namespace utils

#endif //CPP_UTILS_BYTEARRAY_HPP
