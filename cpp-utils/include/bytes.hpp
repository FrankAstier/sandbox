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

#ifndef CPP_UTILS_BYTES_HPP
#define CPP_UTILS_BYTES_HPP

#include <assert.h>
#include <string.h>
#include <iostream>
#include <string>
#include <algorithm>

#include "codecs.hpp"

namespace utils {

template <typename S = size_t, typename C = char, bool own_memory =false>
class Bytes {};

/**
 * Just raw bytes, with size - doesn't own memory or copy anything.
 * This is essentially a pair<size_t, char*>, with convenient methods.
 */
template <typename S, typename C>
class Bytes<S,C,false> {
public:
  typedef S size_type;
  typedef C value_type;
  static_assert(sizeof(value_type) == 1, "Bytes requires size 1 value_type");

  typedef const value_type* pointer_type;
  typedef const value_type* const_iterator;

  Bytes() : _size(0), _bytes(nullptr) {}

  Bytes(pointer_type begin, pointer_type end)
      : _size(end - begin), _bytes(begin)
  {
    assert(begin < end);
    check_invariants();
  }

  Bytes(size_t size, pointer_type data) : _size(size), _bytes(data)
  {
    check_invariants();
  }

  Bytes(const std::string& str)
      : _size(str.size()), _bytes(reinterpret_cast<pointer_type>(str.c_str()))
  {
    check_invariants();
  }

  template <typename S1, typename C1, bool O1>
  Bytes(const Bytes<S1,C1,O1>& other)
      : _size(other.size()),
        _bytes(other.data())
  {
    check_invariants();
  }

  Bytes(const Bytes& other)
      : _size(other._size), _bytes(other._bytes)
  { check_invariants(); }

  Bytes(Bytes&& other) : _size(other._size), _bytes(other._bytes) {
    other._size = 0;
    other._bytes = nullptr;
    check_invariants();
  }

  Bytes& operator=(const Bytes& other) {
    if (&other != this) {
      _size = other._size;
      _bytes = other._bytes;
    }
    check_invariants();
    return *this;
  }

  Bytes& operator=(Bytes&& other) {
    if (&other != this) {
      _size = other._size;
      _bytes = other._bytes;
      other._size = 0;
      other._bytes = nullptr;
    }
    check_invariants();
    return *this;
  }

  operator std::string() const { return std::string(reinterpret_cast<const char*>(_bytes), _size); }

  const_iterator begin() const { return _bytes; }
  const_iterator end() const { return _bytes + _size; }
  size_t size() const { return _size; }
  pointer_type data() const { return _bytes; }

private:
  void check_invariants() {
    assert(_size > 0);
    assert(_bytes);
  }
  size_t _size;
  pointer_type _bytes;
};

/**
 * A class that stores bytes with their size.
 *
 * S is a template parameter to allow using a small number of bytes
 * if the length of the keys is known to be small.
 *
 * C can be char, unsigned char, uint8_t... but needs to be 1 byte wide.
 */
template <typename S, typename C>
class Bytes<S,C,true> {
public:
  typedef S size_type;
  typedef C value_type;
  static_assert(sizeof(value_type) == 1, "Bytes requires size 1 value_type");

  typedef const value_type* pointer_type;
  typedef const value_type* const_iterator;

  Bytes()
  : _size(0), _bytes(nullptr) { }

  Bytes(pointer_type data)
  : _size(strlen(data)), _bytes(new value_type[_size]) {
    std::copy(data, data + _size, _bytes);
    check_invariants();
  }

  Bytes(const std::string& data)
  : _size(data.size()), _bytes(new value_type[_size]) {
    std::copy(data.begin(), data.end(), _bytes);
    check_invariants();
  }

  Bytes(const Bytes& other)
  : _size(other._size), _bytes(new value_type[other._size]) {
    std::copy(other._bytes, other._bytes + _size, _bytes);
    check_invariants();
  }

  template <typename S1, typename C1, bool O1>
  Bytes(const Bytes<S1,C1,O1>& other)
      : _size(other.size()),
        _bytes(new value_type[_size])
  {
    std::copy(other.begin(), other.end(), _bytes);
    check_invariants();
  }

  Bytes(Bytes&& other)
  : _size(other._size), _bytes(other._bytes) {
    other._size = 0;
    other._bytes = nullptr;
    check_invariants();
  }

  Bytes& operator=(const Bytes& other) {
    if (&other != this) {
      check_invariants();
      _size = other._size;
      _bytes = new value_type[_size];
      std::copy(other._bytes, other._bytes + _size, _bytes);
    }
    return *this;
  }

  Bytes& operator=(Bytes&& other) {
    if (&other != this) {
      other.check_invariants();
      _size = other._size;
      _bytes = other._bytes;
      other._size = 0;
      other._bytes = nullptr;
    }
    return *this;
  }

  ~Bytes() {
    delete[] _bytes;
    _size = 0;
    _bytes = nullptr;
  }

  operator std::string() const {
    return std::string(begin(), end());
  }

  const_iterator begin() const { return data(); }
  const_iterator end() const { return data() + size(); }
  size_type size() const { return _size; }
  pointer_type data() const { return _bytes; }

  std::ostream& pretty_print(std::ostream& out) const {
    out << "\"";
    for (size_t i = 0; i < size(); ++i)
      out << _bytes[i]; // use +x.bytes[i] to see integer value of bytes
    out << "\"";
    return out;
  }

private:
  void check_invariants() {
    assert(_size > 0);
    assert(_bytes);
  }

  size_type  _size;
  value_type* _bytes;

  template <typename U, typename V>
  friend std::istream& operator>>(std::istream&, Bytes<U,V,true>&);
};

template <typename S, typename C, bool O1, bool O2>
inline bool operator==(const Bytes<S,C,O1>& a, const Bytes<S,C,O2>& b) {
  return a.size() == b.size() && std::equal(a.data(), a.data() + a.size(), b.data());
}

template <typename S, typename C, bool O1, bool O2>
inline bool operator!=(const Bytes<S,C,O1>& a, const Bytes<S,C,O2>& b) {
  return !(a == b);
}

template <typename S, typename C, bool O1>
inline bool operator==(const Bytes<S,C,O1>& a, typename Bytes<S,C,O1>::pointer_type b) {
  return a.size() == strlen(b) && std::equal(a.data(), a.data() + a.size(), b);
}

template <typename S, typename C, bool O1>
inline bool operator!=(const Bytes<S,C,O1>& a, typename Bytes<S,C,O1>::pointer_type b) {
  return !(a == b);
}

template <typename S, typename C, bool O1>
inline bool operator==(typename Bytes<S,C,O1>::pointer_type b, const Bytes<S,C,O1>& a) {
  return a == b;
}

template <typename S, typename C, bool O1>
inline bool operator!=(typename Bytes<S,C,O1>::pointer_type b, const Bytes<S,C,O1>& a) {
  return !(a == b);
}

template <typename S, typename C, bool O1, bool O2>
inline bool operator<(const Bytes<S,C,O1>& a, const Bytes<S,C,O2>& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

template <typename S, typename C, bool O>
inline std::ostream& operator<<(std::ostream& out, const Bytes<S,C,O>& x) {
  auto it = utils::vbe_encode(std::ostreambuf_iterator<char>(out), x.size());
  std::copy(x.data(), x.data() + x.size(), it);
  return out;
}

template <typename S, typename C>
inline std::istream& operator>>(std::istream& in, Bytes<S,C,true>& x) {
  size_t size = 0;
  auto it = utils::vbe_decode(std::istreambuf_iterator<char>(in), size);
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

#endif //CPP_UTILS_BYTES_HPP
