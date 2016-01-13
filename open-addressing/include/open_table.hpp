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

#ifndef OPEN_ADDRESSING_OPEN_TABLE_HPP
#define OPEN_ADDRESSING_OPEN_TABLE_HPP

#include <stddef.h>
#include <stdint.h>
#include <vector>
#include <iostream>

#include <bytes.hpp>
#include <fast_hash.hpp>

/**
 * This is a specialized open-addressing hash table that use bytearray keys
 * and store uint64_t, which could represent offsets into another data structure
 * that actually stores values corresponding to these offsets.
 *
 * Another feature is that this table does not dynamically resize: it needs to be
 * sized once at construction time, and cannot accommodate more entries when it is
 * full.
 */
class OpenTable {
public:
  OpenTable(size_t N)
      : _table(N)
  {}

  /**
   * Uses linear probing to find a spot in the hash table. The key is stored
   * alongside the value, so that we can check if we have the correct key
   * when we find.
   */
  template <typename S, typename C, bool O>
  void insert(const utils::Bytes<S,C,O>& key, uint64_t value) {
    size_t h = utils::fasthash64(key) % _table.size();
    size_t c = 0;
    while (_table[h].first != 0 && c < _table.size()) {
      h = (h+1) % _table.size();
      ++c;
    }
    if (c == _table.size()) {
      exit(-9); // table full
    }
    if (_table[h].first == 0) {
      _table[h] = {value, key};
    }
  }

  /**
   * Uses same linear probing strategy as insert. The first element of the pair
   * is true only if the key has a corresponding value in the table. It is false
   * if the table does not contain the given key.
   */
  template <typename S, typename C, bool O>
  std::pair<bool,uint64_t> find(const utils::Bytes<S,C,O>& key) {
    size_t h = utils::fasthash64(key) % _table.size();
    size_t c = 0;
    while (_table[h].second != key && c < _table.size()) {
      h = (h+1) % _table.size();
      ++c;
    }
    if (c == _table.size()) {
      return {false, -1};
    }
    if (_table[h].second == key)
      return {true, _table[h].first};
    else
      return {false, 0};
  }

private:
  std::vector<std::pair<uint64_t, utils::Bytes<size_t,char,false>>> _table;
};

#endif //OPEN_ADDRESSING_OPEN_TABLE_HPP
