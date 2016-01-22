#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <set>

#include <fast_hash.hpp>
#include <stl_io.hpp>

using namespace std;
using namespace utils;

using Bucket = list<string>;

vector<string> generate_keys(size_t N) {

  vector<string> keys(N);

  for (size_t i = 0; i < N; ++i) {
    stringstream s;
    s << "key-" << i;
    keys[i] = s.str();
  }

  return keys;
}

/**
 * i: splitting round number
 * p: bucket to split next
 */
uint64_t bucket_index(const string& key, size_t i, size_t p) {
  uint64_t h = fasthash64(key, 0);
  return h;
}


void test_linear_hash() {

  size_t M = 1000, N = 1000;

  for (size_t i = 0; i < M; ++i) {
    vector<string> keys = generate_keys(N);

  }
}

int main() {

  test_linear_hash();
  return 0;
}