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