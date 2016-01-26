#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <set>

#include <stl_io.hpp>
#include <linear-hash.hpp>

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

void test_bucket_index() {

  size_t M = 3, N = 1000;
  LinearHash<string, string> lh;

  for (size_t i = 0; i < M; ++i) {
    vector<string> keys = generate_keys(N);
    cout << lh.bucket_index(keys[i]) << endl;
  }
}

void test_put() {

  {
    LinearHash<string, string> lh;
    lh.put("key-1", "123");
  }
}

int main() {

  test_bucket_index();
  test_put();

  return 0;
}