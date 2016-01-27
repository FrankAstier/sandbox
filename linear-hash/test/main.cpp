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
  LinearHash_LH<string, string> lh;

  for (size_t i = 0; i < M; ++i) {
    vector<string> keys = generate_keys(N);
    cout << lh.bucket_index(keys[i], 0) << endl;
  }
}

void test() {

  {
    LinearHash_LH<string, string> lh(4, 2);
    for (size_t i = 1; i < 100; ++i) {
      stringstream key; key << "key-" << i;
      stringstream value; value << "value-" << i;
      lh.put(key.str(), value.str());
      lh.print();
      pair<bool,string> v = lh.get(key.str());
      assert(v.first);
      assert(v.second == value.str());
    }
  }
}

int main() {

  test_bucket_index();
  test();

  return 0;
}