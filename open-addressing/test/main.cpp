#include <iostream>
#include <sstream>

#include <stl_io.hpp>
#include <open_table.hpp>

using namespace std;
using namespace utils;

vector<Bytes<size_t,char,true>> generate_keys(size_t N) {

  vector<Bytes<size_t,char,true>> keys(N);

  for (size_t i = 0; i < N; ++i) {
    stringstream s;
    s << "key-" << i;
    keys[i] = s.str();
  }

  return keys;
}

void test_open_addressing() {

  size_t N = 20;
  vector<Bytes<size_t,char,true>> keys = generate_keys(N);
  OpenTable table(N);
  for (size_t i = 0; i < N; ++i)
    table.insert(keys[i], i);

  for (size_t i = 0; i < N; ++i) {
    cout << table.find(keys[i]) << endl;
  }
}

int main() {

  test_open_addressing();
  return 0;
}