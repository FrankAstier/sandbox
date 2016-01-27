#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <set>

#include <stl_io.hpp>
#include "linear-hash-LH.hpp"
#include "linear-hash-PS.hpp"

using namespace std;
using namespace utils;

using Data = vector<pair<string,string>>;

Data generate_data(size_t N) {

  Data data(N);

  for (size_t i = 0; i < N; ++i) {
    stringstream k; k << "key-" << i;
    stringstream v; v << "value-" << i;
    data[i] = {k.str(), v.str()};
  }

  return data;
}

void test() {

  Data data = generate_data(100);

  {
    LinearHash_LH<string, string> lh(4, 2);
    for (size_t i = 1; i < data.size(); ++i) {
      lh.put(data[i].first, data[i].second);
      lh.print();
      pair<bool,string> v = lh.get(data[i].first);
      assert(v.first);
      assert(v.second == data[i].second);
    }
  }

  {
    LinearHash_PS<string, string> lh(4, 2);
    for (size_t i = 1; i < data.size(); ++i) {
      lh.put(data[i].first, data[i].second);
      lh.print();
      pair<bool,string> v = lh.get(data[i].first);
      assert(v.first);
      assert(v.second == data[i].second);
    }
  }
}

int main() {

  test();

  return 0;
}