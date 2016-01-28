#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <set>

#include <stl_io.hpp>
#include "linear-hash-LH.hpp"
#include "linear-hash-LH1.hpp"
#include "linear-hash-PS1.hpp"

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

void test_modulo() {
  for (size_t i = 0; i < 20; ++i) {
    cout << (i % 4) << "\t" << (i % 8) << "\t" << i % 16 << endl;
  }
}

void test() {

  Data data = generate_data(10000);

  {
    cout << "Testing LH" << endl;
    LinearHash_LH<string, string> lh(10, 4);
    for (size_t i = 1; i < data.size(); ++i) {
      lh.put(data[i].first, data[i].second);
      if (i == data.size()-1)
        lh.print();
      pair<bool,string> v = lh.get(data[i].first);
      assert(v.first);
      assert(v.second == data[i].second);
    }
  }

  {
    cout << "Testing LH1" << endl;
    LinearHash_LH1<string, string> lh(10, 4);
    for (size_t i = 1; i < data.size(); ++i) {
      lh.put(data[i].first, data[i].second);
      if (i == data.size()-1)
        lh.print();
      pair<bool,string> v = lh.get(data[i].first);
      assert(v.first);
      assert(v.second == data[i].second);
    }
  }

  {
    cout << "Testing PS1" << endl;
    LinearHash_PS1<string, string> lh(20,2);
    for (size_t i = 1; i < data.size(); ++i) {
      lh.put(data[i].first, data[i].second);
      if (i == data.size()-1)
        lh.print();
      pair<bool,string> v = lh.get(data[i].first);
      assert(v.first);
      assert(v.second == data[i].second);
    }
  }
}

int main() {

  //test_modulo();
  test();

  return 0;
}