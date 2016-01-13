#include <assert.h>
#include <sstream>

#include <stl_io.hpp>

using namespace std;
using namespace utils;

void test_std_container_output() {
  cout << "Testing std container output" << endl;
  {
    vector<int> x{{1,2,3,4}};
    stringstream s;
    s << x;
    assert(s.str() == "{1,2,3,4}");
  }
  {
    set<int> x = {1,2,3,4};
    stringstream s;
    s << x;
    assert(s.str() == "{1,2,3,4}");
  }
  {
    map<int,int> m{{ {1,2}, {3,4}, {5,6} }};
    stringstream s;
    s << m;
    assert(s.str() == "{(1,2),(3,4),(5,6)}");
  }
  {
    vector<list<string>> x;
    stringstream s0; s0 << x;
    assert(s0.str() == "{}");
    x.push_back(list<string>());
    stringstream s;
    s << x;
    assert(s.str() == "{{}}");
  }
}

int main() {
  test_std_container_output();
  return 0;
}