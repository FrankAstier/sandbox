#include <assert.h>
#include <sstream>

#include <stl_io.hpp>

using namespace std;
using namespace utils;

void test_std_container_output() {
  cout << "Testing std container output" << endl;
  {
    vector<int> x{{1,2,3,4}};
    stringstream s; s << x;
    assert(s.str() == "{1,2,3,4}");
  }
  {
    set<int> x = {1,2,3,4};
    stringstream s; s << x;
    assert(s.str() == "{1,2,3,4}");
  }
  {
    map<int,int> m{{ {1,2}, {3,4}, {5,6} }};
    stringstream s; s << m;
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

void test_abbreviate_threshold() {
  cout << "Testing abbreviate_threshold" << endl;
  {
    vector<int> x{{1,2,3,4,5}};
    stringstream s; s << set_abbreviate_threshold(3) << x;
    assert(s.str() == "{1,2,3,..}");
  }
  {
    map<int,int> m{{ {1,2}, {3,4}, {5,6} }};
    stringstream s; s << set_abbreviate_threshold(2) << m;
    assert(s.str() == "{(1,2),(3,4),..}");
  }
}

void test_abbreviate_display_n() {
  cout << "Testing abbreviate_display_n" << endl;
  {
    vector<int> x{{1,2,3,4,5,6,7,8}};
    stringstream s; s << set_abbreviate_threshold(3) << set_abbreviate_display_n() << x;
    assert(s.str() == "{1,2,3,..(5)..}");
  }
  {
    map<int,int> m{{ {1,2}, {3,4}, {5,6}, {7,8} }};
    stringstream s; s << set_abbreviate_threshold(2) << set_abbreviate_display_n() << m;
    assert(s.str() == "{(1,2),(3,4),..(2)..}");
  }
}


int main() {
  test_std_container_output();
  test_abbreviate_threshold();
  test_abbreviate_display_n();
  return 0;
}