#include <iostream>
#include <sstream>
#include "bytes.hpp"

using namespace std;
using namespace utils;

void test_Bytes() {
  {
    Bytes<> b;
    assert(b.size() == 0);
    assert(b.begin() == b.end());
    assert(b.data() == nullptr);
    assert(string(b) == "");
    stringstream s;
    s << b;
    Bytes<size_t,char,true> b2;
    s >> b2;
    assert(b2.size() == 0);
    assert(b2.begin() == b2.end());
    assert(b2.data() == nullptr);
    assert(string(b2) == "");
    assert(b2 == b);
  }

  {
    Bytes<size_t,char,true> b("abcdef");
    assert(b.size() == 6);
    assert(b == "abcdef");
    stringstream s;
    s << b;
    Bytes<size_t,char,true> b2;
    s >> b2;
    assert(b2 == "abcdef");
    assert(b2 == b);
  }
}

int main() {
  test_Bytes();
  return 0;
}