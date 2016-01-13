#include <iostream>
#include <sstream>
#include "bytearray.hpp"

using namespace std;
using namespace utils;

void test_vbe() {
  {
    char buf[16];
    unsigned long x = 123123;
    char* p = vbe_encode(buf, x);
    assert(p - buf == vbe_encoded_size(x));
    unsigned long y;
    assert(vbe_decode(buf, y) == p);
    assert(y == x);
  }

}

void test_bytearray() {
  {
    bytearray<> b;
    assert(b.size() == 0);
    assert(b.begin() == b.end());
    assert(b.bytes() == nullptr);
    assert(string(b) == "");
    stringstream s;
    s << b;
    bytearray<> b2;
    s >> b2;
    assert(b2.size() == 0);
    assert(b2.begin() == b2.end());
    assert(b2.bytes() == nullptr);
    assert(string(b2) == "");
    assert(b2 == b);
  }

  {
    bytearray<> b("abcdef");
    assert(b.size() == 6);
    assert(b == "abcdef");
    stringstream s;
    s << b;
    bytearray<> b2;
    s >> b2;
    assert(b2 == "abcdef");
    assert(b2 == b);
  }
}

int main() {
  test_vbe();
  test_bytearray();
  return 0;
}