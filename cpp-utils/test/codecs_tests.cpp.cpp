#include <iostream>
#include <sstream>

#include <codecs.hpp>

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

int main() {
  test_vbe();
  return 0;
}