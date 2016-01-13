#include <iostream>
#include <sstream>

#include <bytes.hpp>
#include <fast_hash.hpp>

using namespace std;
using namespace utils;

void test_fasthash64() {
  {
    char key[] = "abcdef";
    Bytes<size_t,char,true> x(key);
    Bytes<size_t,char,false> y(key);
    assert(fasthash64(x) == fasthash64(y));
  }
}

int main() {
  test_fasthash64();
  return 0;
}