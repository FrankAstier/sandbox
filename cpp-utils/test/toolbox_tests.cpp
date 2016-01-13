#include <toolbox.hpp>

using namespace std;
using namespace utils;

void test_print_binary() {
  cout << "Testing print_binary" << endl;
  stringstream s;
  print_binary(s, 123);
  assert(s.str() == "01111011");
}

void test_print_hex() {
  cout << "Testing print_hex" << endl;
  stringstream s;
  print_hex(s, 123);
  assert(s.str() == "7B");
}

//--------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  test_print_binary();
  test_print_hex();
  //test_extrema();
  return 0;
}

//--------------------------------------------------------------------------------
