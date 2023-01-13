#define DOCTEST_CONFIG_NO_UNPREFIXED_OPTIONS
#define DOCTEST_CONFIG_IMPLEMENT

#include "tests/doctest.h"

int main(int argc, char** argv) {
  doctest::Context context(argc, argv);

  int test_result = context.run();  // run queries, or run tests unless --no-run

  if (context.shouldExit())  // honor query flags and --exit
    return test_result;

  std::cout << "hello world";
}