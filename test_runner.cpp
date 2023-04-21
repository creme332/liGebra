#define DOCTEST_CONFIG_NO_UNPREFIXED_OPTIONS
#define DOCTEST_CONFIG_IMPLEMENT

#include "tests/doctest.h"

int main(int argc, char** argv) {
  doctest::Context context(argc, argv);

  context.setOption("abort-after",
                    1);  // stop test execution after 1 failed assertion
  context.setOption("no-breaks",
                    true);  // don't break in the debugger when assertions fail

  int test_result = context.run();  // run queries, or run tests unless --no-run

  // --- Comment lines below when running tests locally ---
  if (test_result == 1)
    throw std::runtime_error("Test failed");
  // --- Comment lines above when running tests locally ---

  if (context.shouldExit())  // honor query flags and --exit
    return test_result;

  return test_result;
}