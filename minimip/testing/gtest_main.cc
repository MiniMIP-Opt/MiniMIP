#include "gtest/gtest.h"
#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"

// Unit tests typically link against 'gtest_main' and do not need to implement
// their own main. Here, however, we want some custom control over logging.
// Thus, we provide our own `main()`.
int main(int argc, char **argv) {
  // Initialize logging and testing.
  google::InitGoogleLogging(argv[0]);
  testing::InitGoogleTest(&argc, argv);

  // By default we log to stderr only, and equal or above log level 2 (where
  // 0: INFO, 1: WARNING, 2: ERROR, 3: FATAL).
  absl::SetFlag(&FLAGS_logtostderr, true);
  absl::SetFlag(&FLAGS_minloglevel, 2);
  // Parse extra command line arguments that may override, e.g., minloglevel.
  absl::ParseCommandLine(argc, argv);

  return RUN_ALL_TESTS();
}
