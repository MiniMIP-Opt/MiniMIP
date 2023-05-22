// Copyright 2023 the MiniMIP Project
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "gtest/gtest.h"
#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"

// Unit tests typically link against 'gtest_main' and do not need to implement
// their own main. Here, however, we want some custom control over logging.
// Thus, we provide our own `main()`.
int main(int argc, char** argv) {
  // Initialize logging and testing.
  google::InitGoogleLogging(argv[0]);
  testing::InitGoogleTest(&argc, argv);

  // By default, we log to stderr only, and equal or above log level 2 (where
  // 0: INFO, 1: WARNING, 2: ERROR, 3: FATAL).
  absl::SetFlag(&FLAGS_logtostderr, true);
  absl::SetFlag(&FLAGS_minloglevel,
                std::min(absl::GetFlag(FLAGS_minloglevel), 2));
  // Parse extra command line arguments that may override, e.g., minloglevel.
  absl::ParseCommandLine(argc, argv);

  return RUN_ALL_TESTS();
}
